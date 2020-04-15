
use parse_args;
use std::str;
use std::process::{Command, Stdio, exit};
use std::ascii::AsciiExt;
use std::io::{Write, BufReader, BufRead, stderr};

const USAGE: &'static str = "
Usage:
  variant call [options] <genome.fa> <bam_files>...

Options:
  --region=REGION   Genomic region to analyze [default: entire genome]
  --alt-reads=N     Minimum alt read number [default: 5]
  --alt-frac=N      Minimum alt allele fraction [default: 0.1]
  --min-mapq=N      Only count reads with MAPQ > N [default: 0]
";

struct Allele {
	signature: Vec<u8>,
	count: Vec<usize>
}

pub fn main() {
	let args = parse_args(USAGE);
	let genome_path = args.get_str("<genome.fa>");
	let bam_paths = args.get_vec("<bam_files>");
	let alt_reads: usize = args.get_str("--alt-reads").parse().unwrap();
	let alt_frac: f64 = args.get_str("--alt-frac").parse().unwrap();
	let region = args.get_str("--region");
	let min_mapq: usize = args.get_str("--min-mapq").parse().unwrap();
	let S = bam_paths.len();

	print!("CHROM\tPOSITION\tREF\tALT\tNOTES");
	for bam_path in &bam_paths {
		print!("\t{}", bam_path.replace(".bam", ""));
	}
	println!();

	let mut cmd = Command::new("samtools");
	cmd.arg("mpileup");
	cmd.arg("-d"); cmd.arg("1000000");
	cmd.arg("-AxRsB");
	cmd.arg("-q"); cmd.arg(min_mapq.to_string());
	if region.ends_with(".bed") {
		cmd.arg("-l"); cmd.arg(region);
	} else if region != "" {
		cmd.arg("-r"); cmd.arg(region);
	}
	cmd.arg("-f"); cmd.arg(genome_path);
	for bam_path in &bam_paths { cmd.arg(bam_path); }
	let mpileup = cmd.stdout(Stdio::piped()).spawn()
		.expect("Could not start 'samtools mpileup' process.");

	let mpileup_out = BufReader::new(mpileup.stdout.unwrap());
	for l in mpileup_out.lines() {
		let line = l.unwrap();

		if line.is_empty() { continue; }
		let mut cols = line.split('\t');

		let chr = cols.next().unwrap();
		let pos = cols.next().unwrap();
		let ref_allele = cols.next().unwrap();
		assert!(ref_allele.len() == 1);
		let ref_allele = ref_allele.bytes().next().unwrap().to_ascii_uppercase();
		if ref_allele == 'N' as u8 { continue; }

		let mut alleles: Vec<Allele> = Vec::new();
		let mut total_reads = vec!(0; S);
		for s in 0..S {
			let mpileup_read_count: usize = cols.next().unwrap().parse().unwrap();
			let pileup = cols.next().unwrap().as_bytes();
			let baseq = cols.next().unwrap();
			let mapq = cols.next().unwrap();
			if mpileup_read_count == 0 { continue; }

			// Convert ASCII format pileup to a vector of observed alleles.
			let allele_pileup = parse_pileup(&pileup);

			let mut total = 0;   // Count total reads at the site
			'outer: for read in allele_pileup {
				// Don't count indels that contain unknown bases at the site.
				if read.contains(&('N' as u8)) { continue; }

				total += 1;
				for allele in &mut alleles {
					if allele.signature == read {
						allele.count[s] += 1; continue 'outer;
					}
				}
				let mut new = Allele { signature: read, count: vec!(0; S) };
				new.count[s] += 1;
				alleles.push(new);
			}
			total_reads[s] = total;
		}

		for mut allele in alleles {
			// Skip reference alleles '.' and deleted bases '*' unless 
			// we are in "report all" mode (i.e. --min-alt-reads=0).
			if allele.signature.len() == 1 && alt_reads > 0 {
				if allele.signature[0] == '.' as u8 { continue; }
				if allele.signature[0] == '*' as u8 { continue; }
			}

			let mut starred = vec!(false; S);
			for s in 0..S {
				starred[s] = allele.count[s] >= alt_reads &&
					allele.count[s] as f64 / total_reads[s] as f64 >= alt_frac
			}
			if starred.contains(&true) == false { continue; }

			print!("{}\t{}\t", chr, pos);

			// Replace '.' in alternate allele with reference base
			if allele.signature[0] == '.' as u8 {
				allele.signature[0] = ref_allele;
			}

			if allele.signature.len() >= 2 {
				// Convert indels to VCF4 format
				if allele.signature[1] == '+' as u8 {
					print!("{}\t{}{}", ref_allele as char, allele.signature[0] as char, str::from_utf8(&allele.signature[2..]).unwrap());
				} else if allele.signature[1] == '-' as u8 {
					print!("{}{}\t{}", ref_allele as char, str::from_utf8(&allele.signature[2..]).unwrap(), allele.signature[0] as char);
				} else {
					eprintln!("Invalid indel: {}", str::from_utf8(&allele.signature).unwrap());
					exit(-1);
				}
			} else {
				print!("{}\t{}", ref_allele as char, str::from_utf8(&allele.signature).unwrap());
			}
			print!("\t");

			for s in 0..S {
				print!("\t{}:{}", allele.count[s], total_reads[s]);
				if starred[s] { print!(":*"); }
			}
			println!();
		}
	}
}

fn parse_pileup(pileup: &[u8]) -> Vec<Vec<u8>> {
	let mut alleles: Vec<Vec<u8>> = Vec::new();
	let mut i = 0;
	while i < pileup.len() {
		if pileup[i] == b'^' {
			i += 2;
		} else if pileup[i] == b'$' || pileup[i] == b'>' || pileup[i] == b'<' {
			i += 1;     // Ignore
		} else if i + 1 < pileup.len() && (pileup[i + 1] == b'+' || pileup[i + 1] == b'-') {
			// Example of insertion in mpileup output: A+4GGTA
			// Example of deletion in mpileup output: A-3NNN
			let mut indel_len: usize = 0;
			let mut l = i + 2;    // First byte of the length mark
			while pileup[l] >= b'0' && pileup[l] <= b'9' {
				indel_len *= 10;
				indel_len += (pileup[l] - b'0') as usize;
				l += 1;
			}
			let mut indel_seq = vec!(b' '; 2 + indel_len);
			indel_seq[0] = omit_strand(pileup[i]);   // Base before the indel
			indel_seq[1] = pileup[i + 1];            // Plus or minus sign
			for k in 0..indel_len {
				indel_seq[2+k] = omit_strand(pileup[l+k]);
			}
			alleles.push(indel_seq);
			i = l + indel_len;
		} else if pileup[i] == b'N' {   // N must be checked after indels
			i += 1;     // Ignore
		} else if (pileup[i] >= b'A' && pileup[i] <= b'Z') || (pileup[i] >= b'a' && pileup[i] <= b'z') || pileup[i] == b'.' || pileup[i] == b',' || pileup[i] == b'*' {
			alleles.push(vec![omit_strand(pileup[i])]);
			i += 1;
		} else {
			eprintln!("Invalid pileup symbol #{} detected. Pileup looked like this:\n{}", pileup[i], str::from_utf8(pileup).unwrap());
			exit(-1);
		}
	}
	return alleles;
}

fn omit_strand(allele: u8) -> u8 {
	if allele == ',' as u8 { '.' as u8 } else { allele.to_ascii_uppercase() }
}