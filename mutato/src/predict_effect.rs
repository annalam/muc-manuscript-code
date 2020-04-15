use parse_args;

const USAGE: &str = "
Usage:
  mutato predict effect [options] <vcf_file>

Options:
  --genome=GENOME     Genomic region to include in analysis [default: all]
";

struct VcfReader {
	samples: Vec<String>
}

struct VcfRow {
	chromosome: String,
	position: u32,
	ref_allele: String,
	alt_allele: String,
	alt_reads: Vec<u32>,
	total_reads: Vec<u32>
}

pub fn main() {
	let args = parse_args( USAGE);
	let vcf_path = args.get_str("<vcf_file>");
	let annovar_genome = args.get_str("--genome");
	let genome_version = annovar_genome.basename();

	let reformatted_path = format!("annovar-{}.avinput", "46g7");

	let mut vcf_file = VcfReader(vcf_path);
	for variant in vcf_file {
		if variant.ref_allele.len() == 1 && variant.alt_allele.len() > 1 && variant.ref_allele[0] == variant.alt_allele[0] {
			// Simple insertion
		} else if variant.ref_allele.len() > 1 && variant.alt_allele.len() == 1 && variant.ref_allele[0] == variant.alt_allele[0] {
			// Simple deletion
		}

		writeln!(reformatted_file, "");
	}
	reformatted_file.close();

	println!("CHROM\tPOSITION\tREF\tALT\tGENE\tEFFECT\tNOTES\t");

}

fn translate_annovar_effect() {
	// for f in split(func, ';')
	// 	if !in(f, annovar_valid_funcs)
	// 		warn("Unrecognized variant effect '$f'.")
	// 	end
	// end

	// effects = []
	// if contains(func, "splicing"); push!(effects, "Splice site"); end
	// if contains(func, "exonic")
	// 	details = unique(map(m -> m.captures[1],
	// 		eachmatch(r":(p\..+?)(,|$)", aa_change)))
	// 	if exonic_func == "synonymous SNV"
	// 		push!(effects, "Synonymous $(join(details, ", "))")
	// 	elseif exonic_func == "nonsynonymous SNV"
	// 		push!(effects, "Missense $(join(details, ", "))")
	// 	elseif exonic_func == "stopgain"
	// 		push!(effects, "Stopgain $(join(details, ", "))")
	// 	elseif exonic_func == "stoploss"
	// 		push!(effects, "Stoploss $(join(details, ", "))")
	// 	elseif contains(exonic_func, "nonframeshift")
	// 		push!(effects, "Non-frameshift indel $(join(details, ", "))")
	// 	elseif contains(exonic_func, "frameshift")
	// 		push!(effects, "Frameshift $(join(details, ", "))")
	// 	elseif contains(func, "ncRNA_exonic")
	// 		push!(effects, "Exonic (ncRNA)")
	// 	elseif exonic_func == "unknown"
	// 		push!(effects, "Exonic (unknown)")
	// 	else
	// 		error("Unrecognized effect: $([func, exonic_func])")
	// 	end
	// end
	// if contains(func, "UTR3"); push!(effects, "3'-UTR"); end
	// if contains(func, "UTR5"); push!(effects, "5'-UTR"); end
	// if contains(func, "upstream"); push!(effects, "Upstream"); end
	// if contains(func, "downstream"); push!(effects, "Downstream"); end
	// if contains(func, "intronic"); push!(effects, "Intronic"); end
	// if func == "intergenic"; push!(effects, "Intergenic"); gene = ""; end

	// return (gene, join(effects, ". "))
}