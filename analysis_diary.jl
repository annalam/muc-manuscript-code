
# ALIGNMENT STATISTICS
cd ~/datasets/bladder_cfdna/panel_v1/alignments
mkdir ../statistics
echo *.bam | parallel -n10 'sam statistics --on-target=../baits_hg38.bed $x > ../statistics/${x/.bam/.general} && sam coverage histogram --regions=../baits_hg38.bed $x > ../statistics/${x/.bam/.coverage}'

cd ~/datasets/bladder_cfdna/panel_v1/statistics
using QualityMetrics
samples = [replace(s, ".general", "") for s in readdir()
	if endswith(s, ".general")];
plot_quality_metrics(samples)







# IDENTIFY SOMATIC MUTATIONS AND GERMLINE VARIANTS
cd ~/datasets/bladder_cfdna/alignments
mkdir ../mutations
echo X `seq 22` Y | parallel -n10 'mutato call2 --region=chr${x} --alt-reads=5 --alt-frac=0.02 ~/homo_sapiens/hg38.fa *.bam > ../mutations/chr${x}.vcf'

cd ~/datasets/bladder_cfdna/mutations
cat chr1.vcf <(cat chr{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf | grep -v CHROM) > variants.vcf

# Generate list of tumor-normal pairs
cd ~/datasets/bladder_cfdna
samples = [replace(p, ".bam", "") for p in readdir("alignments") if endswith(p, ".bam")];
germline = [s for s in samples if r"-(WBC|[Gg]ermline|[Bb]enign)" in s];
patients = unique(replace.(germline, r"-(WBC|[Gg]ermline|[Bb]enign)", ""));
for s in samples
	if !any(startswith(s, p) for p in patients); println("Orphan: $s"); end
end
out = open("tumor_normal_pairs.txt", "w");
write(out, "TEST\tREF\tALT_FRAC\n");
for patient in patients
	patient_samples = filter(s -> startswith(s, "$(patient)-"), samples)
	normals = filter(s -> s in germline, patient_samples)
	for s in patient_samples
		if s in germline; continue; end
		alt_frac = 0.01
		if !("-cfDNA" in s); alt_frac = 0.08; end
		if isempty(normals)
			write(out, "$(s)\t\t$(alt_frac)\n")
		else
			write(out, "$(s)\t$(normals[1])\t$(alt_frac)\n")
		end
	end
end
close(out)

# Somatic mutations (TERT promoter mutations use special thresholds)
cd ~/datasets/bladder_cfdna/mutations
variant nearby indels variants.vcf | variant somatic --alt-reads=10 --test-ref-ratio=3 --test-bg-ratio=25 --ref-reads=10 --min-sidedness=15 --min-mapq=10 - ../tumor_normal_pairs.txt | variant predict effect - | grep -v TERT > somatic.tmp
variant inside variants.vcf chr5:1295113-1295135 | variant somatic --alt-reads=3 --alt-frac=0 --test-ref-ratio=0 --test-bg-ratio=0 --ref-reads=0 - ../tumor_normal_pairs.txt | variant predict effect - | tail -n +2 >> somatic.tmp    # TERT mutations
cat <(variant protein altering somatic.tmp) <(variant protein altering --invert somatic.tmp | grep -v INDEL | tail -n +2) | sort -k1,1V -k2,2n | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | variant annotate - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jls 0.005 | variant discard if frequency above - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls 0.005 > somatic.vcf
# Correct chr11:108268612:G>C ATM mutation effect to "Splice site"

# Germline variants
variant germline --alt-reads=8 --alt-frac=0.15 --bg-ratio=20 variants.vcf WBC | variant predict effect - | variant protein altering - | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jls 0.005 | variant annotate - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls | variant discard if frequency above - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls 0.005 | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | variant annotate - ~/homo_sapiens/clinvar_2016-07-07_hg38.jls > germline.vcf

# Heterozygous germline SNPs
variant heterozygous snps --min-depth=50 variants.vcf WBC | variant discard indels - | variant annotate - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | egrep 'CHROM|KAVIAR|ExAC' | variant discard blacklisted - ~/hetz_snp_blacklist.tsv > hetz_snps.vcf
variant group by snp profile hetz_snps.vcf > ../snp_profiles.txt





# SUPPLEMENTARY TABLE: SOMATIC MUTATIONS
cd ~/datasets/bladder_cfdna/mutations
using Variant;
vcf = read_vcf(`variant discard blacklisted somatic.vcf blacklist.tsv`);

out = open("/home/annalam/supplementary_table.tsv", "w");
function print_allele(row::Int, col::Int)
	@printf(out, "\t%.1f%% (%d)%s", vcf.alt[row, col] / vcf.total[row, col] * 100, vcf.total[row, col], vcf.star[row, col] ? " *" : "")
end

d = readtsv("../samples_in_muc_cohort.txt"); d = d[2:end, :];
all_normals = unique(filter(s -> s != "", d[:, 2]));
background = indexin(all_normals, vcf.sample);

for normal in all_normals
	tumors = [d[r, 1] for r in 1:size(d, 1) if d[r, 2] == normal]

	write(out, "CHROM\tPOSITION\tREF\tALT\tGENE\tEFFECT")
	for tumor in tumors; write(out, "\t$tumor"); end
	write(out, "\t$normal\tALL WBC\n")

	# Easier to convert directly to sample indexes here
	tumors = indexin(tumors, vcf.sample)
	normal = findone(vcf.sample, normal)

	for r in 1:size(vcf.star, 1)
		if !any(vcf.star[r, tumors]); continue; end
		#if !is_protein_altering(vcf.effect[r]); continue; end  # TODO: All muts
		bg_frac = sum(vcf.alt[r, background]) / sum(vcf.total[r, background])
		@printf(out, "%s\t%d\t%s\t%s\t%s\t%s", vcf.chromosome[r],
				vcf.position[r], vcf.ref_allele[r], vcf.alt_allele[r],
				vcf.gene[r], vcf.effect[r])
		for tumor in tumors; print_allele(r, tumor); end
		print_allele(r, normal)
		@printf(out, "\t%.1f%%\n", bg_frac * 100)
	end
	write(out, '\n')
end
close(out)







# ANALYZE COPY NUMBER ALTERATIONS (PANEL VERSION 1)
cd ~/datasets/bladder_cfdna/panel_v1/alignments
echo *.bam | parallel -n8 '[ ! -e ../coverage/${x/.bam/.tsv} ] && sam count $x ../baits_hg38.bed > ../coverage/${x/.bam/.tsv}'

cd ~/datasets/bladder_cfdna/panel_v1/coverage
copynum call targeted --controls=../cna_controls.txt --gc-fractions=../baits_hg38.gc --report-dir=~/tmp ../baits_hg38.bed ../tumor_normal_pairs.txt > ../gene_cna.tsv



# ANALYZE COPY NUMBER ALTERATIONS (PANEL VERSION 2)
cd ~/datasets/bladder_cfdna/panel_v2/alignments
echo *.bam | parallel -n8 '[ ! -e ../coverage/${x/.bam/.tsv} ] && sam count $x ../baits_hg38.bed > ../coverage/${x/.bam/.tsv}'

# Analyze cfDNA samples (panel version 2)
cd ~/datasets/bladder_cfdna/panel_v2/coverage
copynum call targeted --controls=../cna_controls.txt --gc-fractions=../baits_hg38.gc --report-dir=~/tmp ../baits_hg38.bed ../tumor_normal_pairs.txt > ../gene_cna.tsv








# DETECT CHROMOSOMAL REARRANGEMENTS USING BREAKFAST
cd ~/datasets/bladder_cfdna/alignments
echo *.bam | parallel -n10 'breakfast detect --max-frag-len=1000 --anchor-len=30 $x ~/tools/bowtie-indexes/homo_sapiens/hg38 > ../rearrangements/original/${x/.bam/.sv}'

cd ~/datasets/bladder_cfdna/rearrangements/original
breakfast blacklist *-WBC.sv *enign*.sv > ../blacklist.txt
echo *.sv | parallel -n10 'breakfast filter --merge-duplicates --min-reads=5 --blacklist=../blacklist.txt $x | breakfast annotate - ~/homo_sapiens/ensembl_84/genes.bed > ../annotated/${x}'

cd ~/datasets/bladder_cfdna
breakfast matrix --threads=20 <(cat rearrangements/annotated/*.sv) alignments/*.bam > rearrangements/counts.tsv






# SUPPLEMENTARY TABLE: SOMATIC REARRANGEMENTS
cd ~/datasets/bladder_cfdna/rearrangements/annotated

struct Rearrangement
	sample::String
	patient::String
	chrom_left::String
	strand_left::Char
	position_left::Int
	features_left::String
	chrom_right::String
	strand_right::Char
	position_right::Int
	features_right::String
	best_read::String
end

rearrangements = Vector{Rearrangement}();
for sv_file in [f for f in readdir() if endswith(f, ".sv")]
	sample = replace(sv_file, ".sv", "")
	if "percent" in sample || "cell-line" in sample; continue; end
	m = match(r"^(BC-\d\d\d|BC-T2|d\d+|E-\d+|ER\d+|Ghent-\d+|GU-\d\d-\d+|PS12|T-\d+)",
		sample)
	if m == nothing; error("Cannot infer patient name for $(sample)."); end
	patient = m.captures[1]

	d = readtsv(sv_file)[2:end, :]
	for r in 1:size(d, 1)
		reads = split(d[r, 9], ';')
		min_flank = map(reads) do read
			left, right = split(read, '|')
			min(sum(isuppercase, left), sum(isuppercase, right))
		end
		best_read = reads[argmax(min_flank)]
		push!(rearrangements, Rearrangement(sample, patient,
			d[r, 1], d[r, 2][1], d[r, 3], d[r, 4],
			d[r, 5], d[r, 6][1], d[r, 7], d[r, 8], best_read))
	end
end

matrix = readtsv("../counts.tsv");
samples = matrix[1, 12:end];
counts = matrix[2:end, 12:end];
info = matrix[2:end, 1:7];

wbc = r"WBC|enign" .∈ samples;

out = open("/home/annalam/rearrangements.tsv", "w");
@printf(out, "Patient\tCoordinate 1\tCoordinate 2\tFeatures 1\tFeatures 2\tJunction sequence\tReads\tNotes\n");
for r in rearrangements
	@printf(out, "%s", r.sample)
	@printf(out, "\t%s:%d (%s)", r.chrom_left, r.position_left, r.strand_left)
	@printf(out, "\t%s:%d (%s)", r.chrom_right, r.position_right, r.strand_right)
	@printf(out, "\t%s\t%s", r.features_left, r.features_right)
	@printf(out, "\t%s", r.best_read)

	rows = findall(k -> r.chrom_left == info[k, 1] &&
		r.position_left == info[k, 3] && r.chrom_right == info[k, 5] &&
		r.position_right == info[k, 7], 1:size(info, 1))
	if length(rows) != 1; print(out, "\t\t\t\t\n"); continue; end
	row = rows[1]

	# Raw number of supporting reads
	col = findone(s -> s == r.sample, samples)
	if col == nothing; error("Could not find sample $(r.sample)."); end
	@printf(out, "\t%d", counts[row, col])

	# Add notes if other samples also show evidence for the rearrangement
	notes = fill("", 0)
	for s in findall(counts[row, :] .> 0)
		if startswith(samples[s], r.patient); continue; end
		push!(notes, @sprintf("%s: %d reads", samples[s], counts[row, s]))
	end
	if length(notes) >= 5; insert!(notes, 1, "SUSPICIOUS"); end
	@printf(out, "\t%s", join(notes, ", "))
	@printf(out, "\n")
end
close(out)








# ESTIMATE CTDNA FRACTION
cd ~/datasets/bladder_cfdna/mutations
clonality targeted --report-dir=~/ <(variant discard blacklisted somatic.vcf blacklist.tsv) ../tumor_normal_pairs.txt > ../ctdna_fractions.tsv






# SUPPLEMENTARY TABLE: SOMATIC MUTATION BURDEN
cd ~/datasets/bladder_cfdna
variant discard blacklisted mutations/somatic.vcf mutations/blacklist.tsv | variant mutation rate --detailed -alt-reads=10 - coverage_histograms >> mutations/mutation_rate.tsv











##########################
# WHOLE EXOME SEQUENCING #
##########################

# ALIGNMENT STATISTICS
cd ~/datasets/bladder_cfdna/wxs/alignments
echo *.bam | parallel -n8 '[ ! -e ../statistics/${x/.bam/.general} ] && sam statistics $x > ../statistics/${x/.bam/.general} && sam coverage histogram --regions=/home/annalam/homo_sapiens/ensembl_84/cds.bed $x > ../statistics/${x/.bam/.coverage}'

cd ~/datasets/bladder_cfdna/wxs/statistics
using QualityMetrics
samples = [replace(s, ".general", "") for s in readdir() if endswith(s, ".general")];
plot_quality_metrics(samples)











# IDENTIFY SOMATIC MUTATIONS AND GERMLINE VARIANTS IN WXS
cd ~/datasets/bladder_cfdna/wxs/alignments
mkdir ../mutations
echo X `seq 22` Y | parallel -n10 'mutato call2 --region=chr${x} --alt-reads=5 --alt-frac=0.02 --min-mapq=0 ~/homo_sapiens/hg38.fa *.bam > ../mutations/chr${x}.vcf'

cd ~/datasets/bladder_cfdna/wxs/mutations
cat chr1.vcf <(cat chr{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf | grep -v 'CHROM') > variants.vcf

# Generate list of tumor-normal pairs
cd ~/datasets/bladder_cfdna/wxs
samples = [replace(s, ".bam", "") for s in readdir("alignments")
	if endswith(s, ".bam")];
out = open("tumor_normal_pairs.txt", "w");
write(out, "TEST\tREF\n");
for wbc in filter(s -> "-WBC" in s, samples)
	patient = replace(wbc, r"-WBC.*", "")
	cancers = filter(s -> startswith(s, patient) && !("-WBC" in s), samples)
	for cancer in cancers; write(out, "$(cancer)\t$(wbc)\n"); end
	if isempty(cancers); write(out, "\t$(wbc)\t-\n"); end
end
close(out)

# Somatic mutations
cd ~/datasets/bladder_cfdna/wxs/mutations
variant somatic --alt-reads=8 --alt-frac=0.05 --test-ref-ratio=5 --test-bg-ratio=30 --ref-reads=20 --min-sidedness=15 --min-mapq=10 variants.vcf ../tumor_normal_pairs.txt | variant predict effect - | variant protein altering - > somatic_protein_altering.tmp
variant nearby indels variants.vcf | variant somatic --alt-reads=8 --alt-frac=0.10 --test-ref-ratio=10 --test-bg-ratio=50 --ref-reads=20 --min-sidedness=25 --min-mapq=30 - ../tumor_normal_pairs.txt | variant predict effect - | variant protein altering --invert - | variant discard sketchy silent - > somatic_silent.tmp
cat somatic_protein_altering.tmp <(tail -n +2 somatic_silent.tmp) | sort -k1,1V -k2,2n | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jls 0.005 > somatic.vcf

# Germline variants
cd ~/datasets/bladder_cfdna/wxs/mutations
variant germline --alt-reads=5 --alt-frac=0.15 --bg-ratio=20 variants.vcf WBC | variant predict effect - | variant protein altering - | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jls 0.005 | variant annotate - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls | variant discard if frequency above - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls 0.005 | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jls | variant annotate - ~/homo_sapiens/clinvar-2019-04-08_hg38.jls > germline.vcf


# Heterozygous germline SNPs
variant heterozygous snps --min-depth=50 variants.vcf 'WBC|enign' | variant discard indels - | variant annotate - ~/homo_sapiens/kaviar_2016-02-04_hg38.jls | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jls | egrep 'CHROM|KAVIAR|ExAC' > hetz_snps.vcf




# WHOLE EXOME COPY NUMBER ANALYSIS
cd ~/datasets/bladder_cfdna/wxs/alignments
mkdir ../coverage
coverage grid ~/homo_sapiens/hg38.chrom.sizes 1000 > ../coverage/grid.bed
echo *.bam | parallel -n8 '[ ! -e ../coverage/${x/.bam/.tsv} ] && sam count $x ../coverage/grid.bed > ../coverage/${x/.bam/.tsv}'
fasta gc content ~/homo_sapiens/hg38.fa ../coverage/grid.bed > ../coverage/grid.gc

cd ~/datasets/bladder_cfdna/wxs/coverage
copynum call genomewide --snp-median-decimate=5 --gc-fractions=grid.gc --report-dir=~/tmp/ --hetz-snps=../mutations/hetz_snps.vcf grid.bed ../tumor_normal_pairs.txt







# FIT COPY NUMBER MODEL BASED ON LOGRATIOS
cd ~/datasets/bladder_cfdna/wxs/igv_tracks
clonality try model T-001-1st-cfDNA_logratio.igv 0.3 -0.35









####################################
# COMBINED ANALYSIS OF BOTH PANELS #
####################################

# CREATE UNIFIED FOLDER CONTAINING TARGETED PANEL AND WXS SAMPLES
cd ~/datasets/bladder_cfdna/panel/alignments
for x in *.bam*; do ln -sf ../panel/alignments/${x} ../../alignments/${x/.bam/-Panel.bam}; done
cd ~/datasets/bladder_cfdna/wxs/alignments
for x in *.bam*; do ln -sf ../wxs/alignments/${x} ../../alignments/${x/.bam/-WXS.bam}; done



# CALL SOMATIC MUTATIONS
cd ~/datasets/bladder_cfdna/alignments
echo X `seq 22` Y | parallel -n12 'variant call --alt-reads=5 --alt-frac=0.01 --min-mapq=0 --region=chr${x} ~/homo_sapiens/hg38.fa *.bam > ../mutations/chr${x}.vcf'
cat chr1.vcf <(cat chr{2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf | grep -v 'CHROM') > variants.vcf

# Generate list of tumor-normal pairs
cd ~/datasets/bladder_cfdna/mutations
line = readline(open("variants.vcf")); headers = split(rstrip(line), '\t');
samples = headers[findone(headers .== "NOTES")+1:end];
wxs_out = open("tumor_normal_pairs.wxs.txt", "w");
panel_out = open("tumor_normal_pairs.panel.txt", "w");
for wxs_cfdna in filter(s -> endswith(s, "-cfDNA-WXS"), samples);
	patient = replace(wxs_cfdna, r"-(1st|2nd|3rd|4th)-cfDNA-WXS", "")
	write(wxs_out, "$(wxs_cfdna)\t$(patient)-WBC-WXS\n")
	write(panel_out, "$(replace(wxs_cfdna, "-WXS", "-Panel"))\t$(patient)-WBC-Panel\n")
end
for wbc in filter(s -> endswith(s, "-WBC-WXS"), samples)
	write(wxs_out, "\t$(wbc)\n")
end
for wbc in filter(s -> endswith(s, "-WBC-Panel"), samples)
	write(panel_out, "\t$(wbc)\n")
end
close(wxs_out); close(panel_out);

# Somatic mutations
cd ~/datasets/bladder_cfdna/mutations
variant nearby indels variants.vcf | variant somatic --alt-reads=10 --alt-frac=0.01 --test-ref-ratio=3 --test-bg-ratio=25 --ref-reads=10 - tumor_normal_pairs.panel.txt | variant somatic --keep-old-stars --alt-reads=8 --alt-frac=0.1 --test-ref-ratio=3 --test-bg-ratio=20 --ref-reads=15 - tumor_normal_pairs.wxs.txt | variant predict effect - | variant mappability - ~/homo_sapiens/mappability_170bp_hg38.jld 90 > somatic.tmp
cat <(variant protein altering somatic.tmp) <(variant protein altering --invert somatic.tmp | tail -n +2 | grep -v INDEL | grep -v Mappability) | sort -k1,1V -k2,2n | variant annotate - ~/homo_sapiens/cosmic_77_hg38.jld | variant annotate - ~/homo_sapiens/exac_0.3_hg38.jld | variant discard if frequency above - ~/homo_sapiens/exac_0.3_hg38.jld 0.005 > somatic.vcf











# COMPARE MUTANT ALLELE FRACTIONS BETWEEN THE TWO PANELS
cd ~/datasets/bladder_cfdna/mutations
vcf = read_vcf(pipeline(`variant discard blacklisted somatic.vcf ../panel/mutations/blacklist.tsv`, `variant discard blacklisted - ../wxs/mutations/blacklist.tsv`));
sample = "T-002-3rd";
panel_s = findone(vcf.sample .== "$(sample)-cfDNA-Panel");
wxs_s = findone(vcf.sample .== "$(sample)-cfDNA-WXS");
order = sortperm(vcf.total[:, wxs_s], rev=true);
order = order[vcf.star[order, panel_s]];

x = (vcf.alt[:, panel_s] ./ vcf.total[:, panel_s])[vcf.star[:, panel_s]];
y = (vcf.alt[:, wxs_s] ./ vcf.total[:, wxs_s])[vcf.star[:, panel_s]];
point_size = min(vcf.total[:, wxs_s], vcf.total[:, panel_s])[vcf.star[:, panel_s]];
point_size = sqrt(clamp(point_size, 20, 100) / 15);
figure("~/plot.pdf", size=(2,2));
scatter_plot(x, y, point_size=point_size, square=true, xrange=(0, 0.4), yrange=(0, 0.4));

