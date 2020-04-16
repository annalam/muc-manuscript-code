This repository contains code used to generate the results in study "Plasma ctDNA is a tumor tissue surrogate and enables clinical-genomic stratification of metastatic bladder cancer", currently in review in Nature Communications.

All analyses were run on a computational server running the Linux operating system. No non-standard hardware was utilized.

The somatic mutation analysis code under subfolder "mutato" was written using the Rust programming language, and has been tested to compile with the Rust compiler version 1.39.0. The software can be compiled and installed by typing "cargo install". The compilation process typically takes a few minutes. The software can be tested by running the following command inside the "mutato/example" directory:

	gunzip hg38_chr22.fa.gz
    mutato call2 --alt-reads=20 --alt-frac=0.2 hg38_chr22.fa input.bam > variants.vcf
    diff -u variants.vcf expected_output.vcf

For chromosomal rearrangement analysis, we used our software called "breakfast". The software has been tested to compile with the Rust compiler version 1.39.0. The source code for this software can be found in this public repository:
https://github.com/annalam/breakfast

Software for copy number analysis and post-processing of somatic mutation and germline variant results was written in Julia. This code is found under the subfolder "julia". The code has been tested to work with Julia version 1.3.1, and the following Julia packages:
- Distributions (version 0.22.0)
- HypothesisTests (version 0.8.0)
- PyCall (version 1.91.2)
- KernelDensity (version 0.5.1)
- Loess (version 0.5.0)
- StatsBase (version 0.32.0)
- BioSequences (version 1.1.0)
- BioAlignments (version 1.0.1)

The file "analysis_diary.jl" lists the shell commands and interactive Julia sessions used to carry out genomic analyses of the sequencing data.

All figures for the manuscript were generated in Python 3.7.3 using the packages and libraries listed below. All statistics were performed with the scipy.stats package, except for Cox Proportional-Hazard Models, which were performed with lifelines.

Our Anaconda environment can be recreated by running the following for each package:

    conda install PACKAGE=VERSION

Additionally we used the Spyder IDE beta version 4.0.3b, which can be installed by running:

    conda install -c spyder-ide spyder=4.0.0b3

Python packages used:
- pandas (version 0.24.2)
- matplotlib (version 3.0.3)
- seaborn (version 0.9.0)
- numpy (version 1.16.4)
- scipy (version 1.2.1)
- lifelines (version 0.22.3)
- argparse (version 1.1)

The following Python code was used for determining whether a mutation was protein altering:

	def keepCodingMutations(df_muts):
	    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) |
		(df_muts['EFFECT'].str.contains("Stopgain", regex=False)) |
		(df_muts['EFFECT'].str.contains("Frameshift", regex=False)) |
		(df_muts['EFFECT'].str.contains("Splice", regex=False)) |
		(df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False))]


