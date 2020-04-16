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

The file "analysis_diary.jl" lists the shell commands and interactive Julia sessions used to carry out analyses for the manuscript, allowing replication of the manuscript results.
