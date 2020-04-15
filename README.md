This repository contains code used to generate the results in study "Circulating tumor DNA versus tumor tissue for genomic stratification in metastatic urothelial carcinoma", currently in review in Nature Communications.

The somatic mutation analysis code under subfolder "mutato" was written using the Rust programming language. For the identification of chromosomal rearrangements, we used our own software called "breakfast":
https://github.com/annalam/breakfast

Software for copy number analysis and post-processing of somatic mutation and germline variant results was written in Julia. This code can be found under the subfolder "julia".

The file "analysis_diary.jl" lists the shell commands and interactive Julia sessions that were used to carry out analyses for the manuscript.

