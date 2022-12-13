# assemble-cf-isolates
Snakemake pipeline to assemble ONT reads from P. aeruginosa isolates. 

This workflow filters reads using Filtlong, assembles reads using Flye, polishes assemblies using Medaka and further polishes assemblies using Homopolish.  
The assemblies obtained with Flye, Medaka and Homopolish are compared to a reference P. aeruginosa genome using QUAST. 
