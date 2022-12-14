# This rule further polishes the assemblies created using flye and medaka

rule homopolish:
    output:
        renamed_input   = temp("input_{barcode}.fasta"),
        homopolished    = "results/intermediate/{barcode}/homopolish_assembly/{barcode}_assembly_flye_medaka_homopolish.fasta",
    input:
        medaka          = "results/intermediate/{barcode}/medaka_assembly/{barcode}_assembly_flye_medaka.fasta",
        sketch          = "data/homopolish/bacteria.msh"
    params: 
        outdir = "results/intermediate/{barcode}/homopolish_assembly",
        model  = "R9.4.pkl",
        prefix = "input"
    conda: "homopolish" #this refers to a conda environment that already exists
                        #because this conda environment cannot be produced from a yaml
                        #with strict channel priority
    threads: 4 #set to max threads now to prevent multiple NCBI downloads at the same time, can improve with a resources section for cluster
    shell:
        r"""
            # homopolish outputs yourgenome_homopolished.fasta, where yourgenome is the name of the input file
            # make a temporary copy of the input file, so that the output name is predictable for generic input
            cp {input.medaka} {output.renamed_input}

            # may want to update this to use a local NCBI database, so that all isolates are polishes using the same genomes
            homopolish polish -a {output.renamed_input} -s {input.sketch} -m {params.model} -o {params.outdir} -t {threads}
            mv {params.outdir}/{params.prefix}_{wildcards.barcode}_homopolished.fasta {output.homopolished}
         """