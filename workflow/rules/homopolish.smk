# This rule further polishes the assemblies created using flye and medaka

rule homopolish:
    output:
        renamed_input   = temp("input_{barcode}.fna"),
        homopolished    = "results/main/{barcode}/homopolish_assembly/{barcode}_assembly_flye_medaka_homopolish.fna",
    input:
        medaka          = "results/intermediate/{barcode}/medaka_assembly/{barcode}_assembly_flye_medaka.fna",
        sketch          = "resources/homopolish/bacteria.msh"
    params: 
        outdir = "results/intermediate/{barcode}/homopolish_assembly",
        model  = config.get("homopolish_model", "R9.4.pkl"),
        prefix = "input"
    conda: "homopolish" #specifies the homopolish env, rather than a yml file
    threads: 8
    resources: mem_mb=4000 
    shell:
        r"""
            # homopolish outputs yourgenome_homopolished.fasta, where yourgenome is the name of the input file
            # make a temporary copy of the input file, so that the output name is predictable for generic input
            cp {input.medaka} {output.renamed_input}

            # may want to update this to use a local NCBI database, so that all isolates are polished using the same genomes
            homopolish polish -a {output.renamed_input} -s {input.sketch} -m {params.model} -o {params.outdir} -t {threads}
            mv {params.outdir}/{params.prefix}_{wildcards.barcode}_homopolished.fasta {output.homopolished}
         """
