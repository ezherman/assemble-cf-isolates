# This rule polishes the assemblies created using flye

rule medaka:
    output:
        polished    = "results/intermediate/{barcode}/medaka_assembly/{barcode}_assembly_flye_medaka.fna"
    input:
        assembly    = "results/intermediate/{barcode}/flye_assembly/{barcode}_assembly_flye.fna",
        filtered_fq = "results/intermediate/{barcode}/filtlong/{barcode}_filtered.fq.gz",
    params: 
        outdir = "results/intermediate/{barcode}/medaka_assembly",
        model  = config.get("medaka_model", "r941_min_sup_g507") 
    conda: "../envs/medaka.yml"
    threads: 2
    shell:
        r"""
            export TF_NUM_INTEROP_THREADS={threads} #limit the number of threads that medaka uses

            rm -r {params.outdir} #remove any files from the previous medaka run on this data 
            medaka_consensus -i {input.filtered_fq} -d {input.assembly} -o {params.outdir} -m {params.model}
            mv {params.outdir}/consensus.fasta {output.polished}
         """
