# This rule polishes the assemblies created using flye

rule medaka:
    output:
        polished    = "results/intermediate/{barcode}/medaka_assembly/{barcode}_assembly_flye_medaka.fasta"
    input:
        assembly    = "results/intermediate/{barcode}/flye_assembly/{barcode}_assembly_flye.fasta",
        raw_cat_fq  = "results/intermediate/cat_fastq_files/{barcode}_cat.fq.gz"
    params: 
        outdir = "results/intermediate/{barcode}/medaka_assembly",
        model  = "r941_min_sup_g507" 
    conda: "../envs/medaka.yml"
    threads: 1
    shell:
        r"""
            export OMP_NUM_THREADS={threads} #limit the number of threads that medaka uses
        
            medaka consensus -i {input.raw_cat_fq} -d {input.assembly} -o {params.outdir} -m {params.model}
            mv {params.outdir}/consensus.fasta {output.polished}
         """
