# This rule assembles raw reads using flye

rule flye:
    output:
        assembly    = "results/intermediate/{barcode}/flye_assembly/{barcode}_assembly_flye.fasta",
    input:
        filtered_fq = "results/intermediate/{barcode}/filtlong/{barcode}_filtered.fq.gz"
    params: 
        outdir = "results/intermediate/{barcode}/flye_assembly" 
    conda: "../envs/flye.yml"
    threads: 1
    shell:
        r"""
        flye -o {params.outdir} --threads {threads} --nano-hq {input.filtered_fq}
        mv {params.outdir}/assembly.fasta {output.assembly}
         """

