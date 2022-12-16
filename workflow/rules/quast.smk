# this rule runs quast on the assemblies coming out of flye, medaka, homopolish and polypolish

rule quast:
    output:
        html    = "results/main/quast_reports/html_reports/{barcode}_quast.html"
    input:
        flye        = "results/intermediate/{barcode}/flye_assembly/{barcode}_assembly_flye.fasta",
        medaka      = "results/intermediate/{barcode}/medaka_assembly/{barcode}_assembly_flye_medaka.fasta",
        homopolish  = "results/intermediate/{barcode}/homopolish_assembly/{barcode}_assembly_flye_medaka_homopolish.fasta",
        polypolish  = "results/intermediate/{barcode}/polypolish_assembly/{barcode}_assembly_flye_medaka_polypolish.fasta",
        ref         = "data/reference/GCF_000006765.1_ASM676v1_genomic.fna", #PAO1 reference genome
        annotation  = "data/reference/GCF_000006765.1_ASM676v1_genomic.gff" #PAO1 reference genome annotation
    conda: "../envs/quast.yml"
    params: 
        outdir      = "results/main/quast_reports/full_outputs/{barcode}_quast"
    threads: 1
    shell:
        r"""
            quast -l flye,flye_medaka,flye_medaka_homopolish,flye_medaka_polypolish \
            {input.flye} {input.medaka} {input.homopolish} {input.polypolish} \
            -r {input.ref} -t {threads} --features {input.annotation} -o {params.outdir}

            cp {params.outdir}/report.html {output.html}
         """
