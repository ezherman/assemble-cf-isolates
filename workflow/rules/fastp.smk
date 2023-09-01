rule fastp:
    output:
        fq1         = "results/intermediate/{barcode}/fastp/{barcode}_filt_1.fq.gz",
        fq2         = "results/intermediate/{barcode}/fastp/{barcode}_filt_2.fq.gz",
        unpaired    = "results/intermediate/{barcode}/fastp/{barcode}_unpaired.fq.gz"
    input:
        fq1         = "data/short-read/{barcode}_1.fastq.gz",
        fq2         = "data/short-read/{barcode}_2.fastq.gz",
    conda: "../envs/fastp.yaml"
    threads: 8
    shell:
        r"""
            fastp --in1 {input.fq1} --in2 {input.fq2} \
            --out1 {output.fq1} --out2 {output.fq2} \
            --unpaired1 {output.unpaired} --unpaired2 {output.unpaired} \
            -w {threads} 
         """