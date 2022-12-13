# This rule concatenates the raw fastq files and filters
# the reads using filtlong

rule filtlong:
    output:
        raw_cat_fq  = "results/intermediate/cat_fastq_files/{barcode}_cat.fq.gz",
        filtered_fq = "results/intermediate/{barcode}/filtlong/{barcode}_filtered.fq.gz",
    input:
        raw_fq      = "data/long-read-seq-sup/{barcode}"
    conda: "../envs/filtlong.yml"
    shell:
        r"""
            cat {input.raw_fq}/*.fastq.gz > {output.raw_cat_fq}
            
            filtlong --min_length 1000 --keep_percent 95 {output.raw_cat_fq} | \
            gzip > {output.filtered_fq}
         """

