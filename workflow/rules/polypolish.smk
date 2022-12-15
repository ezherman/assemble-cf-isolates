# This rule further polishes the assemblies created using flye and medaka
# Using short-read sequencing data

rule polypolish:
    output:
        working_dir    = temp(directory("{barcode}_polypolish_working_dir")),
        polypolished   = "results/intermediate/{barcode}/polypolish_assembly/{barcode}_assembly_flye_medaka_polypolish.fasta"
    input:
        medaka          = "results/intermediate/{barcode}/medaka_assembly/{barcode}_assembly_flye_medaka.fasta",
        short_fq1       = "data/short-read/{barcode}_1.fastq.gz",
        short_fq2       = "data/short-read/{barcode}_2.fastq.gz"
    conda: "../envs/polypolish.yml"
    threads: 4
    shell:
        r"""
            cd {output.working_dir}

            bwa index {input.medaka}

            bwa mem -t {threads} -a {input.medaka} > alignments_1.sam
            bwa mem -t {threads} -a {input.medaka} > alignments_2.sam

            polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam \
            --out1 filtered_1.sam --out2 filtered_2.sam

            polypolish {input.medaka} filtered_1.sam filtered_2.sam > {output.polypolished}
         """
