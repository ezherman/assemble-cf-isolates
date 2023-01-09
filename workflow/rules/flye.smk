# This rule assembles raw reads using flye

rule flye:
    output:
        assembly    = "results/intermediate/{barcode}/flye_assembly/{barcode}_assembly_flye.fna",
    input:
        filtered_fq = "results/intermediate/{barcode}/filtlong/{barcode}_filtered.fq.gz"
    params: 
        outdir = "results/intermediate/{barcode}/flye_assembly",
        guppy_v = config.get("guppy_version", "6.3.9") 
    conda: "../envs/flye.yml"
    threads: 8
    shell:
        r"""
        
        # Flye command is conditional on the version of Guppy (whether it's >= v5)
        guppy_v={params.guppy_v}
        if [ ${{guppy_v:0:1}} -lt 5 ]
        then
            flye -o {params.outdir} --threads {threads} --nano-raw {input.filtered_fq}
        else
            flye -o {params.outdir} --threads {threads} --nano-hq {input.filtered_fq}
        fi

        mv {params.outdir}/assembly.fasta {output.assembly}
         """

