# homopolish requires the bacteria sketch file. This rule downloads the file.

rule download_homopolish_sketch:
    output:
        sketch      = "resources/homopolish/bacteria.msh"
    params:
        filename    = "bacteria.msh"
    resources: time_min=180
    shell:
        r"""
            curl http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/{params.filename}.gz -o {params.filename}.gz
            gunzip {params.filename}.gz
            mv {params.filename} {output.sketch}
         """
