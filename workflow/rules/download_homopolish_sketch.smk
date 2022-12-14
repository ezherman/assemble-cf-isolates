# homopolish requires the bacteria sketch file. This rule downloads the file.

rule download_homopolish_sketch:
    output:
        sketch      = "data/homopolish/bacteria.msh"
    params:
        filename    = "bacteria.msh"
    shell:
        r"""
            wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/{params.filename}.gz
            gunzip {params.filename}.gz
            mv {params.filename} {output.sketch}
         """