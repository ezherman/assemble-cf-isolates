# this script converts the filenames of short-read data from NCBI
# to names that match the naming of the long-read sequencing data used in this workflow
# the script requires an input spreadsheet defining the file name conversions
# the first column has RefSeq naming, the second column has names used in the workflow

# run this script from the root directory ('assemble-cf-isolates')

file = "data/isolate_info/short_read_file_name_conversions.csv"

import pandas as pd
import glob
import os

name_conversion = pd.read_csv(file)

#shorten the gcf names
name_conversion["gcf"] = [x.split(".1")[0] for x in name_conversion.gcf]

name_conversion_dict = dict(zip(name_conversion.gcf, name_conversion.workflow_id))

os.chdir("data/short-read")

for file in glob.glob("*.fastq.gz"):
	read_num = file.split(".fastq.gz")[0][-2:] #equals _1 or _2
	os.rename(file, name_conversion_dict[file.split(".1")[0]] + read_num + ".fastq.gz")


print("file names have been converted to:", glob.glob("*.fastq.gz"))
