from json import load
from snakemake import shell
from sys import exit

json_file	= snakemake.input[0]
sra_file	= snakemake.input[1]
cmd			= snakemake.params[0]
fastq_dir	= snakemake.params[1]
fastq_2 	= snakemake.output[1]

json 	= load(open(json_file))
layout	= json["LAYOUT"]

fastq_dump_str = " --split-e --outdir %s %s "

fastq_dump = cmd + fastq_dump_str %(fastq_dir, sra_file)
if layout == "SINGLE":
	dummy_com = " ; touch %s" %fastq_2
	fastq_dump += dummy_com
elif layout != "PAIRED":
	exit("fastq_dump: Layout type could not be interpreted; %s" %layout)

print(fastq_dump)
shell(fastq_dump)