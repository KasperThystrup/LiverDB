from json import load
from snakemake import shell
from sys import exit

json_file 		= snakemake.input[0]
tx_bam			= snakemake.input[1]
cmd				= snakemake.params[0]
tx_sorted		= snakemake.output[0]
threads 		= snakemake.threads

json 			= load(open(json_file))
layout			= json["LAYOUT"]

sort_str = " -n -@ %s -o %s -O bam %s"
sort = cmd + sort_str %(threads, tx_sorted, tx_bam)

if layout == "SINGLE":
	sort = "touch %s" %tx_sorted
elif layout != "PAIRED":
	exit("tx_sort: Layout type could not be interpreted; %s" %layout)

print(sort)
shell(sort)