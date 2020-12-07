from json import load
from snakemake import shell

json_file = snakemake.input[0]
cmd  = snakemake.params[0]
sra_file  = snakemake.output[0]

print("DEBUG:", json_file, cmd, sra_file, sep = "\n")

meta 	= load(open(json_file))
print("DEBUG:", type(meta[0]).__name__, meta[[0]], sep = "\n")
sra_run = meta[0]["RUN_accession"]

prefetch_str = " %s -o %s"

prefetch = cmd + prefetch_str %(sra_run, sra_file)

print(prefetch)
shell(prefetch)