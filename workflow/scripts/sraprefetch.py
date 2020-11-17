from json import load
from snakemake import shell

json_file = snakemake.input[0]
cmd  = snakemake.params[0]
sra_file  = snakemake.output[0]

json 	= load(open(json_file))
sra_run = json["RUN_accession"]

prefetch = cmd %(sra_run, sra_file)

print(prefetch)
shell(prefetch)