from pandas import read_csv
from snakemake import shell

csv_file = snakemake.input[0]
cmd  = snakemake.params[0]
sra_file  = snakemake.output[0]

# print("DEBUG:", csv_file, cmd, sra_file, sep = "\n")
meta = read_csv(csv_file)

# print("DEBUG:", meta["Run"][0], meta["LibraryLayout"][0], sep = "\n")
sra_run = meta["Run"][0]

prefetch_str = " %s -o %s"

prefetch = cmd + prefetch_str %(sra_run, sra_file)

print(prefetch)
shell(prefetch)