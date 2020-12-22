from os.path import dirname
from snakemake import shell

email = snakemake.params[0]
samples = snakemake.params[1]
output_dir = snakemake.params[2]

for sample in samples:
	metadata_file = output_dir + sample + ".csv"
	print("DEBUG:" + metadata_file)
	shell("esearch -db sra -query %s -email %s | efetch -format runinfo > %s" %(sample, email, metadata_file))
