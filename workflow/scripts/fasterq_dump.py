from pandas import read_csv
from snakemake import shell

metadata_files = snakemake.input[0]
sra_files = snakemake.input[1]
cmd = snakemake.params[0]
fastq_dir = snakemake.params[1]
fastq_2 = snakemake.output[1]  ## Dummy input if sample are SINGLE-end
threads = snakemake.threads

meta = read_csv(metadata_files)

layout = meta["LibraryLayout"][0]

fasterq_dump_str = " --split-files --outdir %s %s --threads %s"

fasterq_dump = cmd + fasterq_dump_str %(fastq_dir, sra_files, threads)
if layout == "SINGLE":
	dummy_com = " ; touch %s" %fastq_2
	fasterq_dump += dummy_com
elif layout != "PAIRED":
	exit("fasterq_dump: Layout type could not be interpreted; %s" %layout)

shell(fasterq_dump)
