from pandas import read_csv
from os.path import dirname
from snakemake import shell

metadata_files = snakemake.input[0]
sra_files = snakemake.input[1]
cmd = snakemake.params[0]
fastq_1 = snakemake.output[0]
fastq_2 = snakemake.output[1]
threads = snakemake.threads

#print("DEGBUG:", metadata_files, sra_files, cmd, fastq_1, fastq_2, threads)
meta = read_csv(metadata_files)

layout = meta["LibraryLayout"][0]

fasterq_dump_str = " --force --split-files --temp %s --threads %s --outdir %s %s"

temp_dir = dirname(fastq_1)
if layout == "SINGLE":
	dummy_com = " ; touch %s"
	fasterq_dump_str += dummy_com
	fasterq_dump = cmd + fasterq_dump_str %(temp_dir, threads, temp_dir, sra_files, fastq_2)
elif layout != "PAIRED":
	exit("fasterq_dump: Layout type could not be interpreted; %s" %layout)
else:
	fasterq_dump = cmd + fasterq_dump_str %(temp_dir, threads, temp_dir, sra_files)

print("DEBUG:" + fasterq_dump)
shell(fasterq_dump)
