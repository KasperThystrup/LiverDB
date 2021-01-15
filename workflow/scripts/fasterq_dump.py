from pandas import read_csv
from os.path import basename, dirname, exists
from snakemake import shell

metadata_file = snakemake.input[0]
sra_file = snakemake.input[1]
cmd = snakemake.params[0]
fastq_1 = snakemake.output[0]
fastq_2 = snakemake.output[1]
unpaired = snakemake.output[2]
threads = snakemake.threads


if threads > 10:
	threads = 10  ## More cores dirsupts the run 
meta = read_csv(metadata_file)

layout = meta["LibraryLayout"][0]

fasterq_dump_str = " --force --split-files --temp %s --threads %s --outdir %s %s"

if layout == "SINGLE":
	mate_fix = " ; mv %s %s" %(fastq, fastq_1)
	fasterq_dump_str += mate_fix
	fasterq_dump = cmd + fasterq_dump_str %(fastq_dir, threads, fastq_dir, sra_file)
elif layout != "PAIRED":
	exit("fasterq_dump: Layout type could not be interpreted; %s" %layout)
else:
	fasterq_dump = cmd + fasterq_dump_str %(temp_dir, threads, temp_dir, sra_file, fastq)

shell(fasterq_dump)

