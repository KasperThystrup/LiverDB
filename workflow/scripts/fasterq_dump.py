from pandas import read_csv
from os.path import basename, dirname, exists
from snakemake import shell

metadata_file = snakemake.input[0]
sra_file = snakemake.input[1]
cmd = snakemake.params[0]
fastq = snakemake.output[0]
fastq_1 = snakemake.output[1]
fastq_2 = snakemake.output[2]
threads = snakemake.threads


if threads > 10:
	threads = 10  ## More cores dirsupts the run 
meta = read_csv(metadata_file)

layout = meta["LibraryLayout"][0]

fasterq_dump_str = " --force --split-3 --include-technical --temp %s --threads %s --outdir %s %s"

temp_dir = dirname(fastq_1)

if layout == "SINGLE":
	dummy_com = " ; mv %s %s ; touch %s ; echo discarded > %s"
	fasterq_dump_str += dummy_com
	fasterq_dump = cmd + fasterq_dump_str %(temp_dir, threads, temp_dir, sra_file, fastq, fastq_1, fastq_2, fastq)
elif layout != "PAIRED":
	exit("fasterq_dump: Layout type could not be interpreted; %s" %layout)
else:
	fasterq_dump_str += " ; echo discarded > %s "
	fasterq_dump = cmd + fasterq_dump_str %(temp_dir, threads, temp_dir, sra_file, fastq)

print("DEBUG:" + fasterq_dump)
shell(fasterq_dump)

if layout == "PAIRED" and not exists(fastq_1):
	sample_id = basename(sra_file).replace("\.sra", "")
	print("Paired end dataset:", id, "")

	shell("touch %s && touch %s && echo %s >> unpaired_samples.txt" %(fastq_1, fastq_2, sample_id))