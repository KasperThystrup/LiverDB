from pandas import read_csv
from snakemake import shell

metadata_file = snakemake.input[0]
idx_dir = snakemake.input[1]
fastq_1 = snakemake.input[2]
fastq_2 = snakemake.input[3]
cmd = snakemake.params[0]
rsem_prefix = snakemake.params[2]
rsem_genes = snakemake.output[0]
rsem_isoforms = snakemake.output[1]
threads = snakemake.threads

meta = read_csv(metadata_file)

layout = meta["LibraryLayout"][0]

rsem_str = " %s \
            --star \
            --num-threads %s \
            --star-output-genome-bam \
            --star-gzipped-read-file \
            %s \
            %s/rsem \
            %s" # Prefix


if layout == "SINGLE":
	layout_option = ""
	inputs = "%s" %fastq_1
elif layout != "PAIRED":
	exit("rsem: Layout type could not be interpreted; %s" %layout)
else:
	layout_option = "--paired-end"
	inputs = "%s %s" %(fastq_1, fastq_2)

rsem = cmd + rsem_str %(layout_option, threads, inputs, idx_dir, rsem_prefix)

print("DEBUG:" + rsem)
shell(rsem)