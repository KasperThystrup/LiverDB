from pandas import read_csv
from snakemake import shell

# print(snakemake.input, snakemake.output, snakemake.params, sep  = "\n")

metadata_file = snakemake.input[0]
idx_dir = snakemake.input[1]
fastq_1 = snakemake.input[2]
fastq_2 = snakemake.input[3]
cmd = snakemake.params[0]
rsem_prefix = snakemake.params[1]
rsem_genes = snakemake.output[0]
rsem_isoforms = snakemake.output[1]
threads = snakemake.threads

meta = read_csv(metadata_file)

layout = meta["LibraryLayout"][0]

# --star-gzipped-read-file \
rsem_str = " %s \
            --star \
            --num-threads %s \
            --star-output-genome-bam \
            %s \
            %s/rsem \
            %s" # Prefix

with open(fastq_1, "r") as fq_1:
      header = fq_1.readline().strip()
      print(header)
      if header != "discarded":
            if layout == "SINGLE":
            	layout_option = ""
            	inputs = "%s" %fastq_1
            elif layout != "PAIRED":
            	exit("rsem: Layout type could not be interpreted; %s" %layout)
            else:
            	layout_option = "--paired-end"
            	inputs = "%s %s" %(fastq_1, fastq_2)

            rsem = cmd + rsem_str %(layout_option, threads, inputs, idx_dir, rsem_prefix)
      else:
            rsem = "echo discarded > %s ; echo discarded > %s" %(rsem_genes, rsem_isoforms)

print("DEBUG:" + rsem)
shell(rsem)