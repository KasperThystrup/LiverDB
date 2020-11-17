from json import load
from snakemake import shell
from sys import exit

json_file		= snakemake.input[0]
fastq_1			= snakemake.input[1]
fastq_2			= snakemake.input[2]
cmd				= snakemake.params[0]
ref_dir 		= snakemake.params[1]
threads			= snakemake.params[2]
star_prefix		= snakemake.params[3]
bam				= snakemake.output[0]
tx_bam			= snakemake.output[1]
chim			= snakemake.output[2]
log_final		= snakemake.output[3]
log_out			= snakemake.output[4]
log_progress	= snakemake.output[5]
out_tab			= snakemake.output[6]

json 	= load(open(json_file))
tax_id	= json["TAXON_ID"]
layout	= json["LAYOUT"]

star_idx = "%s/%s/STAR" %(ref_dir, tax_id)

star_str = """ \
--runThreadN			%s						\
--runMode		 		alignReads				\
--outSAMtype 			BAM SortedByCoordinate	\
--genomeDir 			%s						\
--outFileNamePrefix 	%s						\
--chimSegmentMin 		20						\
--quantMode 			TranscriptomeSAM		\
--readFilesIn  			%s %s """


if layout == "SINGLE":
	fastq_2 = " ; touch %s " %fastq_2
elif layout != "PAIRED":
	exit("star_align: Layout type could not be interpreted; %s" %layout)

star = cmd + star_str %(threads, star_idx, star_prefix, fastq_1, fastq_2)

print(star.replace("\t", " "))
shell(star)