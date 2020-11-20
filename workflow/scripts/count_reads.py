from json import load
from snakemake import shell
from sys import exit

json_file 			= snakemake.input[0]
gene_bam			= snakemake.input[1]
tx_bam				= snakemake.input[2]
cmd					= snakemake.params[0]
ref_dir		 		= snakemake.params[1]
strandedness		= snakemake.params[2]
gene_counts			= snakemake.output[0]
transcript_counts	= snakemake.output[1]

if strandedness:
	strandedness = "yes"
elif not strandedness:
	strandedness = "no"
else:
	strandedness = "reverse"

json 	= load(open(json_file))
tax_id	= json["TAXON_ID"]
layout	= json["LAYOUT"]

HTSeq_idx = "%s/%s/gtf/%s.gtf" %(ref_dir, tax_id, tax_id)

HTSeq_str = """ 	\
--order		pos		\
--stranded	%s		\
--format	bam		\
--type		exon	\
--idattr %s			\
%s %s > %s ; """

HTSeq_count = cmd + HTSeq_str %(strandedness, "gene_id", gene_bam, HTSeq_idx, gene_counts)
if layout == "SINGLE":
	dummy_cmd = "touch %s " %transcript_counts
	HTSeq_count += dummy_cmd
elif layout == "PAIRED":
	HTSeq_count += cmd + HTSeq_str %(strandedness, "transcript_id", tx_bam, HTSeq_idx, transcript_counts)
else:
	exit("count_reads: Layout type could not be interpreted; %s" %layout)

shell(HTSeq_count)