from json import load
from snakemake import shell
from sys import exit

json_file		= "results/Metadata/SRR12845350.json"
fastq_1			= "results/Fastq/SRR12845350_1.fastq"
fastq_2			= "results/Fastq/SRR12845350_2.fastq"
com				= """/disk2/Resources/software/STAR/bin/Linux_x86_64/STAR \
    --runMode alignReads \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir %s \
    --outFileNamePrefix %s \
    --chimSegmentMin 20 \
    --quantMode TranscriptomeSAM \
    --readFilesIn  %s %s """
# idx_dir 		= snakemake.params[1]
star_prefix		= "results/STAR/SRR12845350_"
bam				= "results/STAR/SRR8094770_Aligned.sortedByCoord.out.bam"
tx_bam			= "results/STAR/SRR8094770_Aligned.toTranscriptome.out.bam"
chim			= "results/STAR/SRR8094770_Chimeric.out.junction"
log_final		= "results/STAR/SRR8094770_Log.final.out"
log_out			= "results/STAR/SRR8094770_Log.out"
log_progress	= "results/STAR/SRR8094770_Log.progress.out"
out_tab			= "results/STAR/SRR8094770_SJ.out.tab"

json 	= load(open(json_file))
tax_id	= json["TAXON_ID"]
layout	= json["LAYOUT"]

star_idx = "/disk2/Resources/seqDB/ensembl/release-100/9606/STAR"


if layout == "SINGLE":
	fastq_2 = " ; touch %s" %fastq_2
elif layout != "PAIRED":
	exit("Layout type could not be interpreted; %s" %layout)

star = com %(star_idx, star_prefix, fastq_1, fastq_2)

shell(star)