from json import load
from snakemake import shell
from sys import exit

json_file	= "results/Metadata/SRR8094770.json"
sra_file	= "results/SRA/SRR8094770.sra"
com			= "/disk2/Resources/software/sratoolkit.2.10.5-centos_linux64/bin/fastq-dump.2.10.5 --split-files -O %s %s"
fq_dir		= "results/Fastq/"
fastq_dummy = fq_dir + "SRR8094770.fastq"

json 	= load(open(json_file))
layout	= json["LAYOUT"]

fastq_dump = com %(fq_dir, sra_file)
if layout == "SINGLE":
	dummy_com = " ; touch %s" %fastq_dummy
	fastq_dump += dummy_com
elif layout != "PAIRED":
	exit("Layout type could not be interpreted; %s" %layout)

shell(fastq_dump)