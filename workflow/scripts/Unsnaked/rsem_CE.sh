sample=GSM3960736
species=9606
cores=15

rsem-calculate-expression --paired-end \
                           --star \
                           --star-gzipped-read-file \
                           --star-output-genome-bam \
                           --paired-end \
                           -p $cores \
                           Results/Rawdata/Fastq/$sample'_1.fastq.gz' \
                           Results/Rawdata/Fastq/$sample'_2.fastq.gz' \
                           Results/References/$species/indices/rsem/ \
                           Results/Counts/$species/$sample_quant

