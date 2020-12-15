species=9606
gtf=Results/References/$species/$species.gtf
fa=Results/References/$species/$species.fa
cores=15

rsem-prepare-reference --gtf $gtf \
--star \
--num-threads $cores \
$fa \
Results/References/$species/indices/rsem/
