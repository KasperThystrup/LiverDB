cores=36
species=10090
avg_length=100

mkdir Results/References/$species/STAR

STAR --runThreadN $cores \
--runMode genomeGenerate \
--genomeDir Results/References/$species/STAR \
--genomeFastaFiles Results/References/$species/$species*.fa \
--sjdbGTFfile Results/References/$species/$species.gtf \
--sjdbOverhang $avg_length \
--outFileNamePrefix Results/References/$species/STAR/$species_
