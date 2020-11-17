from scripts import sradata
# from srahelp import getftpfilesize
import json

# Setup snakemake background parameters
configfile: "config/config.yaml"
workdir: config["workdir"]

# Defining software
prefetch                    = config["prefetch"]
fastq_dump                  = config["fastq_dump"]
star                        = config["star"]
samtools_sort               = config["samtools_sort"]
htseq_count                 = config["htseq_count"]


# Rule paths and file names --------------------------------------------------
ref_dir                     = config["ref_dir"]
max_threads                 = config["max_threads"]
sub_threads                 = config["sub_threads"]
threshold                   = config["threshold"]
strandedness                     = config["strandedness"]
species                     = config["species"]
samples                     = config["samples"]

output_dir                  = config["output_dir"] + "results/"

log_dir                    = config["output_dir"] + "info/logs/"
benchmark_dir              = config["output_dir"] + "info/benchmark/"

metadata_files              = output_dir + "Metadata/{sample}.json"
sra_files                   = output_dir + "Rawdata/SRA/{sample}.sra"
fastq_dir                   = output_dir + "Rawdata/Fastq/"
fastq_1                     = fastq_dir  + "{sample}_1.fastq"
fastq_2                     = fastq_dir  + "{sample}_2.fastq"
align_dir                   = output_dir + "Rawdata/Aligned/{sample}_"
gene_bam                    = align_dir  + "Aligned.sortedByCoord.out.bam"
tx_bam                      = align_dir  + "Aligned.toTranscriptome.out.bam"
chim                        = align_dir  + "Chimeric.out.junction"
log_final                   = align_dir  + "Log.final.out"
log_out                     = align_dir  + "Log.out"
log_progress                = align_dir  + "Log.progress.out"
out_tab                     = align_dir  + "SJ.out.tab"
tx_bam_sorted               = align_dir  + "Aligned.toTranscriptome.sortedByCoord.out.bam"
# gene_files                  = output_dir + "Rawdata/HTSeq_count/{sample}_genes.tsv"
# transcript_files            = output_dir + "Rawdata/HTSeq_count/{sample}_transcripts.tsv"
gene_counts_files           = output_dir + "Rawdata/Counts/{species}/{sample}_genes.tsv"
transcript_counts_files     = output_dir + "Rawdata/Counts/{species}/{sample}_transcripts.tsv"
gene_filtered_files         = output_dir + "Filetered_counts/{species}/Samples/{sample}_gene_count.tsv"
transcript_filtered_files   = output_dir + "Filetered_counts/{species}/Samples/{sample}_transcript_count.tsv"
gene_normalized_files       = output_dir + "Counts/{species}/gene_normalized.tsv"
transcript__normalized_files   = output_dir + "Counts/{species}/transcript_normalized.tsv"

# Start rules ----------------------------------------------------------------
# All file rules are included here for trouble shooting purposes. 
#In the cases of single end make blank transcripts file.
rule all:
    input:
        # expand(gene_files, sample = samples),
        # expand(transcript_files, sample = samples),
        expand(gene_filtered_files, sample = samples, species = species),
        expand(transcript_filtered_files, sample = samples, species = species),
        # expand(fastq_1, sample=samples),
        # expand(fastq_2, sample=samples),
        expand(gene_bam, sample=samples),
        expand(tx_bam_sorted, sample=samples)
        # expand(tx_bam,  sample=samples),
        # expand(chim, sample=samples),
        # expand(log_final, sample=samples),
        # expand(log_out, sample=samples),
        # expand(log_progress, sample=samples),
        # expand(out_tab, sample=samples)
        

"""
The getmetadata rule downloads various metadata from NCBI regarding the SRA in question. 
This is used for determining various future steps in analysis such as to align the sample 
as a paired end of single end, etc. 
"""
rule metadata:
    output:
        metadata_files
    params:
        "{sample}"
    threads:
        1  ## No multi-threading support, however simultanous samples may lead to timeout errors.
    log:
        log_dir + "{sample}.metadata.log"
    benchmark:
        benchmark_dir + "{sample}.metadata.bench"
    script:
        "scripts/sradata.py"
        
"""
The download sequence data rule downloads the SRA using the link stored in the METADATA. 
"""
rule srafetch: 
    input: 
        json_file   = metadata_files
    output:
        sra_file    = temp(sra_files)
    params:
        cmd         = prefetch
    conda:
      "envs/sra_tools.yaml"
    threads:
        1  ## No multi-threading support, instead samples are run simultanously
    log:
        log_dir + "{sample}.srafetch.log"
    benchmark:
        benchmark_dir + "{sample}.srafetch.bench"
    script:
        "scripts/sraprefetch.py"

"""
The fastqdump rule simply extracts the fastq files from the SRA files. 
Fastq files are large, so these are deleted when the rules calling fastqdump are finished.

Links trouble shooting the problem: 
https://edwards.sdsu.edu/research/fastq-dump/
https://www.biostars.org/p/156909/
I have no idea why -M 0 keeps the lines equal but it does. 
Nevermind -M 0 allows for sequences short enought to be illegal in STAR, FML

@TODO: Think about adding output as gziped fastq. 
fastq-dump --gzip --skip-technical  --readids --dumpbase --split-files --clip sra_filename
"""
rule fastqdump:
    input:
        json_file   = metadata_files,
        sra_file    = sra_files
    output:
        fastq_1     = temp(fastq_1),
        fastq_2     = temp(fastq_2)
    params:
        cmd         = fastq_dump,
        fastq_dir   = fastq_dir
    conda:
      "envs/sra_tools.yaml"
    threads:
        1  ## No multi-threading support, instead samples are run simultanously
    log:
        log_dir + "{sample}.fastqdump.log"
    benchmark:
        benchmark_dir + "{sample}.fastqdump.bench"
    script:
        "scripts/fastq_dump.py"


"""
STAR_align runs guided alignments on the fastq files extracted from SRAs. 
It needs an index file for the guided alignment. So these must be built beforehand.
See the STAR documentation for details. 
"""
rule STAR_align:
    input:
        json_file   = metadata_files,
        fastq_1     = fastq_1,
        fastq_2     = fastq_2
    output:
       gene_bam     = gene_bam,
       # gene_bam     = temp(gene_bam),
       tx_bam       = temp(tx_bam),
       chim         = chim,
       log_final    = log_final,
       log_out      = log_out,
       log_progress = temp(log_progress),
       out_tab      = out_tab
    params:
        cmd         = star,
        ref_dir     = ref_dir,
        threads     = sub_threads,
        star_prefix = align_dir
    conda:
      "envs/star.yaml"
    threads:
        max_threads
    log:
        log_dir + "{sample}.star_align.log"
    benchmark:
        benchmark_dir + "{sample}.star_align.bench"
    script:
        "scripts/star_align.py"

"""
In case of paired end data, Transcriptome output must be sorted by coordinate,
as it is a requirement for HTSeq-counts to have data sorted by reads or coordinates
"""
rule tx_sort:
    input:
        json_file       = metadata_files,
        tx_bam          = tx_bam
    output:
        tx_bam_sorted   = tx_bam_sorted
        # tx_bam_sorted   = temp(tx_bam_sorted)
    params:
        cmd             = samtools_sort
    conda:
      "envs/samtools.yaml"
    threads:
        sub_threads
    log:
        log_dir + "{sample}.tx_sort.log"
    benchmark:
        benchmark_dir + "{sample}.tx_sort.bench"
    script:
        "scripts/tx_sort.py"


"""
makeHTSeqTables extracts counts from the STAR alignment for genes if a sigle end SRA is used, 
or genes and transcripts if a paired end SRA is used. 
"""
rule count_reads:
    input:
        json_file               = metadata_files,
        gene_bam                = gene_bam,
        tx_bam_sorted           = tx_bam_sorted
    output:
       gene_counts_files        = gene_counts_files,
       transcript_counts_files  = transcript_counts_files
    params:
        cmd                     = htseq_count,
        ref_dir                 = ref_dir,
        strandedness            = strandedness
    conda:
      "envs/samtools.yaml"
    threads:
        sub_threads
    log:
        log_dir + "{species}.{sample}.count_reads.log"
    benchmark:
        benchmark_dir + "{species}.{sample}.count_reads.bench"
    script:
        "scripts/count_reads.py"


"""
Samples with high zero inflation could skewer the normalization.
Zero contents in the samples are screened for zero contents,
and samples with a higher zero contents than a set threshold are dismissed.
"""
rule zero_intolerance:
    input:
        gene_counts_files           = gene_counts_files,
        transcript_counts_files     = transcript_counts_files
    output:
        gene_filtered_files         = gene_filtered_files,
        transcript_filtered_files   = transcript_filtered_files
    params:
        threshold                   = threshold,
        threads                     = sub_threads
    conda:
        "envs/R.yaml"
    threads:
        max_threads
    script:
        "scripts/zero_intolerance.R"

# """
# Normalize samples counts
# """
# rule normalize_counts:
#     input:
#         gene_filtered_files         = gene_filtered_files,
#         transcript_filtered_files   = transcript_filtered_files
#     output:
#         gene_normalized_counts      = gene_normalized_counts,
#         transcript_normalized_files = transcript_normalized_files
#     params:
#         threads                     = sub_threads
#     conda:
#         "envs/R.yaml"
#     threads:
#         round(max_threads / sub_threads) ## Allows simultanous sample processing
#     script:
#         "scripts/normalize_counts.R"
