from snakemake import shell
configfile: "workflow/config/config.yaml" # Remember to direct to correct config file

sample_file = "accepted_samples.txt"  ## Insert correct path!
Results = config["output_dir"]

with open (samples_files, 'r') as sample_read:
    samples = sample_read.read().splitlines()
[x for x in samples if x.strip()]

compress = config["compress"]
threshold = config["threshold"]


"""
Sample fastq files are compressed, to save disk space.
"""
rule Compress:
    input:
        fastq_1 = rules.Fasterqdump.output.fastq_1,
        fastq_2 = rules.Fasterqdump.output.fastq_2
    output:
        fastq_1 = temp(Results + "Rawdata/Fastq/{sample}_1.fastq.gz"),  ##should be temp
        fastq_2 = temp(Results + "Rawdata/Fastq/{sample}_2.fastq.gz")  ##should be temp
    threads:
        workflow.cores / parallel_processes
    conda:
      "envs/decompress.yaml"
    group:
        "Samples"
    shell:
        "pigz \
            --force \
            --processes {threads} \
            {input.fastq_1} {input.fastq_2}"


"""
Gene and transcript quantifications are calculated using RSEM.
"""
rule Rsem_expression:
    input:
        metadata_file = rules.Metadata.output.metadata_file,
        idx_dir = rules.Rsem_idx.output.idx_dir,
        fastq_1 = rules.Fasterqdump.output.fastq_1,
        fastq_2 = rules.Fasterqdump.output.fastq_2,
        unpaired = rules.Checkpaired.output.unpaired
    params:
        cmd = "rsem-calculate-expression",
        sample = "{sample}",
        rsem_prefix = Results + "Counts/{taxid}/{sample}_rsem"
    output:
        rsem_genes = Results + "Counts/{taxid}/RSEM/{sample}_rsem.genes.results",
        rsem_isoforms = Results + "Counts/{taxid}/RSEM/{sample}_rsem.isoforms.results"
    params:
        rsem_idx = directory(rules.Rsem_idx.output.idx_dir)
    conda:
        "envs/rsem.yaml"
    group:
        "Quantification"
    threads:
        workflow.cores / parallel_processes
    script:
        "scripts/rsem.py"

"""
Merge count matrices and dismiss samples with zero contents which exceeds the threshold percentage.
Samples which has an zero count average which exceeds the threshold set in config, will be sorted out.
"""
rule Zero_intolerance:
    input:
        genes_count_files = expand(rules.Extract_counts.output.genes_count_file, taxid = taxids, sample = samples), ## Redirect to single sample
        transcripts_count_files = expand(rules.Extract_counts.output.transcripts_count_file, taxid = taxids, sample = samples) ## Redirect to single sample
    params:
        threshold = threshold
    output:
        genes_filtered_counts = Results + "Counts/{taxid}/genes_filtered.tsv", ## Redirect to single sample
        transcripts_filtered_counts = Results + "Counts/{taxid}/transcripts_filtered.tsv" ## Redirect to single sample
    conda:
        "envs/R.yaml"
    threads:
        1
    script:
        "scripts/zero_intolerance.R"

# rule TXIMPORT AND DESEQ2:
    # input:
        # {sample}_files
    # output:
        # {sample}.normalized

# rule concatinate counts:
    # input:
        # normalzied {sample} counts
    # output:
        # {taxid}_{seq_type}.matrix


