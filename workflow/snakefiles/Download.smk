

# Set input file locations
samples_files = config["samples"]
taxids_file = config["taxids"]
Results = config["output_dir"]

# Set user parameters
parallel_processes = config["parallel_processes"]
email = config["email"]
ensembl_release = config["ensembl_release"]

with open (samples_files, 'r') as sample_read:
    samples = sample_read.read().splitlines()

[x for x in samples if x.strip()]


"""
Sample file metadata are collected through Entrez-direct, by querying sample IDs against the sra database.
RunInfo data are requested and saved as .csv files
"""
rule Metadata:
    params:
        sample = "{sample}",
        email = email
    output:
        metadata_file = Results + "Metadata/{sample}.csv"
    conda:
        "../envs/entrez_direct.yaml"
    group:
        "Metadata"
    threads:
        workflow.cores # Simultanous jobs causes issues
    shell:
        "esearch -db sra -query {params.sample} -email {params.email} | efetch -format runinfo > {output.metadata_file}"


"""
As the user must provide taxids in the taxids.txt file, there is a chance of human errors.
Sample taxids and Scientific namesare collected and compared against the user provided taxids.
As output, a tab-separated file containing taxids, along with the scientific names are provided.
"""
rule Checkspecies:
    input:
        metadata_files = expand(rules.Metadata.output.metadata_file, sample = samples)
    params:
        taxids = taxids
    output:
        detected_species = Results + "Metadata/detectedSpecies.txt"
    group:
        "Genome"
    threads:
        1
    script:
        "../scripts/checkspecies.R"


"""
Using the Scientific names provided from the Checkspecies rule, Ensembl gtf and fasta files are downloaded for each species.
"""
rule References:
    input:
        detected_species = rules.Checkspecies.output.detected_species
    params:
        release = ensembl_release,
    output:
        gtf_file = Results + "References/{taxid}.gtf.gz",   # Remove Expansion
        fa_file = Results + "References/{taxid}.fa.gz"   # Remove Expansion
    conda:
        "../envs/references.yaml"
    group:
        "Genome"
    threads:
        workflow.cores / parallel_processes  ## No multi-threading support
    script:
        "../scripts/references.R"


"""
Reference files are decompressed using pigz, in order to enable index generation with STAR and possibly also RSEM
"""
rule Decompress:
    input:
        gtf_file = rules.References.output.gtf_file,
        fa_file = rules.References.output.fa_file
    output:
        gtf_file = temp(Results + "References/{taxid}.gtf"),
        fa_file = temp(Results + "References/{taxid}.fa")
    wildcard_constraints:
        taxid = "\d+"
    conda:
        "../envs/decompress.yaml"
    group:
        "Genome"
    threads:
        workflow.cores / parallel_processes
    shell:
        "unpigz --force --processes {workflow.cores} {input.gtf_file} {input.fa_file}"


"""
Genome and transcript indices are generated for enabling downstream gene and transcript quantification.
"""
rule Rsem_idx:
    input:
        gtf_file = rules.Decompress.output.gtf_file,
        fa_file = rules.Decompress.output.fa_file      
    output:
        idx_dir = directory(Results + "References/{taxid}/Indices/rsem/"),
        idx_tx = Results + "References/{taxid}/Indices/rsem/rsem.transcripts.fa",
        idx_fa = Results + "References/{taxid}/Indices/rsem/rsem.idx.fa",
        idx_n2g = Results + "References/{taxid}/Indices/rsem/rsem.n2g.idx.fa"        
    conda:
        "../envs/rsem.yaml"
    group:
        "Genome"
    threads:
        workflow.cores
    shell:
        "rsem-prepare-reference \
            --gtf {input.gtf_file} \
            --star \
            --num-threads {threads} \
            {input.fa_file} \
            {output.idx_dir}/rsem"


"""
Sample metadata are used to optain SRR numbers for each sample. SRA-tools prefetch are used to download .sra files.
"""
rule Srafetch:
    input:
        metadata_file = rules.Metadata.output.metadata_file
    output:
        sra_file = temp(Results + "Rawdata/SRA/{sample}.sra")
    params:
        cmd = "prefetch"
    conda:
        "../envs/sra_tools.yaml"
    group:
        "Samples"
    threads:
        workflow.cores
    script:
        "../scripts/prefetch.py"


"""
SRA-tools Fasterq-dump are used to extract Read information into paired-end mate fastq files, enabling multithreaded processing.
"""
rule Fasterqdump:
    input:
        metadata_file = rules.Metadata.output.metadata_file,
        sra_file = rules.Srafetch.output.sra_file
    output:
        fastq_1 = Results + "Rawdata/Fastq/{sample}_1.fastq",
        fastq_2 = Results + "Rawdata/Fastq/{sample}_2.fastq"
    params:
        cmd = "fasterq-dump"
    conda:
        "../envs/sra_tools.yaml"
    group:
        "Samples"
    threads:
        min(workflow.cores / parallel_processes, 10) # Simultanous jobs may cause issues
    script:
        "../scripts/fasterq_dump.py"

rule Checkpaired:
    input:
        all_fastq_1 = expand(rules.Fasterqdump.output.fastq_1, sample = samples),
        all_fastq_2 = expand(rules.Fasterqdump.output.fastq_2, sample = samples)
    output:
        unpaired = temp(Results + "Rawdata/Fastq/unpaired.txt")
    run:
        for fastq_1 in input.all_fastq_1:

            fastq_dir = dirname(fastq_1)
            fastq = fastq_dir + basename(fastq_1).replace("_1.fastq", ".fastq")
            fastq_2 = fastq_dir + basename(fastq_1).replace("_1.fastq", "_2.fastq")
            msg = ""
            sample_id = basename(sra_file).replace("\.sra", "")

            if exists(fastq):
                print("Warning, unintended file has been detected and will be removed: %s" %fastq)
                shell("rm %s" %fastq)
            if layout == "PAIRED" and (os.stat(fastq_1).st_size == 0):
                msg += "WARNING: Read mate 1 of Paried end sample: %s is empty!\n" %sample_id
            if layout == "PAIRED" and (os.stat(fastq_2).st_size == 0):
                msg += "WARNING: Read mate 2 of Paried end sample: %s is empty!\n" %sample_id

            if msg != "":
                print(msg, "Sample will be omitted and recorded in %s" %unpaired)
                omit = "echo %s >> %s " %(sample_id, unpaired)
                shell(omit)
            else:
                shell("touch %s" %unpaired)