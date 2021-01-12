genes_filtered_files <- snakemake@input[["genes_filtered_files"]]
transcripts_filtered_files <- snakemake@input[["transcripts_filtered_files"]]
genes_normalized_file <- snakemake@output[["genes_normalized_file"]]
transcripts_normalized_files <- snakemake@output[["transcripts_normalized_files"]]
threads <- snakemake@threads

# genes_filtered_files <- "/disk3/Data/RNAseq/dev20.004_LiverDB/Results/Counts/9606/genes_filtered.tsv"
# transcript_filtered_files <- "/disk3/Data/RNAseq/dev20.004_LiverDB/Results/Counts/9606/transcripts_filtered.tsv"
# genes_normalized_file <- "/disk3/Data/RNAseq/dev20.004_LiverDB/Results/Counts/9606/genes_normalized.tsv"
# transcripts_normalized_files <- "/disk3/Data/RNAseq/dev20.004_LiverDB/Results/Counts/9606/transcripts_normalized.tsv"
logger::log_debug("Importing filederd gene counts")
genes_filtered <- readr::read_tsv(file = genes_filtered_files)
genes <- tibble::column_to_rownames(genes_filtered, var = "EnsemblID")

samples <- colnames(genes)
metadata <- data.frame(condition = rep(x = 1, length(samples)), row.names = samples)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = genes,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq2::estimateSizeFactors(dds)
counts_raw <- DESeq2::counts(object = dds, normalized = TRUE)

counts <- tibble::rownames_to_column(as.data.frame(counts), var = "EnsemblID")

readr::write_tsv(x = counts, file = genes_normalized_file)

###########
transcripts_filtered <- readr::read_tsv(file = transcripts_filtered_files)
transcripts <- tibble::column_to_rownames(transcripts_filtered, var = "EnsemblID")

samples <- colnames(transcripts)
metadata <- data.frame(condition = rep(x = 1, length(samples)), row.names = samples)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = transcripts,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq2::estimateSizeFactors(dds)
counts_raw <- DESeq2::counts(object = dds, normalized = TRUE)

counts <- tibble::rownames_to_column(as.data.frame(counts), var = "EnsemblID")

readr::write_tsv(x = counts, file = transcripts_normalized_file)