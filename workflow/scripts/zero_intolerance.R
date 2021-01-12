# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tibble)
# })
`%>%` <- dplyr::`%>%`
options(readr.num_columns = 0)

logger::log_threshold(level = logger::DEBUG)

genes_count_files <- snakemake@input[["genes_count_files"]]
transcripts_count_files <- snakemake@input[["transcripts_count_files"]]
threshold <- snakemake@params[["threshold"]]
genes_filtered_files <- snakemake@output[["genes_filtered_files"]]
transcripts_filtered_files <- snakemake@output[["transcripts_filtered_files"]]
threads <- snakemake@threads

logger::log_debug("Importing sample gene count files")
gene_counts_raw <- lapply(genes_count_files, function(gene_file) {
  name_raw <- basename(gene_file)
  name <- gsub(name_raw, pattern = "_\\w+\\.tsv", replacement = "", )
  genes_raw <- readr::read_tsv(
    file = gene_file
  )
  genes <- tibble::column_to_rownames(genes_raw, "gene_id")
  
  colnames(genes) <- name
  
  return(genes)
}) %>% do.call(what = cbind)

logger::log_debug("Removing samples with a zero content higher than the threshold")
idx <- colMeans(gene_counts_raw == 0) <= threshold
gene_counts <- gene_counts_raw[idx]

dismissed <- colnames(gene_counts_raw[!idx])
if (length(dismissed) > 0) {
  dismissed_file <- file.path(dirname(genes_count_files), "dismissed_sample_genes.txt")
  logger::log_info(paste("Dismissed samples written to:", dismissed_file))
  readr::write_lines(x = dismissed, file = dismissed_file)
} else {
  logger::log_info(paste0("No samples had a zero-inflation above: ", threshold*100, "%"))
}

logger::log_debug("Writing filtered gene count matrix")
readr::write_tsv(x = tibble::rownames_to_column(gene_counts, "EnsemblID"),
                 file = genes_filtered_files)


logger::log_debug("Importing sample transcript count files")
transcript_counts_raw <- lapply(transcripts_count_files, function(transcript_file) {
  name_raw <- basename(transcript_file)
  name <- gsub(name_raw, pattern = "_\\w+\\.tsv", replacement = "")
  transcripts_raw <- readr::read_tsv(
    file = transcript_file
  )
  transcripts <- tibble::column_to_rownames(transcripts_raw, "transcript_id")
  colnames(transcripts) <- name
  
  return(transcripts)
}) %>% do.call(what = cbind)

logger::log_debug("Removing samples with a zero content higher than the threshold")
idx <- colMeans(transcript_counts_raw == 0) <= threshold
transcript_counts <- transcript_counts_raw[idx]

dismissed <- colnames(transcript_counts_raw[!idx])
if (length(dismissed) > 0) {
  dismissed_file <- file.path(dirname(transcripts_count_files), "dismissed_sample_transcripts.txt")
  logger::log_info(paste("Dismissed samples written to:", dismissed_file))
  readr::write_lines(x = dismissed, file = dismissed_file)
} else {
  logger::log_info(paste0("No samples had a zero-inflation above: ", threshold*100, "%"))
}

logger::log_debug("Writing filtered transcript count matrix")
readr::write_tsv(x = tibble::rownames_to_column(transcript_counts, "EnsemblID"),
                 file = transcripts_filtered_files)
