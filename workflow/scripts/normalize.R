genes_filtered_files <- snakemake@input[["genes_filtered_files"]]
transcripts_filtered_files <- snakemake@input[["transcripts_filtered_files"]]
genes_normalized_file <- snakemake@output[["genes_normalized_file"]]
transcripts_normalized_files <- snakemake@output[["transcripts_normalized_files"]]
threads <- snakemake@threads

logger::log_debug("Importing filederd gene counts")
genes_filtered <- readr::read_tsv(file = genes_filtered_files)
transcripts_filtered <- readr::read_tsv(file = transcripts_filtered_files)
## Copy pasta code!!!
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
#readr::write_tsv(x = tibble::rownames_to_column(gene_counts, "EnsemblID"),
                 file = genes_filtered_files)