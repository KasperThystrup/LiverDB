suppressPackageStartupMessages(library(dplyr))

logger::log_threshold(level = logger::DEBUG)

importCountMatrix <- function(count_file) {
  logger::log_debug("Importing count matrix")
  suppressMessages(readr::read_tsv(file = count_file, col_names = FALSE))
}


screenZeroContents <- function(counts, threshold) {
  logger::log_debug("Screening for zero inlfation")
  zero_counts <- colSums(counts == 0) / nrow(counts)
  logger::log_debug("Zero contents:")
  print(zero_counts)
  zero_counts  > threshold
}


zeroIntolerance <- function(count_file, threshold, filtered_file){
  counts <- importCountMatrix(count_file)
  
  logger::log_debug("Identifying samples with zero contents higher than: ", threshold * 100, "%")
  idx <- screenZeroContents(counts, threshold)
  
  logger::log_debug("Removing zero inflated samples")
  counts_filtered <- removeZeroInflated(counts, idx)
  
  if (ncol(counts_filtered) < 2) {
    warning("No samples passed the threshold for gene count matrices")
    counts_filtered <- data.frame()
  }
  
  readr::write_tsv(x = counts_filtered, file = filtered_file, col_names = FALSE)
}
gene_counts_files         <- snakemake@input[["gene_counts_files"]]
transcript_counts_files   <- snakemake@input[["transcript_counts_files"]]
threshold                 <- snakemake@params[["threshold"]]
threads                   <- snakemake@params[["threads"]]
gene_filtered_files       <- snakemake@output[["gene_filtered_files"]]
transcript_filtered_files <- snakemake@output[["transcript_filtered_files"]]

files <- c(gene_counts_files, transcript_counts_files)
files_filt <- c(gene_filtered_files, transcript_filtered_files)

invisible({
  # parallel::mclapply(X = seq_along(files), mc.cores = threads, FUN = function(i, threshold = threshold) {
  lapply(X = seq_along(files), FUN = function(i) {
    count_file <- files[i]
    
    counts <- importCountMatrix(count_file)
    
    idx <- screenZeroContents(counts = counts, threshold = threshold)
    
    logger::log_debug("Removing samples with high zero contents")
    counts_filtered <- counts[, !idx]
  
    
    if (ncol(counts_filtered) < 2) {
      warning("No samples passed the threshold for gene count matrices")
      counts_filtered <- data.frame()
    }
    
    filtered_file <- files_filt[i]
    readr::write_tsv(x = counts_filtered, file = filtered_file, col_names = FALSE)
  })
})
