logger::log_threshold(level = logger::DEBUG)

importMetadata <- function(json) {
  json_raw <- jsonlite::read_json(path = json, simplifyVector = FALSE)
  json_vector <- unlist(json_raw)
  json_matrix <- t(json_vector)
  as.data.frame(json_matrix)
}

makeCountMatrix <- function(json_files, count_tsvs, unit, organism) {
  if (!(unit %in% c("gene", "transcript")))
      warning("`unit` should be either 'gene' or 'transcript'!")
  
  count_matrices <- lapply(X = seq_along(count_tsvs), FUN = function(i) {
    tsv <- count_tsvs[i]
    json <- json_files[i]
    meta <- importMetadata(json = json)
    species <- meta$TAXON_ID
    
    output <- FALSE
    if (species == organism)
      output <- readr::read_tsv(
        file = tsv, col_names = c(paste0(unit, "ID"), meta$RUN_accession), col_types = "cd"
      )
    return(output)
  })
  
  if (isFALSE(count_matrices)) {
    warning(paste("Species: ", organism, " not found in metadata of all samples"))
    return(FALSE)
  }
  # Remove all empty object
  count_matrices <- count_matrices[lapply(count_matrices, length) > 0]
  
  if (identical(count_matrices, list())){
    warning(paste("No count matrices found for ", unit, " files"))
    return(FALSE)
  }
  
  plyr::join_all(dfs = count_matrices, by = paste0(unit, "ID"))
}


writeOutput <- function(matrix, unit, path) {
  if (!isFALSE(matrix)){
    logger::log_debug("Writing ", unit, " count matrices to: ", path)
    readr::write_tsv(x = matrix, file = path, col_names = TRUE)
  }
}

json_file               <- snakemake@input[["json_file"]]
gene_files              <- snakemake@input[["gene_files"]]
transcript_files        <- snakemake@input[["transcript_files"]]
organisms               <- snakemake@params[["species"]]
threads                 <- snakemake@params[["threads"]]
gene_counts_files       <- snakemake@output[["gene_counts_files"]]
transcript_counts_files <- snakemake@output[["transcript_counts_files"]]

logger::log_info("Generating and exporting Count matrices for each speices")
# parallel::mclapply(X = organisms, mc.cores = threads, FUN = function(organism) {
done <- lapply(X = organisms, FUN = function(organism) {

  logger::log_debug("Generating count matrices for species: ", organism)
  gene <- makeCountMatrix(json_files = json_file, count_tsvs = gene_files, unit = "gene", organism = organism)
  transcript <- makeCountMatrix(json_files = json_file, count_tsvs = transcript_files, unit = "transcript", organism = organism)

  writeOutput(matrix = gene, unit = "gene", path = gene_counts_files)
  writeOutput(matrix = transcript, unit = "transcript", path = transcript_counts_files)
})
