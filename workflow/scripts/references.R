`%>%` <- dplyr::`%>%`  ## Circumventing library() calls for the pipe operator
options(readr.num_columns = 0)
logger::log_threshold(logger::INFO)

detected_species <- snakemake@input[["detected_species"]]
release <- snakemake@params[["release"]]
gtf_file <- snakemake@output[["gtf_file"]]
fa_file <- snakemake@output[["fa_file"]]


locateEnsemblResources <- function(base, release, type, species) {

  # Set up FTP location
  base_url = file.path(base, paste0("release-", release))
  
  if (type == "gtf") {
    url <- file.path(base_url, type, tolower(species))
    pattern <- paste0(species, "\\..+\\.", release, "\\.gtf\\.gz")
  } else if (type == "fa") {
    url <- file.path(base_url, "fasta", tolower(species), "dna")
    pattern <- paste0(species, "\\..+\\.dna\\.primary_assembly(\\.\\w+)?\\.fa\\.gz")
  } else {
    stop("Argument: type = '", type, "' is not correctly set.\n",
         "  Use either type = `gtf` or `fasta`!")
  }
  
  logger::log_debug("Fetching links for ", species)
  ftp <- RCurl::getURL(url = paste0(url, "/"), dirlistonly = TRUE)
  files <- strsplit(x = ftp, split = "\n")[[1]]
  
  logger::log_debug(
    "Identififying target ", toupper(type), " file from:\n",
    paste(files, collapse = "\n")
  )
  
  target <- grep(pattern, x = files, value = TRUE)
  
  # Ensure that at least one gtf file is located
  if (length(target) < 1) {
    stop(
      "No ", toupper(type), " files were found at:\n", url,
      "  Ensure that the base link and the url is correct:\n", base
    )
  }
  
  file.path(url, target)
}

downloadEnsemblResources <- function(url_file, output) {

  # Define destination for target files  
  logger::log_debug(paste("Downloading from url:", url_file, "\nSaving file as:", output))
  download.file(url = url_file, output, method = "wget", quiet = TRUE)
  
}

# ## Development options
# metadata_files = c("Results/Metadata/GSM3960735.csv", "Results/Metadata/GSM3960736.csv",
#                    "Results/Metadata/GSM3960737.csv", "Results/Metadata/GSM4735163.csv")
# release = 102
# taxids = c(10090, 9606)
# output_dir = "Results/References"
# setwd("~/tmp/test")
# ##

species <- readr::read_tsv(detected_species)

base <- "ftp://ftp.ensembl.org/pub"

for (i in 1:nrow(species)) {
  tax <- dplyr::pull(species[i, ], TaxID)
  name <- dplyr::pull(species[i, ], ScientificName)

  logger::log_debug("Detected following species:", tax, " - ", name)
  logger::log_info("Downloading GTF and FASTA files for ", name)
  
  url_gtf <- locateEnsemblResources(base, release, type = "gtf", species = name)
  downloadEnsemblResources(url_file = url_gtf, output = gtf_file)

  url_fa <- locateEnsemblResources(base, release, type = "fa", species = name)
  downloadEnsemblResources(url_file = url_fa, output = fa_file)

}
