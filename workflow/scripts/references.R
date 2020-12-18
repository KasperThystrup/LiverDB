`%>%` <- dplyr::`%>%`  ## Circumventing library() calls for the pipe operator
options(readr.num_columns = 0)
logger::log_threshold(logger::INFO)

metadata_files <- snakemake@input[["metadata_files"]]
release <- snakemake@params[["release"]]
taxids <- snakemake@params[["taxids"]]
output_dir <- snakemake@params[["output_dir"]]


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
  
  logger::log_debug("Fetching links for GTF files for ", species)
  ftp <- RCurl::getURL(url = paste0(url, "/"), dirlistonly = TRUE)
  files <- strsplit(x = ftp, split = "\n")[[1]]
  
  logger::log_debug(
    "Identififying target", toupper(type), "file from:\n",
    paste(files, collapse = "\n")
  )
  
  target <- grep(pattern, x = files, value = TRUE)
  
  # Ensure that at least one gtf file is located
  if (length(target) < 1) {
    stop(
      "No ", toupper(type), " files were found.\n",
      "  Ensure that the base link is correct:\n", base
    )
  }
  
  file.path(url, target)
}

downloadEnsemblResources <- function(url_file, taxid, type, path) {

  # Define destination for target files
  destfile <- as.character(file.path(path, paste(taxid, type, "gz", sep = ".")))
  
  logger::log_debug(paste("Downloading", type, "from url:", url_file, "\nSaving file as:", destfile))
  download.file(url = url_file, destfile, method = "wget", quiet = TRUE)
  
}

# ## Development options
# metadata_files = c("Results/Metadata/GSM3960735.csv", "Results/Metadata/GSM3960736.csv",
#                    "Results/Metadata/GSM3960737.csv", "Results/Metadata/GSM4735163.csv")
# release = 102
# taxids = c(10090, 9606)
# output_dir = "Results/References"
# setwd("~/tmp/test")
# ##


species_full <- lapply(metadata_files, function(meta_file) {
  meta = readr::read_csv(meta_file)
  dplyr::select(meta, TaxID, ScientificName)
}) %>%
  do.call(what = rbind)

# Reduce table to contain unique entries
species <- dplyr::summarise(species_full,
                            TaxID = unique(TaxID),
                            
                            # Substitute spaces to underscores
                            ScientificName = unique(
                              gsub(x = ScientificName, pattern = " ", replacement = "_")
                            )
)
base <- "ftp://ftp.ensembl.org/pub"

missing <- c()
for (i in 1:nrow(species)) {
  browser()
  tax <- dplyr::pull(species[i, ], TaxID)
  name <- dplyr::pull(species[i, ], ScientificName)
  
  logger::log_debug("Detected following species:", tax, " - ", name)
  logger::log_info("Downloading GTF and FASTA files for ", name)
  if (tax %in% taxids) {
    url_gtf <- locateEnsemblResources(base, release, type = "gtf", species = name)
    downloadEnsemblResources(url_file = url_gtf, taxid = tax, type = "gtf", path = output_dir)
    
    url_fa <- locateEnsemblResources(base, release, type = "fa", species = name)
    downloadEnsemblResources(url_file = url_fa, taxid = tax, type = "fa", path = output_dir)
  } else {
    missing <- c(missing, taxid)
  }
}
if (length(missing) > 0)
  stop(paste(
    'The following TaxIDs detected in the samples metadata was not specified in the "taxid.txt" file:', missing,
    'Please update the "taxids.txt" file, or check your sample selection for unintended speciments'
  ))
