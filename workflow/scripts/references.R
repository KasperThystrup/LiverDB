`%>%` <- dplyr::`%>%`  ## Circumventing library() calls for the pipe operator

species_file = snakemake@input[["species_file"]]
release = snakemake@params[["release"]]
gtf_files = snakemake@output[["gtf_files"]]
fa_files = snakemake@output[["fa_files"]]

### Development options
# species_file = "~/tmp/gsm/Results/Species.txt"
# release = 102
# gtf_files = "~/tmp/gtf_files.txt"
# fa_files = "~/tmp/fa_files.txt"
###

locateEnsemblResources <- function(base, release, type, species) {
  
  # Set up FTP location
  base_url = file.path(base, paste0("release-", release))
  
  if (type == "gtf") {
    url <- file.path(base_url, type, tolower(species))
    pattern <- paste0(species, "\\..+\\.", release, "\\.gtf\\.gz")
  } else if (type == "fasta") {
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

downloadEnsemblResources <- function(url_file, taxid, type, path, output) {

  # Define destination for target files
  destfile <- as.character(file.path(path, paste(taxid, type, "gz", sep = ".")))
  download.file(url = url_file, destfile, method = "wget", quiet = TRUE)

  # Write downloaded files
  readr::write_lines(destfile, output, append = TRUE)

}

path <- dirname(species_file)

logger::log_info("Reading: ", species_file)
species_data <- readr::read_tsv(species_file)

logger::log_debug("Species_data implemented:", species_data, sep = "\n")

species_names <- dplyr::pull(species_data, ScientificName)
taxids <- dplyr::pull(species_data, TaxID)

base <- "ftp://ftp.ensembl.org/pub"

readr::write_lines("PIK", fa_files)
for (i in seq_along(species_names)) {
  species <- species_names[i]
  tax <- taxids[i]
  logger::log_debug("Locating GTF file URL")
  url_gtf <- locateEnsemblResources(base, release, type = "gtf", species = species)
  
  logger::log_info("Downloading GTF files for ", species,
                   " writing file paths to:\n  ", gtf_files)
  downloadEnsemblResources(url_file = url_gtf, taxid = tax, type = "gtf", path, output = gtf_files)
  
  logger::log_debug("Locating FASTA file URL")
  url_fasta <- locateEnsemblResources(base, release, type = "fasta", species = species)
  
  logger::log_info("Downloading FASTA files for ", species,
                   " writing file paths to:\n  ", fa_files)
  downloadEnsemblResources(url_file = url_fasta, taxid = tax, type = "fa", path, output = fa_files)
}
