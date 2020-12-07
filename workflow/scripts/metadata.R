# Define script options
`%>%` <- dplyr::`%>%`  ## Circumventing library() calls for the pipe operator
logger::log_threshold(logger::INFO)
options(readr.num_columns = 0)

# Define snakemake objects
runinfo_files <- snakemake@input[["runinfo_files"]]
metadata_files <- snakemake@output[["metadata_files"]]
species_file <- snakemake@output[["species_file"]]

# ### Development variables
# gsm <- c("GSM3960735", "GSM3960736", "GSM3960737", "GSM4735163")
# 
# runinfo_files <- paste0("~/tmp/gsm/Results/Metadata/", gsm, "_runinfo.csv")
# metadata_files <- paste0("~/tmp/gsm/Results/Metadata/", gsm, ".json")
# species_file <- "~/tmp/gsm/Species.txt"
# ###

# Record species/taxonomy table
species <- lapply(seq_along(runinfo_files), function(i) {
  
  
  data_file <- runinfo_files[i]
  logger::log_debug("Reading: ", data_file)
  meta <- readr::read_csv(data_file)
  
  meta_file <- metadata_files[i]
  logger::log_debug("Writting metadata to: ", meta_file)
  dplyr::select(meta,
    "Run_accession" = Run, "TAXON_ID" = TaxID, "TUMOR" = Tumor, "LAYOUT" = LibraryLayout, "AVGLENGTH" = avgLength) %>%
    jsonlite::toJSON() %>%
    jsonlite::write_json(path = meta_file)

  # Prepare Taxonomy and Species rows
  dplyr::select(meta, TaxID, ScientificName)
}) %>% do.call(what = rbind) %>%
  
  # Reduce table to unique entries
  dplyr::summarise(
    TaxID = unique(TaxID),
    
    # Substitue spaces to underscores
    ScientificName = unique(
      gsub(x = ScientificName, pattern = " ", replacement = "_")
    )
  )

logger::log_info("Writing species table to: ", species_file)
readr::write_tsv(x = species, file = species_file, col_names = FALSE)
