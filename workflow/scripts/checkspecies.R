`%>%` <- dplyr::`%>%`
options(readr.num_columns = 0)
	
metadata_files <- snakemake@input[["metadata_files"]]
taxids <- snakemake@params[["taxids"]]
detected_species <- snakemake@output[["detected_species"]]

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

taxs <- dplyr::pull(species, TaxID)
nams <- dplyr::pull(species, ScientificName)

missing <- c()
for (tax in taxids){
	predicted <- tax %in% taxs
	if (!predicted)
		missing <- c(missing, tax)
}

unexpected <- c()
for (i in 1:nrow(species)) {
  tax <- taxs[i]
  name <- nams[i]

  registered <- tax %in% taxids
  if (!registered)
  	unexpected <- c(missing, tax)
}

if (length(missing) > 0)
	stop('TaxID(s) in the "taxids.txt" file differs from the sample taxids:\n', paste(missing, collapse = " , "), '\nPlease update the "taxids.txt" file!')
if (length(unexpected) > 0)
	stop('TaxID(s) in the samples differs from the taxids file.:\n', paste(unexpected, collapse = ", "), '\nPlease update the "taxids.txt" file!')

readr::write_tsv(species, file = detected_species)	