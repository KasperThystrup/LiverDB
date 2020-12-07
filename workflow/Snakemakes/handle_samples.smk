# Setup run environment
from math import ceil

# Define output parameters
runinfo_files = output_dir + "Metadata/{sample}_runinfo.csv"
metadata_files = output_dir + "Metadata/{sample}.json"
species_file = output_dir + "Species.txt"

rule all:
	input:
		species_file


rule runinfo:
	params:
		email = email,
		gsm_id = "{sample}"
	output:
		runinfo_files = runinfo_files
	conda:
		"envs/entrez_direct.yaml"
	threads:
		max(ceil(workflow.cores / 3), 3)  ## At most 3 request pr second, due to Entrez-Eutils restrictions
	shell:
		"echo {params.gsm_id}; esearch -db sra -query {params.gsm_id} | efetch -format runinfo > {output.runinfo_files} -email {params.email}"

rule metadata:
	input:
		runinfo_files = expand(runinfo_files, sample = samples)
	output:
		metadata_files = expand(metadata_files, sample = samples),
		species_file = species_file
	conda:
		"envs/metadata.yaml"
	threads:
		1
	script:
		"scripts/metadata.R"
		

