#! /bin/bash
#BSUB -o log_files/manifest_fetch-%J-output.log
#BSUB -e log_files/manifest_fetch-%J-error.log
#BSUB -q normal
#BSUB -n 4
#BSUB -G team222
#BSUB -M 2000
#BSUB -R "select[mem>2000] rusage[mem=2000]" 


# Fetch all BIOSCAN manifests in Sanger STS 
# Run from ~/aw43/bioscan_processing

echo "Starting BIOSCAN manifest fetch job"
echo "Activating virtual python environment"
module load python-3.12.0/
# python -m venv bioscan_manifests
source /lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/bioscan_manifests/bin/activate

# pip install hatchling
# pip install --extra-index-url https://pypi.org/simple -i https://gitlab.internal.sanger.ac.uk/api/v4/projects/3429/packages/pypi/simple tol-sdk
# pip install -i https://gitlab.internal.sanger.ac.uk/api/v4/projects/3429/packages/pypi/simple tol-sdk did not work... 
pip install -U --extra-index-url https://gitlab.internal.sanger.ac.uk/api/v4/projects/3429/packages/pypi/simple tol-sdk
touch .env.dev

# Generate today's date in YYYYMMDD format
TODAY=$(date +%Y%m%d)

tol data --source=portal --operation=list --type=sample \
--filter='{"and_":{"sts_project":{"eq":{"value":"BIOSCAN"}}}}' \
--fields='sts_specimen.id,sts_scientific_name,sts_col_date,sts_CATCH_SOLUTION,sts_col_time,sts_COLLECTION_METHOD,sts_habitat,sts_organism_part,sts_CATCH_LOT,sts_preservative_solution,sts_BOTTLE_DIRECTION,sts_COUNTRY_OF_COLLECTION,sts_DURATION_OF_COLLECTION,sts_collection_locality,sts_collection_country,sts_longitude,sts_latitude,bold_bin_created_date,bold_bin_uri,bold_class,bold_collectors,bold_family,bold_genus,bold_kingdom,bold_nuc,bold_nuc_basecount,bold_order,bold_phylum,bold_sequence_run_site,bold_sequence_upload_date,bold_species,bold_subfamily,bold_tribe,bioscan_extra_adult_feeding_guild,bioscan_extra_associations,bioscan_extra_broad_biotope,bioscan_extra_conservation_status,bioscan_extra_current_conservation_status,bioscan_extra_designation_summary,bioscan_extra_family,bioscan_extra_habitat,bioscan_extra_habitat_score,bioscan_extra_larval_feeding_guild,bioscan_extra_order,bioscan_extra_resources,bioscan_extra_specific_assemblage_type,bioscan_extra_sqs,bioscan_extra_vernacular,sts_AMOUNT_OF_CATCH_PLATED,sts_CATCH_BOTTLE_TEMPERATURE_STORAGE,sts_colleqct_name,sts_DATE_OF_PLATING,sts_elevation,bioscan_qc_uksi_name_status,bioscan_qc_sanger_qc_description,bioscan_qc_sanger_qc_result' \
--output=tsv > /lustre/scratch126/tol/teams/lawniczak/projects/bioscan/processing/sts_manifests_${TODAY}.tsv
# All columns included
echo "Job finished."
