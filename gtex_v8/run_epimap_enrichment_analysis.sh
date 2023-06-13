#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)




tgfm_results_dir="$1"
gtex_susie_gene_models_dir="$2"
liftover_directory="$3"
gtex_pseudotissue_file="$4"
trait_list_file="$5"
epimap_data_dir="$6"
epimap_enrichment_raw_data_dir="$7"
epimap_enrichment_dir="$8"

if false; then
source ~/.bash_profile
fi
#################################################################
# First create list of epimap samples and list of epimap tracks
#################################################################
epimap_sample_file=${epimap_enrichment_dir}"epimap_samples.txt"
epimap_track_file=${epimap_enrichment_dir}"epimap_tracks.txt"
if false; then
python3 create_epimap_sample_and_track_files.py $epimap_data_dir $epimap_sample_file $epimap_track_file
fi


#################################################################
# Liftover epimap bed files from hg19 to hg38
#################################################################
if false; then
python3 liftover_epimap_bed_files_from_hg19_to_hg38.py $epimap_track_file $epimap_data_dir $epimap_enrichment_raw_data_dir $liftover_directory
fi


#################################################################
# Extract high confidence tissue specific, trait mediating variants
#################################################################
if false; then
python3 extract_high_confidence_tissue_specific_trait_mediating_variants.py $tgfm_results_dir $gtex_susie_gene_models_dir $gtex_pseudotissue_file $trait_list_file $epimap_enrichment_dir
fi

#################################################################
# Intersect variants and bed files
#################################################################
python3 intersect_variants_with_epimap_bed_files.py $gtex_pseudotissue_file $epimap_track_file $epimap_enrichment_raw_data_dir $epimap_enrichment_dir 




