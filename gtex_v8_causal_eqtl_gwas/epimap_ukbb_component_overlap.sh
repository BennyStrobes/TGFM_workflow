#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1




global_ukbb_component_file="$1"
epimap_input_dir="$2"
epimap_ukbb_genome_wide_susie_overlap_dir="$3"
liftover_directory="$4"

if false; then
python3 liftover_epimap_from_hg19_to_hg38.py $epimap_input_dir $epimap_ukbb_genome_wide_susie_overlap_dir $liftover_directory
fi

if false; then
python3 tally_number_of_epimap_links_at_various_thresholds.py $epimap_input_dir $epimap_ukbb_genome_wide_susie_overlap_dir
fi


correlation_threshold=".75"
if false; then
python3 epimap_ukbb_component_overlap.py $global_ukbb_component_file $epimap_input_dir $epimap_ukbb_genome_wide_susie_overlap_dir $correlation_threshold
fi