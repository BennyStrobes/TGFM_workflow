#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)






global_ukbb_component_file="$1"
abc_input_dir="$2"
abc_ukbb_genome_wide_susie_overlap_dir="$3"
liftover_directory="$4"

if false; then
python3 liftover_abc_from_hg19_to_hg38.py $abc_input_dir $abc_ukbb_genome_wide_susie_overlap_dir $liftover_directory
fi

if false; then
python3 reorganize_abc_file_by_cell_type.py $abc_ukbb_genome_wide_susie_overlap_dir
fi

if false; then
python3 tally_number_of_abc_links_at_various_thresholds.py $abc_ukbb_genome_wide_susie_overlap_dir
fi

correlation_threshold=".015"
python3 abc_ukbb_component_overlap.py $global_ukbb_component_file $abc_ukbb_genome_wide_susie_overlap_dir $correlation_threshold


correlation_threshold=".1"
python3 abc_ukbb_component_overlap.py $global_ukbb_component_file $abc_ukbb_genome_wide_susie_overlap_dir $correlation_threshold
