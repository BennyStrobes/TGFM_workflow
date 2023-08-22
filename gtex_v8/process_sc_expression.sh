#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=90GB     


input_h5py_file="$1"
input_h5py_pb_file="$2"
processed_sc_expression_dir="$3"
individual_info_file="$4"
genotype_data_dir="$5"
visualize_processed_sc_expression_dir="$6"
processed_genotype_dir="$7"

source ~/.bash_profile



# SUBSET to cells from EUROPEAN INDIVIDUALS with genotype data
updated_individual_info_file=$processed_sc_expression_dir"individual_info_european_rna_and_dna.txt"
python3 make_list_of_individuals_with_rna_seq_and_of_european_ancestry.py $updated_individual_info_file $individual_info_file $genotype_data_dir $input_h5py_file


















##################
# OLD (no longer used)
###################


# SUBSET to cells from EUROPEAN INDIVIDUALS with genotype data
if false; then
python3 process_sc_expression.py $input_h5py_file $input_h5py_pb_file $updated_individual_info_file $processed_sc_expression_dir
fi

if false; then
python3 generate_single_cell_clusters_at_various_resolutions.py $updated_individual_info_file $input_h5py_file $processed_sc_expression_dir
fi



if false; then
module load R/3.5.1
Rscript visualize_processed_sc_expression.R $processed_sc_expression_dir $visualize_processed_sc_expression_dir
fi
