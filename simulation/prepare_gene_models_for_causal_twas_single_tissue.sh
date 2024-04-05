#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=50GB                         # Memory total in MiB (for all cores)



input_causal_eqtl_effects_file="$1"
simulated_learned_gene_models_stem="$2"
genotype_bim_file="$3"
pos_file="$4"
tmp_pos_file="$5"
gene_model_output_root="$6"
gene_model_base_root="$7"
eqtl_sample_size="$8"
tissue_number="$9"


python3 prepare_gene_models_for_causal_twas_pt1_single_tissue.py $input_causal_eqtl_effects_file $simulated_learned_gene_models_stem $genotype_bim_file $tmp_pos_file $gene_model_output_root $eqtl_sample_size $tissue_number

Rscript prepare_gene_models_for_causal_twas_pt2.R $tmp_pos_file $gene_model_base_root

python3 prepare_gene_models_for_causal_twas_pt3.py $tmp_pos_file $pos_file $gene_model_base_root
