#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
simulated_gene_expression_dir="$8"
simulated_gwas_dir="$9"
eqtl_sample_size="${10}"
ge_h2="${11}"
eqtl_architecture="${12}"
per_element_heritability="${13}"
total_heritability="${14}"
fraction_expression_mediated_heritability="${15}"
gene_trait_architecture="${16}"
simulated_trait_dir="${17}"
simulated_smr_data_dir="${18}"
simulated_smr_results_dir="${19}"


date

echo $simulation_name_string
echo $eqtl_sample_size


source /home/bes710/.bash_profile
#######################################################
# Step 1: Run GWAS on simulated trait on only snps in TGFM windows w/o standardized genotypes
#######################################################
echo "Simulation Step 1"
global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"
if false; then
python3 run_gwas_on_simulated_trait_at_snps_in_tgfm_windows_w_o_standardizing_genotype.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $global_window_file $simulated_smr_data_dir
fi


# Make learned gene models output root seperate for each simulatino
mkdir ${simulated_smr_data_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_smr_data_dir}"simulation_"${simulation_number}"/"


#######################################################
# Step 2: Get eQTL summary stats for smr
#######################################################
module load R/4.0.1
python3 get_eqtl_sumstats_for_smr.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $ge_h2 $eqtl_architecture $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $gene_trait_architecture $eqtl_sample_size




date