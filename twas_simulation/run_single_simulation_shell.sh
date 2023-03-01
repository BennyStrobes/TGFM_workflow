#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=23GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
per_element_heritability="$8"
n_non_mediated_variants_per_gene="$9"
fraction_causal_genes="${10}"
simulated_gene_expression_dir="${11}"
simulated_learned_gene_models_base_dir="${12}"
simulated_trait_dir="${13}"
simulated_gwas_dir="${14}"
simulated_twas_dir="${15}"
simulated_expr_snp_corr_dir="${16}"

source ~/.bash_profile
module load R/4.0.1

echo "Simulation"$simulation_number
date


mkdir ${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}"/"

#######################################################
# Step 1: Simulate gene expression and fit gene models
#######################################################
echo "Simulation Step 1"
python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir
if false; then

#######################################################
# Step 2: Simulate trait values
#######################################################
echo "Simulation Step 2"
gene_effect_version="random"
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $per_element_heritability $n_non_mediated_variants_per_gene $fraction_causal_genes $simulated_trait_dir $n_gwas_individuals $processed_genotype_data_dir $gene_effect_version
gene_effect_version="fixed" 
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $per_element_heritability $n_non_mediated_variants_per_gene $fraction_causal_genes $simulated_trait_dir $n_gwas_individuals $processed_genotype_data_dir $gene_effect_version


#######################################################
# Step 3: Run GWAS on all snps
#######################################################
echo "Simulation Step 3"
gene_effect_version="random"
python3 run_gwas_on_simulated_trait_at_all_snps.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulated_gwas_dir $simulated_gene_expression_dir $gene_effect_version

gene_effect_version="fixed"
python3 run_gwas_on_simulated_trait_at_all_snps.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulated_gwas_dir $simulated_gene_expression_dir $gene_effect_version



#######################################################
# Step 4: Run TWAS
#######################################################
echo "Simulation Step 4"
gene_effect_version="random"
python3 run_various_twas_methods.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_twas_dir $gene_effect_version
gene_effect_version="fixed"
python3 run_various_twas_methods.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_twas_dir $gene_effect_version


#######################################################
# Step 5: Run snp correlation analysis
#######################################################
fi
echo "Simulation Step 5"
python3 run_predicted_expression_snp_correlation_analysis.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_expr_snp_corr_dir




date