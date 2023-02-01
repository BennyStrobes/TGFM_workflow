#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=50GB                         # Memory total in MiB (for all cores)


#############################
# Simulation of:
####### 1. TGFM-SLDSC
####### 2. TGFM fine-mapping
##############################




############################
# Input data
############################

# Directory containing ldsc scripts
ldsc_code_dir="/n/groups/price/ben/tools/ldsc/"

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/scratch3/users/j/jz286/imp_geno/"

# Ldsc weights for hg19
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# LDSC baseline annotations (hg19)
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"

# LDSC baselineLD annotations (hg19)
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

# Directory containing LDSC Summary statistics (This will be used to create a realistic gene model)
ldsc_summary_stats_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"



############################
# Output directories
############################
# Output roots (one temporary on scratch and one permanent)
temp_output_root="/n/scratch3/users/b/bes710/causal_eqtl_gwas/simulation/"
perm_output_root="/n/groups/price/ben/causal_eqtl_gwas/simulation/"

# Directory containing processed genotype data
ldsc_real_data_results_dir=$temp_output_root"ldsc_real_data_results/"

# Directory containing processed genotype data
processed_genotype_data_dir=$temp_output_root"processed_genotype/"

# Directory containing simulated gene positions
simulated_gene_position_dir=$temp_output_root"simulated_gene_positions/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_gene_expression_dir=$temp_output_root"simulated_gene_expression/"

# Directory containing simulated gene expression (ie causal eQTL effect sizes in each tissue)
simulated_learned_gene_models_dir=$temp_output_root"simulated_learned_gene_models/"

# Directory containing simulated trait values
simulated_trait_dir=$temp_output_root"simulated_trait/"




############################
# Simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# Chromosome to simulate on 
chrom_num="21"

# cis window arround genes to define eQTLs
cis_window="25000"

# Per genetic-element heritabilities
per_element_heritability="0.0001"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"








############################
# Run S-LDSC on real data.
## This is done in order to have a realistic heritability model (taus)
############################
if false; then
sh run_ldsc_on_real_trait_data_to_get_realistic_heritability_model.sh $ldsc_code_dir $ldsc_baseline_hg19_annotation_dir $ldsc_weights_dir $ldsc_summary_stats_dir $kg_genotype_dir $ldsc_real_data_results_dir
fi

############################
# Prepare genotype data for analysis:
## 1. Filter number of individuals in original data
## 2. Filter sites to be those in LDSC annotation file
## 3. Convert to plink bed files
############################
if false; then
sh prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir
fi

############################
# Prepare gene file for simulation:
# Genes are defined by actual tss
# Limit to protein coding genes
# In simulation, I will assume gene is expressed in each tissue
############################
simulated_gene_position_file=${simulated_gene_position_dir}"gene_positions_chr"${chrom_num}".txt"
if false; then
sh prepare_simulated_gene_position_list.sh $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file
fi





############################
# Run single simulation
############################
# Iteration of simulation (also works as seed)
simulation_number="1" 
# Simulation string used for output file
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}
if false; then
sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir
fi













