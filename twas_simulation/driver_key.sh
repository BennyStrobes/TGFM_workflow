#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=50GB                         # Memory total in MiB (for all cores)


#############################
# Simulation of:
####### 1. TWAS
##############################



# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/scratch3/users/j/jz286/imp_geno/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"

# LDSC baseline annotations (hg19)
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"



############################
# Output directories
############################
# Output roots (one temporary on scratch and one permanent)
temp_output_root="/n/scratch3/users/b/bes710/causal_eqtl_gwas/twas_simulation/"

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

# Directory contaiing simulated gwas results
simulated_gwas_dir=$temp_output_root"simulated_gwas/"

# Directory contaiing simulated twas results
simulated_twas_dir=$temp_output_root"simulated_twas/"

# Directory contaiing simulated twas results
simulated_organized_results_dir=$temp_output_root"simulated_organized_results/"

# Directory visualizing twas results
visualize_twas_results_dir=$temp_output_root"simulated_visualization/"


############################
# Simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# Chromosome to simulate on 
chrom_num="21"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Fraction of causal genes
fraction_causal_genes="0.5"

# Total heritability
n_non_mediated_variants_per_gene="2"



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
############################
simulated_gene_position_file=${simulated_gene_position_dir}"gene_positions_chr"${chrom_num}".txt"
if false; then
sh prepare_simulated_gene_position_list.sh $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file
fi




############################
# Run single simulation
###########################
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}
if false; then
sh run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $per_element_heritability $n_non_mediated_variants_per_gene $fraction_causal_genes $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_twas_dir
fi


if false; then
for simulation_number in $(seq 2 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $per_element_heritability $n_non_mediated_variants_per_gene $fraction_causal_genes $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_twas_dir
done
fi


############################
# Organize simulation results
###########################
simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
if false; then
sh organize_simulated_twas_results.sh $simulated_twas_dir $simulated_organized_results_dir $simulation_name_string
fi


############################
# Visualize simulation results
###########################
if false; then
source ~/.bash_profile
module load R/3.5.1
fi
if false; then
Rscript visualize_simulation_twas.R $simulation_name_string $simulated_organized_results_dir $visualize_twas_results_dir
fi
