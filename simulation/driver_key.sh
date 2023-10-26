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

# Mod-ldscore regression code
curr_dir=`pwd`
mod_ldsc_code_dir=$curr_dir"/modified_sldsc/"

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"


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
alkes_temp_output_root="/n/scratch3/users/a/ap92/ben/tgfm/simulation/"
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

# Directory contaiing simulated gwas results
simulated_gwas_dir=$temp_output_root"simulated_gwas/"


# Directory containing simulated coloc results
simulated_coloc_results_dir=$temp_output_root"simulated_coloc_results/"

# Directory containing organized simulation results
simulated_organized_results_dir=$perm_output_root"simulated_organized_results/"

# Directory containing simulated tgfm input data
simulated_tgfm_input_data_dir=$alkes_temp_output_root"simulated_tgfm_input/"

# Directory containing simulated tgfm results
simulated_tgfm_results_dir=$alkes_temp_output_root"simulated_tgfm_results/"

# Directory containing simulated focus input data
simulated_focus_input_dir=$temp_output_root"simulated_focus_input/"

# Directory containing simulated focus results dir
simulated_focus_results_dir=$temp_output_root"simulated_focus_results/"

# Directory containing visualizations of simulated results
visualize_simulated_results_dir=$perm_output_root"visualize_simulated_results/"



############################
# Simulation parameters
############################
# Chromosome to simulate on 
chrom_num="1"





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
# NOTE: THERE IS CURRENTLY A HACK IN HERE TO REMOVE 3 variants (out of 500000) on chrom 1 that have no variance across the 100-sample eqtl data set.
############################
n_gwas_individuals="100000"
if false; then
sbatch prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir
fi

# Need to run these 2
if false; then
n_gwas_individuals="200000"
sbatch prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir

n_gwas_individuals="50000"
sbatch prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir
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
# Main simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# cis window arround genes to define eQTLs
cis_window="100000"

# Per genetic-element heritabilities
per_element_heritability="0.0005"

# Total heritability
total_heritability="0.3"

# Fraction of heritability mediated by gene expression
fraction_expression_mediated_heritability="0.1"

# cis window arround genes to define eQTLs
cis_window="100000"

# Gene expression heritability
ge_h2="075"


############################
# Run single simulation of processing
############################
# Iteration of simulation (also works as seed)
if false; then
for simulation_number in $(seq 1 100); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir
done
fi


if false; then
n_gwas_individuals="50000"
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir

n_gwas_individuals="200000"
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir


ge_h2="05"
n_gwas_individuals="100000"
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir


ge_h2="1"
n_gwas_individuals="100000"
simulation_number="1"
simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
sbatch run_single_simulation_shell_one_ss.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir"gwas_sample_size_"${n_gwas_individuals}"/" $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $ge_h2 $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $simulated_coloc_results_dir
fi



############################
# Run single simulation of TGFM
############################
# cis window arround genes to define eQTLs
cis_window="100000"

if false; then
for simulation_number in $(seq 1 20); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}

	eqtl_sample_size="300"
	sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
	
	eqtl_sample_size="500"
	sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	eqtl_sample_size="1000"
	sh run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi


if false; then
for simulation_number in $(seq 21 100); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}

	eqtl_sample_size="300"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
	
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size

	eqtl_sample_size="1000"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi


if false; then
for simulation_number in $(seq 1 20); do 
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	eqtl_sample_size="500"
	sbatch run_single_tgfm_simulation_shell.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_tgfm_results_dir $eqtl_sample_size
done
fi





# Run FOCUS simulation
if false; then
for simulation_number in $(seq 1 20); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_focus_simulation_shell_outer.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_focus_input_dir $simulated_focus_results_dir
done
fi

if false; then
for simulation_number in $(seq 21 100); do 
	# Simulation string used for output file
	simulation_name_string="simulation_"${simulation_number}"_chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
	sbatch run_single_focus_simulation_shell_outer.sh $simulation_number $chrom_num $cis_window $n_gwas_individuals $simulation_name_string $simulated_gene_position_file $processed_genotype_data_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_gwas_dir $simulated_tgfm_input_data_dir $simulated_focus_input_dir $simulated_focus_results_dir
done
fi




# Organize simulation results across parallel simulations
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
if false; then
sh organize_simulation_results_across_parallel_simulations.sh $chrom_num $cis_window $n_gwas_individuals $global_simulation_name_string $total_heritability $fraction_expression_mediated_heritability $simulated_organized_results_dir $simulated_tgfm_results_dir $simulated_trait_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $simulated_gene_position_file $processed_genotype_data_dir $simulated_focus_results_dir $simulated_coloc_results_dir
fi


############################
# Visualize single simulation
############################
if false; then
source ~/.bash_profile
module load R/3.5.1
fi
if false; then
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}
global_simulation_name_string="chrom"${chrom_num}"_cis_window_"${cis_window}"_ss_"${n_gwas_individuals}"_ge_h2_"${ge_h2}
Rscript visualize_single_simulation.R $global_simulation_name_string $simulated_organized_results_dir $visualize_simulated_results_dir
fi


