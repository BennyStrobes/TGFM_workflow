#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-41:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=25GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
ldsc_real_data_results_dir="$8"
per_element_heritability="$9"
total_heritability="${10}"
fraction_expression_mediated_heritability="${11}"
ge_h2="${12}"
simulated_gene_expression_dir="${13}"
simulated_learned_gene_models_base_dir="${14}"
simulated_trait_dir="${15}"
simulated_gwas_dir="${16}"
simulated_tgfm_input_data_dir="${17}"
simulated_tgfm_results_dir="${18}"
simulated_coloc_results_dir="${19}"
gene_trait_architecture="${20}"
eqtl_architecture="${21}"
eqtl_sample_size="${22}"
simulated_focus_input_dir="${23}"
simulated_focus_results_dir="${24}"
simulated_causal_twas_input_data_dir="${25}"
simulated_causal_twas_gene_models_dir="${26}"
simulated_causal_twas_results_dir="${27}"
processed_ctwas_genotype_data_dir="${28}"
n_bootstraps="${29}"



source /home/bes710/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date
echo $simulation_name_string
echo $eqtl_architecture
echo $eqtl_sample_size


# Make learned gene models output root seperate for each simulatino
mkdir ${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}"/"


global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"


#######################################################
# Step 2: Preprocess data for TGFM
#######################################################
source /home/bes710/.bash_profile
echo "Simulation Step 2"
date
eqtl_type="susie"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $global_window_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir $n_bootstraps



##################################
# TGFM 
##################################
source /home/bes710/.bash_profile
module load R/4.0.1

# TGFM parameters
init_method="best"
est_resid_var="False"
gene_type="component_gene"

# File summarizing TGFM input
tgfm_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_nsamp_"${n_bootstraps}"_bootstrapped_tgfm_input_data_summary.txt"

# File summarizing simulated gene expression
sim_gene_expression_summary=${simulated_gene_expression_dir}${simulation_name_string}"_causal_eqtl_effect_summary.txt"
# File summarizing simulated gene-trait effect sizes
sim_gene_trait_effect_size_file=${simulated_trait_dir}${simulation_name_string}"_expression_mediated_gene_causal_effect_sizes.txt"

# Add tissue prefix to simulation_name_string
tgfm_simulation_name_string=${simulation_name_string}"_nsamp_"${n_bootstraps}"_"${gene_type}


echo "Part 5: TGFM with Uniform prior and only PMCES"
date
# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $gene_type


echo "Part 6: TGFM with Uniform prior and sampling"
date
# Uniform (sampler)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $gene_type


echo "Part 7: TGFM tissue specific prior"
date
# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version $n_bootstraps


echo "Part 8: TGFM with tissue-specific prior and sampling"
date
version="pmces"
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $gene_type








if false; then
echo "Part 10: Delete unnessary gene models"
python3 delete_unnessary_learned_gene_model_files.py $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $eqtl_sample_size
fi







date
