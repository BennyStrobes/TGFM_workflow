#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-45:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
simulated_gene_expression_dir="$8"
simulated_learned_gene_models_dir="$9"
simulated_trait_dir="${10}"
simulated_gwas_dir="${11}"
simulated_tgfm_input_data_dir="${12}"
simulated_tgfm_results_dir="${13}"
tgfm_tissues="${14}"
simulated_best_tagging_gt_dir="${15}"
gene_type="${16}"




source /home/bes710/.bash_profile
module load R/4.0.1

date
echo $simulation_number

# TGFM parameters
init_method="best"
est_resid_var="False"

# File summarizing TGFM input
tgfm_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_susie_bootstrapped_tgfm_input_data_summary.txt"
tgfm_lasso_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_susie_lasso_tgfm_input_data_summary.txt"
# File summarizing simulated gene expression
sim_gene_expression_summary=${simulated_gene_expression_dir}${simulation_name_string}"_causal_eqtl_effect_summary.txt"
# File summarizing simulated gene-trait effect sizes
sim_gene_trait_effect_size_file=${simulated_trait_dir}${simulation_name_string}"_expression_mediated_gene_causal_effect_sizes.txt"

# Add tissue prefix to simulation_name_string
simulation_name_string=${simulation_name_string}"_"${tgfm_tissues}"_"${gene_type}



echo $simulation_name_string




echo "Part 1: Uniform PMCES"
# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type


echo "Part 1: Uniform LASSO"
# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_lasso_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_lasso_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues "lasso"


echo "Part 2: Uniform Sampler"
# Uniform (sampler)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type



echo "Part 3: iterative prior"
# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_susie_pmces_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version $tgfm_tissues



echo "Part 4: prior - sampler"
version="pmces"
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type









if false; then
echo "Part 5: default susie"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_susie_variant_only_"${ln_pi_method}
python3 run_default_susie.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi


























date