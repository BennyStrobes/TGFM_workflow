#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-50:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)



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
simulated_two_step_tgfm_results_dir="${18}"
gene_trait_architecture="${19}"
eqtl_architecture="${20}"
eqtl_sample_size="${21}"


echo $eqtl_sample_size
echo $simulated_two_step_tgfm_results_dir



source /home/bes710/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date
echo $simulation_name_string
echo $eqtl_architecture
echo $eqtl_sample_size
echo "TWO STEP"




global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"


##################################
# TGFM 
##################################
source /home/bes710/.bash_profile
module load R/4.0.1

# TGFM parameters
init_method="best"
est_resid_var="False"
gene_type="component_gene"
n_bootstraps="100"


# File summarizing TGFM input
tgfm_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_nsamp_"${n_bootstraps}"_bootstrapped_tgfm_input_data_summary.txt"


# Add tissue prefix to simulation_name_string
#tgfm_simulation_name_string=${simulation_name_string}"_nsamp_"${n_bootstraps}"_"${gene_type}
tgfm_simulation_name_string=${simulation_name_string}

echo "Part 1: TGFM with Uniform prior and only PMCES"
date
# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_two_step_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $gene_type


echo "Part 2: TGFM tissue specific prior"
date
# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_two_step_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version $n_bootstraps



echo "Part 3: Extract tissues to run two-step fine-mapping on"
date
tgfm_iterative_prior_file=$tgfm_output_stem"_iterative_variant_gene_prior_pip_level_bootstrapped.txt"
two_step_fine_mapping_tissues_file=${simulated_two_step_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}"_two_step_tissues.txt"
python3 extract_tissues_to_run_two_step_fine_mapping_on.py $tgfm_iterative_prior_file $two_step_fine_mapping_tissues_file $gene_trait_architecture

# Loop through tissues to run two step fine-mapping on 
sed 1d $two_step_fine_mapping_tissues_file | while read tissue_name tissue_pvalue; do

	echo "Part 4: Run step 1a of two-step fine-mapping on "$tissue_name
	date
	# TGFM with uniform prior and PMCES
	ln_pi_method="uniform"
	tgfm_output_stem=${simulated_two_step_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}"_two_step_"$tissue_name
	python3 run_tgfm_pmces_two_step.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $gene_type $tissue_name


	echo "Part 5: Run step 1b of two-step fine-mapping on "$tissue_name
	date
	# TGFM iterative prior
	version="pmces"
	ln_pi_method="uniform"
	n_bootstraps="100"
	tgfm_output_stem=${simulated_two_step_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}"_two_step_"$tissue_name
	python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped_two_step.py $tgfm_input_summary_file $tgfm_output_stem $version $n_bootstraps $tissue_name


	echo "Part 6: Run step 2 of two-step fine-mapping on "$tissue_name
	date
	version="pmces"
	ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
	tgfm_output_stem=${simulated_two_step_tgfm_results_dir}${tgfm_simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_iterative_two_step_"$tissue_name
	python3 run_tgfm_sampler_two_step.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $gene_type $tissue_name
done









date