#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
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
ldsc_weights_dir="${12}"
simulated_ld_scores_dir="${13}"
mod_ldsc_code_dir="${14}"
simulated_sldsc_results_dir="${15}"
simulated_tgfm_input_data_dir="${16}"
simulated_tgfm_results_dir="${17}"
eqtl_sample_size="${18}"
parr_version="${19}"



source ~/.bash_profile
module load R/4.0.1

date
echo $simulation_number
echo $eqtl_sample_size
echo $parr_version

# TGFM parameters
init_method="best"
est_resid_var="False"

# File summarizing TGFM input
tgfm_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_bootstrapped_tgfm_input_data_summary.txt"



if [ $parr_version == "parallel_1" ]; then

echo "Part 1"
# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method



# Uniform (sampler)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi



if [ $parr_version == "parallel_2" ]; then

# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
if false; then
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version
fi

ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi


if [ $parr_version == "parallel_3" ]; then

# Uniform (sampler)
ln_pi_method="uniform"
init_method="standard"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_bayesian_"${init_method}"_"${ln_pi_method}
python3 run_tgfm_bayesian.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method


fi






























if [ $parr_version == "parallel_8" ]; then

echo "PART 2"

# variant-gene (PMCES)
# Extract variant-gene log prior information
ln_pi_method="tglr_variant_gene"
tmp_anno="pmces"
python3 extract_tglr_variant_gene_log_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $tmp_anno
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$tmp_anno


# variant-gene (sampler)
# Extract variant-gene log prior information
ln_pi_method="tglr_variant_gene"
tmp_anno="sampler"
python3 extract_tglr_variant_gene_log_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $tmp_anno
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$tmp_anno




# Iterative prior (PMCES)
version="pmces"
ln_pi_method="tglr_variant_gene"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior.py $tgfm_input_summary_file $tgfm_output_stem $version
iterative_prior_summary_file=$tgfm_output_stem"_iterative_emperical_distribution_variant_gene_prior.txt"
python3 extract_iterative_variant_gene_tissue_log_prior_info_for_tgfm.py $tgfm_input_summary_file $iterative_prior_summary_file $version
ln_pi_method="iterative_variant_gene_tissue"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$version



# Iterative prior (sampler)
version="sampler"
ln_pi_method="tglr_variant_gene"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior.py $tgfm_input_summary_file $tgfm_output_stem $version
iterative_prior_summary_file=$tgfm_output_stem"_iterative_emperical_distribution_variant_gene_prior.txt"
python3 extract_iterative_variant_gene_tissue_log_prior_info_for_tgfm.py $tgfm_input_summary_file $iterative_prior_summary_file $version
ln_pi_method="iterative_variant_gene_tissue"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$version


fi


if [ $parr_version == "parallel_20" ]; then

# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version



# Iterative prior (Sampler)
version="sampler"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version
ln_pi_method="iterative_variant_gene_tissue_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$version

fi

if [ $parr_version == "parallel_25" ]; then

# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
if false; then
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version
fi
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method


# Iterative prior (Sampler)
version="sampler"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
if false; then
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version
fi
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi

if [ $parr_version == "parallel_5" ]; then
bs_nn_tglr_file=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_tglr_bootstrapped_nonnegative_per_element_h2s.txt"
python3 extract_bootstrapped_nonnegative_tglr_variant_gene_tissue_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $bs_nn_tglr_file


ln_pi_method="tglr_bootstrapped_nonnegative_sampler"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi

if [ $parr_version == "parallel_6" ]; then
# NOTE: REQUIRES PRECOMPUTED ITERATIVE PRIOR AND TGLR PRIOR

ln_pi_method="tglr_bootstrapped_nonnegative_pmces"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method


ln_pi_method="pmces_uniform_iterative_variant_gene_prior_pip_level_pmces"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method

fi




date