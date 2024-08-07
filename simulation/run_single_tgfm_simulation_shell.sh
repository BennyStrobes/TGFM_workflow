#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-25:00                         # Runtime in D-HH:MM format
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
eqtl_sample_size="${14}"
tgfm_tissues="${15}"
simulated_best_tagging_gt_dir="${16}"
gene_type="${17}"




source /home/bes710/.bash_profile
module load R/4.0.1

date
echo $simulation_number
echo $eqtl_sample_size

# TGFM parameters
init_method="best"
est_resid_var="False"

# File summarizing TGFM input
tgfm_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_bootstrapped_tgfm_input_data_summary.txt"
# File summarizing simulated gene expression
sim_gene_expression_summary=${simulated_gene_expression_dir}${simulation_name_string}"_causal_eqtl_effect_summary.txt"
# File summarizing simulated gene-trait effect sizes
sim_gene_trait_effect_size_file=${simulated_trait_dir}${simulation_name_string}"_expression_mediated_gene_causal_effect_sizes.txt"

# Add tissue prefix to simulation_name_string
simulation_name_string=${simulation_name_string}"_"${tgfm_tissues}"_"${gene_type}



echo $simulation_name_string


echo "Part 0: Get best tagging gene tissue pairs"
best_tagging_gt_output_stem=${simulated_best_tagging_gt_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_best_tagging_gt_pairs"
if false; then
python3 get_best_tagging_gene_tissue_pairs.py $tgfm_input_summary_file $tgfm_tissues $best_tagging_gt_output_stem $sim_gene_expression_summary $processed_genotype_data_dir $sim_gene_trait_effect_size_file ${eqtl_sample_size}
fi

echo "Part 1: Uniform PMCES"
# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
if false; then
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type
fi

if false; then
echo "Part 2: Uniform Sampler"
# Uniform (sampler)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type
fi

echo "Part 3: iterative prior"
# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
if false; then
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version $tgfm_tissues
fi

echo "Part 4: prior - sampler"
version="pmces"
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
if false; then
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type
fi

echo "Part 3_v2: iterative prior"
# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
if false; then
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped_v2.py $tgfm_input_summary_file $tgfm_output_stem $version $tgfm_tissues
fi

echo "Part 4_v2: prior - sampler"
version="pmces"
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped_v3"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
if false; then
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type
fi


if false; then
echo "Part 3: iterative prior w prior"
# Iterative prior (PMCES)
version="pmces"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped_w_prior.py $tgfm_input_summary_file $tgfm_output_stem $version $tgfm_tissues


echo "Part 4: prior - sampler with prior-prior"
version="pmces"
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_w_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type
fi






if false; then
echo "Part 4: prior - sampler"
version="pmces"
ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped_cg"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method $tgfm_tissues $gene_type
fi


if false; then
echo "Part 5: default susie"
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_variant_only_"${ln_pi_method}
python3 run_default_susie.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi














################
# OLD
################






if false; then
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
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $tgfm_input_summary_file $tgfm_output_stem $version


ln_pi_method=${version}"_uniform_iterative_variant_gene_prior_pip_level_bootstrapped"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
if false; then
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi
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


fi

date