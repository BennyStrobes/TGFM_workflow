#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-14:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=25GB                         # Memory total in MiB (for all cores)



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



source ~/.bash_profile
module load R/4.0.1

date
echo $simulation_number
echo $eqtl_sample_size

# TGFM parameters
init_method="best"
est_resid_var="False"

# File summarizing TGFM input
tgfm_input_summary_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_bootstrapped_tgfm_input_data_summary.txt"


# Uniform (PMCES)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
if false; then
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi

# Uniform (sampler)
ln_pi_method="uniform"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
if false; then
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method
fi

# variant-gene (PMCES)
# Extract variant-gene log prior information
ln_pi_method="tglr_variant_gene"
tmp_anno="pmces"
if false; then
python3 extract_tglr_variant_gene_log_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $tmp_anno
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$tmp_anno
fi



# variant-gene (sampler)
# Extract variant-gene log prior information
ln_pi_method="tglr_variant_gene"
tmp_anno="sampler"
if false; then
python3 extract_tglr_variant_gene_log_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $tmp_anno
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$tmp_anno
fi

# sparse-variant-gene tissue (PMCES)
ln_pi_method="tglr_sparse_variant_gene_tissue"
tmp_anno="pmces"
if false; then
python3 extract_tglr_sparse_variant_gene_tissue_log_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $tmp_anno
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$tmp_anno
fi

if false; then
# sparse-variant-gene tissue (sampler)
ln_pi_method="tglr_sparse_variant_gene_tissue"
tmp_anno="sampler"
python3 extract_tglr_sparse_variant_gene_tissue_log_prior_info_for_tgfm.py $tgfm_input_summary_file $eqtl_sample_size $simulation_name_string $simulated_sldsc_results_dir $tmp_anno
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py $tgfm_input_summary_file $tgfm_output_stem $init_method $est_resid_var $ln_pi_method"_"$tmp_anno
fi

# ITerative prior

# iterative prior (PMCES)

# Iterative prior (sampler)






date