#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)




preprocessed_tgfm_sldsc_data_dir="$1"
full_sumstat_dir="$2"
ldsc_code_dir="$3"
sldsc_h38_weights_dir="$4"
ref_1kg_genotype_dir="$5"
tgfm_sldsc_results_dir="$6"
trait_name="$7"
mod_ldsc_code_dir="$8"
quasi_independent_dir="$9"




trait_file=$full_sumstat_dir"UKB_460K."$trait_name".sumstats"




variant_models=( "genotype_intercept" "baseline_no_qtl" "baselineLD_no_qtl" "baseline_plus_LDanno" "LDanno_only")

gene_types=( "component_gene")
gene_models=( "pmces_gene_adj_ld_scores")
tissue_versions=( "no_testis")
for variant_model in "${variant_models[@]}"; do
for gene_type in "${gene_types[@]}"; do
for gene_model in "${gene_models[@]}"; do
for tissue_version in "${tissue_versions[@]}"; do

	data_version=${variant_model}"_"${gene_type}"_"${tissue_version}"_"${gene_model}
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	if false; then
	python ${mod_ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${preprocessed_tgfm_sldsc_data_dir}"regression_weights." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
	fi
	# Bootstrapped version
	if false; then
	python ${mod_ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${preprocessed_tgfm_sldsc_data_dir}"regression_weights." --bootstrap --nonnegative-coefficient-file ${preprocessed_tgfm_sldsc_data_dir}${variant_model}"_"${tissue_version}"_nonnegative_coefficients.txt" --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_nonnegative_eqtl_bootstrapped_"
	fi
	source ~/.bash_profile
	python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}
done
done
done
done




















