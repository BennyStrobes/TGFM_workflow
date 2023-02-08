#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)




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
gene_types=( "component_gene" "cis_heritable_gene")
gene_models=( "pmces_gene_adj_ld_scores")
tissue_versions=( "all_tissues" "non_sex_tissues")

if false; then
for variant_model in "${variant_models[@]}"; do
for gene_type in "${gene_types[@]}"; do
for gene_model in "${gene_models[@]}"; do
for tissue_version in "${tissue_versions[@]}"; do


	data_version=${variant_model}"_"${gene_type}"_"${tissue_version}"_"${gene_model}
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
	source ~/.bash_profile
	python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}

done
done
done
done
fi

variant_model="baseline_no_qtl"
gene_type="component_gene"
gene_model="pmces_gene_adj_ld_scores"
tissue_version="non_sex_tissues"
	data_version=${variant_model}"_"${gene_type}"_"${tissue_version}"_"${gene_model}

	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	python ${mod_ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"

#source ~/.bash_profile
#python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}






















if false; then
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="genotype_intercept_gene_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="genotype_intercept_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baselineLD_no_qtl_gene_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baselineLD_no_qtl_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baselineLD_no_qtl_pmces_gene_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}
fi








if false; then
#source /n/groups/price/ben/environments/sldsc/bin/activate
#module load python/2.7.12
data_version="baselineLD_no_qtl_pmces_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
#python ${mod_ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --w-gene-chr ${preprocessed_tgfm_sldsc_data_dir}"pmces_gene_weights." --bootstrap_window_file $quasi_independent_dir"large_10_quasi_independent_ld_blocks_hg38_bed.txt" --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_gene_weighted_"


#source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_gene_weighted_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}

#module load gcc/6.2.0
#module load python/3.6.0
#source /n/groups/price/ben/environments/tensor_flow_cpu/bin/activate
#python3 organize_tgfm_sldsc_results_tf.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}
fi















if false; then
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baseline_no_qtl_gene_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baseline_no_qtl_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baseline_no_qtl_pmces_gene_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baseline_no_qtl_pmces_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}
fi
