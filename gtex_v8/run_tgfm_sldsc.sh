#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=30G                         # Memory total in MiB (for all cores)




preprocessed_tgfm_sldsc_data_dir="$1"
full_sumstat_dir="$2"
ldsc_code_dir="$3"
sldsc_h38_weights_dir="$4"
ref_1kg_genotype_dir="$5"
tgfm_sldsc_results_dir="$6"
trait_name="$7"




trait_file=$full_sumstat_dir"UKB_460K."$trait_name".sumstats"

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
fi
#source /n/groups/price/ben/environments/sldsc/bin/activate
#module load python/2.7.12
data_version="baselineLD_no_qtl_gene_adj_ld_scores"
#python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}
if false; then
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baselineLD_no_qtl_pmces_gene_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
data_version="baselineLD_no_qtl_pmces_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-delete-vals --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
source ~/.bash_profile
python3 organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}


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
