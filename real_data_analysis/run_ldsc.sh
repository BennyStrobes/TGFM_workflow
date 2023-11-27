#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=30G                         # Memory total in MiB (for all cores)





trait_name="$1"
sumstat_dir="$2"
ldsc_code_dir="$3"
ldsc_baseline_ld_hg19_annotation_dir="$4"
ldsc_weights_dir="$5"
ldsc_genotype_dir="$6"
standard_sldsc_processed_data_dir="$7"
standard_sldsc_results="$8"



source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12


if false; then
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${standard_sldsc_processed_data_dir}"baseline_intercept_only." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.QC." --out ${standard_sldsc_results}${trait_name}"_sldsc_res_intercept_only_"
fi


trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${standard_sldsc_processed_data_dir}"baselineld_non_eqtl_cpp_removed." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.QC." --out ${standard_sldsc_results}${trait_name}"_sldsc_res_baselineld_non_eqtl_cpp_removed_"


if false; then
trait_file=$sumstat_dir$trait_name".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_hg19_annotation_dir}"baselineLD." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.QC." --out ${standard_sldsc_results}${trait_name}"_sldsc_res_"
fi









if false; then
python organize_standard_ldsc_results.py ${standard_sldsc_processed_data_dir}"baseline_intercept_only." ${standard_sldsc_results}${trait_name}"_sldsc_res_intercept_only_"
fi


if false; then
python organize_standard_ldsc_results.py ${ldsc_baseline_ld_hg19_annotation_dir}"baselineLD." ${standard_sldsc_results}${trait_name}"_sldsc_res_"
fi