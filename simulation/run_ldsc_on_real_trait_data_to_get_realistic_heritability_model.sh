#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)



ldsc_code_dir="$1"
ldsc_baseline_hg19_annotation_dir="$2"
ldsc_weights_dir="$3"
ldsc_summary_stats_dir="$4"
genotype_dir="$5"
ldsc_real_data_results_dir="$6"




source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12



trait_name="blood_WHITE_COUNT"
trait_file=$ldsc_summary_stats_dir"UKB_460K."${trait_name}".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_hg19_annotation_dir}"baseline." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${genotype_dir}"1000G.EUR.QC." --out ${ldsc_real_data_results_dir}${trait_name}"_sldsc_baseline"