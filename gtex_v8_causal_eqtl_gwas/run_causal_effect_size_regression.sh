#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)





trait_name="$1"
gtex_gene_file="$2"
gtex_susie_results_dir="$3"
ukbb_susie_results_dir="$4"
gtex_pseudotissue_file="$5"
gtex_onto_trait_causal_effect_size_regression_dir="$6"
gtex_onto_trait_causal_effect_size_regression_viz_dir="$7"


source ~/.bash_profile


########## NEED TO UPDATE SO OUTPUT FILE HAS TRAIT NAME IN IT


# Number of parallel jobs
total_jobs="20"

if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
	echo $job_number
	sbatch run_causal_effect_size_regression_in_parallel.sh $trait_name $gtex_gene_file $gtex_susie_results_dir $ukbb_susie_results_dir $gtex_pseudotissue_file $gtex_onto_trait_causal_effect_size_regression_dir $job_number $total_jobs
done
fi

if false; then
python3 merge_causal_effect_size_regressions.py $trait_name $gtex_onto_trait_causal_effect_size_regression_dir $total_jobs
fi

tiffany_dir="/n/groups/price/tiffany/subpheno/AllGTExTissues_restore/Marginal_alphas_sumstats_1KG_v8_320EUR/UKB_460K.blood_WHITE_COUNT/"
if false; then
python3 preprocess_tiffanys_marginal_twas_results.py $gtex_pseudotissue_file $gtex_onto_trait_causal_effect_size_regression_dir $tiffany_dir
fi

module load R/3.5.1

Rscript visualize_causal_effect_size_regression.R $trait_name $gtex_onto_trait_causal_effect_size_regression_dir $gtex_onto_trait_causal_effect_size_regression_viz_dir

