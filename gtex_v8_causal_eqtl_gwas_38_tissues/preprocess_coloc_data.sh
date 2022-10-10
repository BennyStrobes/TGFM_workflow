#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)



gtex_susie_input_data_file="$1"
gtex_fusion_processed_intermediate_data="$2"
coloc_data_dir="$3"
trait_file="$4"
coloc_output_dir="$5"
gtex_tissue_file="$6"

if false; then
sed 1d $trait_file | while read trait_name study_file trait_sample_size bolt_lmm_h2; do
	echo $trait_name"_"$trait_sample_size
	sumstat_file=$gtex_fusion_processed_intermediate_data$trait_name"_fusion_processed_sumstats_hg38.txt"
	sbatch preprocess_coloc_data_for_single_trait.sh $gtex_susie_input_data_file $trait_name $sumstat_file $trait_sample_size $coloc_data_dir
done
fi


total_jobs="10"
if false; then
sed 1d $trait_file | while read trait_name study_file trait_sample_size bolt_lmm_h2; do
	gene_file=$coloc_data_dir$trait_name"_processed_gene_list.txt"
	for job_number in $(seq 0 `expr $total_jobs - "1"`); do
		sbatch run_adaptive_prior_coloc_first_pass.sh $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $job_number $total_jobs
	done
done
fi
