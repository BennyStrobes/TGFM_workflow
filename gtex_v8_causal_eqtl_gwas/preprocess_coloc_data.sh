#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)



gtex_susie_input_data_file="$1"
gtex_fusion_processed_intermediate_data="$2"
coloc_data_dir="$3"



# Multivariate susie expression TWAS using SuSiE
trait_name="blood_WHITE_COUNT"
trait_sample_size="326723"
sumstat_file=$gtex_fusion_processed_intermediate_data$trait_name"_fusion_processed_sumstats_hg38.txt"
if false; then
sbatch preprocess_coloc_data_for_single_trait.sh $gtex_susie_input_data_file $trait_name $sumstat_file $trait_sample_size $coloc_data_dir
fi

trait_name="lung_FEV1FVCzSMOKE"
trait_sample_size="274172"
sumstat_file=$gtex_fusion_processed_intermediate_data$trait_name"_fusion_processed_sumstats_hg38.txt"
if false; then
sbatch preprocess_coloc_data_for_single_trait.sh $gtex_susie_input_data_file $trait_name $sumstat_file $trait_sample_size $coloc_data_dir
fi
trait_name="body_WHRadjBMIz"
trait_sample_size="336847"
sumstat_file=$gtex_fusion_processed_intermediate_data$trait_name"_fusion_processed_sumstats_hg38.txt"
if false; then
sbatch preprocess_coloc_data_for_single_trait.sh $gtex_susie_input_data_file $trait_name $sumstat_file $trait_sample_size $coloc_data_dir
fi
trait_name="bp_DIASTOLICadjMEDz"
trait_sample_size="310831"
sumstat_file=$gtex_fusion_processed_intermediate_data$trait_name"_fusion_processed_sumstats_hg38.txt"
if false; then
sbatch preprocess_coloc_data_for_single_trait.sh $gtex_susie_input_data_file $trait_name $sumstat_file $trait_sample_size $coloc_data_dir
fi