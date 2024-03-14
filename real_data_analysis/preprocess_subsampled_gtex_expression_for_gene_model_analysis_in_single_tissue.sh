#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4GB                         # Memory total in MiB (for all cores)



original_tissue_name="$1"
tissue_name="$2"
xt_pc_gene_list_file="$3"
gtex_expression_dir="$4"
gtex_covariate_dir="$5"
gtex_genotype_dir="$6"
gtex_processed_expression_dir="$7"
gtex_processed_genotype_dir="$8"
ukbb_preprocessed_for_genome_wide_susie_dir="$9"
subsample_size="${10}"


source ~/.bash_profile

# Expression
tissue_gtex_expression_data_dir=$gtex_processed_expression_dir$tissue_name"/"
mkdir $tissue_gtex_expression_data_dir
python3 preprocess_subsampled_gtex_data_for_fusion_weights_analysis.py $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $tissue_gtex_expression_data_dir $subsample_size $original_tissue_name
# List of subsampled individuals
subsampled_individuals_file=$tissue_gtex_expression_data_dir${tissue_name}"_subsampled_individuals.txt"




# Genotype
tissue_gtex_genotype_data_dir=$gtex_processed_genotype_dir$tissue_name"/"
mkdir $tissue_gtex_genotype_data_dir
for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	sh filter_gtex_variants_to_those_in_ukbb_and_subsample.sh $chrom_num $tissue_name $gtex_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir $tissue_gtex_genotype_data_dir $subsampled_individuals_file $original_tissue_name
done










