#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)

trait_name="$1"
ukbb_genome_wide_susie_organized_results_dir="$2"
gtex_pseudotissue_file="$3"
pseudotissue_gtex_susie_pmces_fusion_weights_dir="$4"
gtex_fusion_processed_intermediate_data="$5"
samp_size="$6"
pseudotissue_gtex_rss_multivariate_twas_data_dir="$7"
pseudotissue_gtex_rss_multivariate_twas_dir="$8"
gene_version="$9"
coloc_results_dir="${10}"



if false; then
for chrom_num in {1..22}; do 
	sbatch organize_rss_twas_data_on_single_chromosome.sh $chrom_num $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $gene_version
done
fi

if false; then
for chrom_num in {1..22}; do 
	component_data_file=${pseudotissue_gtex_rss_multivariate_twas_data_dir}${trait_name}"_"${gene_version}"_"${chrom_num}"_component_rss_multivariate_twas_data_organized.txt"
	sbatch run_rss_twas_on_single_chromosome.sh $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version
done
fi

gene_count_method="count_genes_once"
init_version="trained_init"
if false; then
sbatch run_rss_twas_tissue_specific_prior_inference.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version
fi
gene_count_method="count_all_genes"
init_version="trained_init"
if false; then
sbatch run_rss_twas_tissue_specific_prior_inference.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version
fi

gene_count_method="count_genes_once"
init_version="null_init"
if false; then
sbatch run_rss_twas_tissue_specific_prior_inference.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version
fi
gene_count_method="count_all_genes"
init_version="null_init"
if false; then
sbatch run_rss_twas_tissue_specific_prior_inference.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version
fi


# Gonna move forward with following versions (gene_count_method="count_genes_once" and init_version="null_init")
if false; then
for chrom_num in {1..22}; do 
	component_data_file=${pseudotissue_gtex_rss_multivariate_twas_data_dir}${trait_name}"_"${gene_version}"_"${chrom_num}"_component_rss_multivariate_twas_data_organized.txt"
	sbatch run_rss_twas_on_single_chromosome_with_tissue_specific_prior.sh $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version
done
fi




#####################
# ROBUST TWAS
#####################
if false; then
for chrom_num in {1..22}; do 
	component_data_file=${pseudotissue_gtex_rss_multivariate_twas_data_dir}${trait_name}"_"${gene_version}"_"${chrom_num}"_component_rss_multivariate_twas_data_organized.txt"
	sbatch run_robust_rss_twas_on_single_chromosome.sh $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version
done
fi

gene_count_method="count_genes_once"
init_version="null_init"
if false; then
sbatch run_robust_rss_twas_tissue_specific_prior_inference.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version
fi


if false; then
for chrom_num in {1..22}; do 
	component_data_file=${pseudotissue_gtex_rss_multivariate_twas_data_dir}${trait_name}"_"${gene_version}"_"${chrom_num}"_component_rss_multivariate_twas_data_organized.txt"
	sbatch run_robust_rss_twas_on_single_chromosome_with_tissue_specific_prior.sh $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version
done
fi

if false; then
for chrom_num in {1..22}; do 
	component_data_file=${pseudotissue_gtex_rss_multivariate_twas_data_dir}${trait_name}"_"${gene_version}"_"${chrom_num}"_component_rss_multivariate_twas_data_organized.txt"
	sbatch run_robust_rss_twas_on_single_chromosome_with_pleiotropic_only_prior.sh $chrom_num $trait_name $gtex_pseudotissue_file $component_data_file $ukbb_genome_wide_susie_organized_results_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version
done
fi






gene_count_method="count_genes_once"
init_version="null_init"
if false; then
sbatch run_robust_rss_twas_pleiotropic_effect_only_prior_inference.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_count_method $init_version
fi













if false; then
sbatch organize_rss_twas_results.sh $trait_name $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $coloc_results_dir
fi