#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=40GB                         # Memory total in MiB (for all cores)



chrom_num="$1"
trait_name="$2"
ukbb_genome_wide_susie_organized_results_dir="$3"
gtex_tissue_file="$4"
gtex_fusion_weights_dir="$5"
kg_genotype_dir="$6"
gtex_fusion_associations_dir="$7"
samp_size="$8"
gtex_fusion_multivariate_associations_dir="$9"


echo $trait_name
echo $chrom_num

module load R/4.0.1

if false; then
Rscript run_multivariate_twas_from_summary_stats_with_susie_for_single_trait_and_chromosome.R $chrom_num $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_tissue_file $gtex_fusion_weights_dir $kg_genotype_dir $gtex_fusion_associations_dir $samp_size $gtex_fusion_multivariate_associations_dir
fi
Rscript run_multivariate_twas_from_summary_stats_with_susie_for_single_trait_and_chromosome_debug.R $chrom_num $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_tissue_file $gtex_fusion_weights_dir $kg_genotype_dir $gtex_fusion_associations_dir $samp_size $gtex_fusion_multivariate_associations_dir
