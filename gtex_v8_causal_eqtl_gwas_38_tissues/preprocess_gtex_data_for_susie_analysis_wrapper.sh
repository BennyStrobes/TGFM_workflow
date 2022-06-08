#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)



gtex_tissue_file="$1"
xt_pc_gene_list_file="$2"
eqtl_summary_stats_dir="$3"
ref_1kg_genotype_dir="$4"
ukbb_sumstats_hg38_dir="$5"
gtex_preprocessed_for_susie_dir="$6"

if false; then
for chrom_num in $(seq 1 22); do 
	sbatch preprocess_gtex_data_for_susie_analysis.sh $chrom_num $gtex_tissue_file $xt_pc_gene_list_file $eqtl_summary_stats_dir $ref_1kg_genotype_dir $ukbb_sumstats_hg38_dir $gtex_preprocessed_for_susie_dir
done
fi

source ~/.bash_profile
python3 merge_susie_input_gene_file_across_chromosomes.py $gtex_preprocessed_for_susie_dir
