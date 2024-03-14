#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)



chrom_num="$1"
tissue_name="$2"
gtex_genotype_dir="$3"
ukbb_preprocessed_for_genome_wide_susie_dir="$4"
tissue_gtex_genotype_data_dir="$5"
subsampled_individuals_file="$6"
original_tissue_name="$7"



ukbb_gtex_formatted_bim=$tissue_gtex_genotype_data_dir"ukbb_chrom_"${chrom_num}"_snps_gtex_formatted.bim"
python3 convert_ukbb_sumstat_file_to_gtex_snp_id_format.py $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_gtex_formatted_bim $chrom_num


plink --bfile ${gtex_genotype_dir}${original_tissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num} --keep ${subsampled_individuals_file} --extract $ukbb_gtex_formatted_bim --make-bed --keep-allele-order --out ${tissue_gtex_genotype_data_dir}${tissue_name}"_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_"${chrom_num}


python3 check_genotype_ordering.py ${tissue_gtex_genotype_data_dir}${tissue_name}"_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_"${chrom_num}".fam" ${subsampled_individuals_file}