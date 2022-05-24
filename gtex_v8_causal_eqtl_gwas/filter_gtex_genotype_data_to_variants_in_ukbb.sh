#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)





chrom_num="$1"
tissue_name="$2"
gtex_genotype_dir="$3"
ukbb_sumstats_hg38_dir="$4"
tissue_gtex_fusion_weights_data_dir="$5"



echo $chrom_num
ukbb_gtex_formatted_bim=$tissue_gtex_fusion_weights_data_dir"ukbb_chrom_"${chrom_num}"_snps_gtex_formatted.bim"
python3 convert_ukbb_sumstat_file_to_gtex_snp_id_format.py $ukbb_sumstats_hg38_dir $ukbb_gtex_formatted_bim $chrom_num


plink --bfile ${gtex_genotype_dir}${tissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num} --extract $ukbb_gtex_formatted_bim --make-bed --keep-allele-order --out ${tissue_gtex_fusion_weights_data_dir}${tissue_name}"_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_"${chrom_num}
