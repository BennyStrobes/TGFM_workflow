#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40GB                         # Memory total in MiB (for all cores)



hapmap3_rsid_file="$1"
ldsc_baseline_annotation_dir="$2"
ldsc_baseline_ld_annotation_dir="$3"
ref_1kg_genotype_dir="$4"
chrom_num="$5"
preprocessed_tgfm_sldsc_data_dir="$6"




source ~/.bash_profile

echo $chrom_num

python3 create_variant_ld_scores_for_tgfm_sldsc.py $hapmap3_rsid_file $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $chrom_num $preprocessed_tgfm_sldsc_data_dir
