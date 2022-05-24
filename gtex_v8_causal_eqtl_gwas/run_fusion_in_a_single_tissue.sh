#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




tissue_name="$1"
gtex_fusion_weights_data_dir="$2"
gtex_fusion_weights_dir="$3"
gtex_fusion_processed_intermediate_data="$4"
gtex_fusion_associations_dir="$5"



module load R/3.5.1


if false; then
study_name="blood_WHITE_COUNT"
samp_size="326723"
for chrom_num in {1..22}; do 
	echo $study_name"_"$chrom_num
	Rscript FUSION.assoc_test.R \
	--sumstats ${gtex_fusion_processed_intermediate_data}${study_name}"_fusion_processed_sumstats_hg38.txt" \
	--weights ${gtex_fusion_weights_dir}${tissue_name}"/"${tissue_name}".pos" \
	--weights_dir "" \
	--ref_ld_chr ${gtex_fusion_processed_intermediate_data}"1000G.EUR.gtex_formatted.hg38." \
	--chr ${chrom_num} \
	--out ${gtex_fusion_associations_dir}${tissue_name}"."${study_name}"_"${chrom_num}".dat"
done


study_name="body_WHRadjBMIz"
samp_size="336847"
for chrom_num in {1..22}; do 
	echo $study_name"_"$chrom_num
	Rscript FUSION.assoc_test.R \
	--sumstats ${gtex_fusion_processed_intermediate_data}${study_name}"_fusion_processed_sumstats_hg38.txt" \
	--weights ${gtex_fusion_weights_dir}${tissue_name}"/"${tissue_name}".pos" \
	--weights_dir "" \
	--ref_ld_chr ${gtex_fusion_processed_intermediate_data}"1000G.EUR.gtex_formatted.hg38." \
	--chr ${chrom_num} \
	--out ${gtex_fusion_associations_dir}${tissue_name}"."${study_name}"_"${chrom_num}".dat"
done


study_name="bp_DIASTOLICadjMEDz"
samp_size="310831"
for chrom_num in {1..22}; do 
	echo $study_name"_"$chrom_num
	Rscript FUSION.assoc_test.R \
	--sumstats ${gtex_fusion_processed_intermediate_data}${study_name}"_fusion_processed_sumstats_hg38.txt" \
	--weights ${gtex_fusion_weights_dir}${tissue_name}"/"${tissue_name}".pos" \
	--weights_dir "" \
	--ref_ld_chr ${gtex_fusion_processed_intermediate_data}"1000G.EUR.gtex_formatted.hg38." \
	--chr ${chrom_num} \
	--out ${gtex_fusion_associations_dir}${tissue_name}"."${study_name}"_"${chrom_num}".dat"
done



study_name="lung_FEV1FVCzSMOKE"
samp_size="274172"
for chrom_num in {1..22}; do 
	echo $study_name"_"$chrom_num
	Rscript FUSION.assoc_test.R \
	--sumstats ${gtex_fusion_processed_intermediate_data}${study_name}"_fusion_processed_sumstats_hg38.txt" \
	--weights ${gtex_fusion_weights_dir}${tissue_name}"/"${tissue_name}".pos" \
	--weights_dir "" \
	--ref_ld_chr ${gtex_fusion_processed_intermediate_data}"1000G.EUR.gtex_formatted.hg38." \
	--chr ${chrom_num} \
	--out ${gtex_fusion_associations_dir}${tissue_name}"."${study_name}"_"${chrom_num}".dat"
done
fi
