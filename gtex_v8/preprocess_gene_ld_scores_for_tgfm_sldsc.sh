#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-40:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=40GB                         # Memory total in MiB (for all cores)



ldsc_code_dir="$1"
hapmap3_rsid_file="$2"
ldsc_baseline_annotation_dir="$3"
ldsc_baseline_ld_annotation_dir="$4"
ref_1kg_genotype_dir="$5"
pseudotissue_name="$6"
chromosome_group="$7"
gtex_susie_gene_models_dir="$8"
preprocessed_tgfm_sldsc_data_dir="$9"
gene_type="${10}"


source ~/.bash_profile

echo $pseudotissue_name
echo $chromosome_group

num_chrom="22"
for chrom_num in $(seq 1 $(($num_chrom))); do 
	# Residual (used to calculate even vs odd)
	rs=`expr $chrom_num % 2`
	
	# Even
	if [ "$chromosome_group" == "even" ]; then
		if [ $rs == 0 ]; then
			echo $chrom_num
			variant_level_ld_score_file=${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}.l2.ldscore.gz
			gene_level_sldsc_output_root=${preprocessed_tgfm_sldsc_data_dir}tissue_eqtl.${chrom_num}.
			python3 create_gene_level_ld_scores.py $variant_level_ld_score_file ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num} $pseudotissue_name $gtex_susie_gene_models_dir $chrom_num $gene_type $gene_level_sldsc_output_root
		fi
	fi
	# Odd
	if [ "$chromosome_group" == "odd" ]; then
		if [ $rs == 1 ]; then
			echo $chrom_num
			variant_level_ld_score_file=${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}.l2.ldscore.gz
			gene_level_sldsc_output_root=${preprocessed_tgfm_sldsc_data_dir}tissue_eqtl.${chrom_num}.
			python3 create_gene_level_ld_scores.py $variant_level_ld_score_file ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num} $pseudotissue_name $gtex_susie_gene_models_dir $chrom_num $gene_type $gene_level_sldsc_output_root
		fi
	fi
done

