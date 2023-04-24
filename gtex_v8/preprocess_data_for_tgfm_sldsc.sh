#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




ldsc_code_dir="$1"
hapmap3_rsid_file="$2"
ldsc_baseline_annotation_dir="$3"
ldsc_baseline_ld_annotation_dir="$4"
ref_1kg_genotype_dir="$5"
gtex_pseudotissue_file="$6"
gtex_susie_gene_models_dir="$7"
preprocessed_tgfm_sldsc_data_dir="$8"
gene_type="${9}"

chrom_num="21"

num_chrom="22"



for chrom_num in $(seq 1 $(($num_chrom))); do 
	echo $chrom_num
	if false; then
	##################################
	# Create variant annotation files
	###################################
	# input file (baselineld annotation file)
	baseline_ld_annotation_stem=${ldsc_baseline_ld_annotation_dir}baselineLD.${chrom_num}
	# output file (baselineld annotation file without any qtl maxcpp annotations)
	baseline_ld_no_qtl_annotation_stem=${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}
	# Perform filtering
	source ~/.bash_profile
	python3 remove_qtl_annotations_from_baselineld_annotation_file.py $baseline_ld_annotation_stem $baseline_ld_no_qtl_annotation_stem

	# input file (baselineld annotation file)
	baseline_ld_annotation_stem=${ldsc_baseline_ld_annotation_dir}baselineLD.${chrom_num}
	# output file (baselineld annotation file without any baselineLD functional annotations)
	baseline_ld_no_funct_annotation_stem=${preprocessed_tgfm_sldsc_data_dir}baseline_plus_LDanno.${chrom_num}
	# Perform filtering
	source ~/.bash_profile
	python3 remove_functional_annotations_from_baselineld_annotation_file.py $baseline_ld_annotation_stem $baseline_ld_no_funct_annotation_stem

	# input file (baselineld annotation file)
	baseline_ld_annotation_stem=${ldsc_baseline_ld_annotation_dir}baselineLD.${chrom_num}
	# output file (baselineld annotation file without any functional annotations)
	ld_anno_only_annotation_stem=${preprocessed_tgfm_sldsc_data_dir}LDanno_only.${chrom_num}
	# Perform filtering
	source ~/.bash_profile
	python3 remove_all_functional_annotations_from_baselineld_annotation_file.py $baseline_ld_annotation_stem $ld_anno_only_annotation_stem


	# Create baseline annotation file
	cp ${ldsc_baseline_annotation_dir}baseline.${chrom_num}.annot.gz ${preprocessed_tgfm_sldsc_data_dir}baseline_no_qtl.${chrom_num}.annot.gz



	##################################
	# Create variant level LD-scores
	###################################
	# Load in LDSC module
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12

	# Run standard S-LDSC preprocessing on baselineLD annotations
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
		--ld-wind-cm 1\
		--annot ${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}.annot\
		--out ${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}\
		--print-snps ${hapmap3_rsid_file}
	# Run standard S-LDSC preprocessing on baseline annotations
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
		--ld-wind-cm 1\
		--annot ${preprocessed_tgfm_sldsc_data_dir}baseline_no_qtl.${chrom_num}.annot.gz\
		--out ${preprocessed_tgfm_sldsc_data_dir}baseline_no_qtl.${chrom_num}\
		--print-snps ${hapmap3_rsid_file}
	# Run standard S-LDSC preprocessing on baseline_plus_LDanno
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
		--ld-wind-cm 1\
		--annot ${preprocessed_tgfm_sldsc_data_dir}baseline_plus_LDanno.${chrom_num}.annot\
		--out ${preprocessed_tgfm_sldsc_data_dir}baseline_plus_LDanno.${chrom_num}\
		--print-snps ${hapmap3_rsid_file}
	# Run standard S-LDSC preprocessing on LDanno_only
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
		--ld-wind-cm 1\
		--annot ${preprocessed_tgfm_sldsc_data_dir}LDanno_only.${chrom_num}.annot\
		--out ${preprocessed_tgfm_sldsc_data_dir}LDanno_only.${chrom_num}\
		--print-snps ${hapmap3_rsid_file}

	fi

	source ~/.bash_profile
	# Filter genotype data to just regression snps in 1KG
	plink2 --bfile ${ref_1kg_genotype_dir}"1000G.EUR.hg38."${chrom_num} --extract ${hapmap3_rsid_file} --threads 1 --make-bed --keep-allele-order --out ${preprocessed_tgfm_sldsc_data_dir}"_100G_regression_snps_only."${chrom_num}
	# Create regression snp weights
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	python ${ldsc_code_dir}ldsc.py\
		--l2\
		--bfile ${preprocessed_tgfm_sldsc_data_dir}"_100G_regression_snps_only."${chrom_num}\
		--ld-wind-cm 1\
		--out ${preprocessed_tgfm_sldsc_data_dir}"regression_weights".${chrom_num}\
		--print-snps ${hapmap3_rsid_file}

	# Delete uncessary plink file
	rm ${preprocessed_tgfm_sldsc_data_dir}"_100G_regression_snps_only."${chrom_num}*


done















##################################
# No longer used in this analysis
###################################


##################################
# Create gene-level LD-scores
###################################
if false; then
source ~/.bash_profile
variant_level_ld_score_file=${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}.l2.ldscore.gz
gene_level_sldsc_output_root=${preprocessed_tgfm_sldsc_data_dir}tissue_eqtl.${chrom_num}.
python3 create_gene_level_ld_scores.py $variant_level_ld_score_file ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num} $gtex_pseudotissue_file $gtex_susie_gene_models_dir $chrom_num $gene_type $gene_level_sldsc_output_root
fi



















##############
# OLD
##############


if false; then
source ~/.bash_profile
python3 tmp_debug_generate_ld_scores.py ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num} ${hapmap3_rsid_file}
fi




