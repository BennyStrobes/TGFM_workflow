#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-38:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)


pseudotissue_name="$1"
composit_tissue_string="$2"
gtex_pseudotissue_gene_model_input_dir="$3"
gtex_processed_expression_dir="$4"
gtex_processed_genotype_dir="$5"
gtex_susie_gene_models_dir="$6"
job_number="$7"




module load R/4.0.1



echo $pseudotissue_name"_"$composit_tissue_string"_"$job_number


# Using Tiffany's paths because I didn't feel like downloading myself
gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
gemma_path="/n/groups/price/tiffany/subpheno/gemma-0.98.1-linux-static"


# File containing names of genes to loop through
gene_summary_file=$gtex_pseudotissue_gene_model_input_dir$pseudotissue_name"_gene_summary_"$job_number".txt"


# Output stem
tissue_gtex_fusion_weights_dir=$gtex_susie_gene_models_dir$pseudotissue_name"/"
mkdir $tissue_gtex_fusion_weights_dir

# Window size around tss
cis_window_size="500000"






# Loop through lines (genes) of gene summary file while skipping header
sed 1d $gene_summary_file | while read gene_id chrom_num tss composit_gene_pheno_file composit_covariate_file composit_tissue_string; do
	# Get stand and end position of the cis window for this gene
	p0="$(($tss - $cis_window_size))"
	p1="$(($tss + $cis_window_size))"

	if (( $p0 < 0 )); then
		p0="0"
	fi



	# Set OUT directory for this gene
	OUT=$tissue_gtex_fusion_weights_dir$pseudotissue_name"_"$gene_id"_1KG_only_fusion_input"

	# Split composit tissue string into an array where each element is a composit tissue
	IFS=', ' read -r -a composit_tissue_arr <<< "$composit_tissue_string"
	# Split phenotype file string into an array where each element is phenotype file in a composit tissue
	IFS=', ' read -r -a gene_pheno_files <<< "$composit_gene_pheno_file"	


	# Loop through composit tissues
	for tissue_index in "${!composit_tissue_arr[@]}"; do
		# Name of composit tissue
		composit_tissue=${composit_tissue_arr[tissue_index]}
		# Pheno (expression file) corresponding to this tissue
		gene_pheno_file=${gene_pheno_files[tissue_index]}

		# Run PLINK to set up gene, tissue to get genotype of cis snps in this tissue
		plink --bfile ${gtex_processed_genotype_dir}${composit_tissue}"/"${composit_tissue}"_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_"${chrom_num} --keep-allele-order --pheno $gene_pheno_file --make-bed --out $OUT"_"$composit_tissue --keep $gene_pheno_file --chr $chrom_num --from-bp $p0 --to-bp $p1 --allow-no-sex
	done
	

	# Set FINAL_OUT directory for this gene
	FINAL_OUT=$tissue_gtex_fusion_weights_dir$pseudotissue_name"_"$gene_id"_1KG_only_fusion_output"

	# Run SuSiE analysis to create this gene model
	Rscript create_susie_gene_model_in_a_single_pseudotissue_for_one_gene.R $pseudotissue_name $composit_tissue_string $gene_id $chrom_num $composit_covariate_file $OUT $FINAL_OUT $gcta_path $gemma_path

	# Remove unnecessary files
	rm $OUT*

done









