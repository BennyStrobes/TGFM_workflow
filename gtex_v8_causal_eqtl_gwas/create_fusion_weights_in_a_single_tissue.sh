#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)



tissue_name="$1"
xt_pc_gene_list_file="$2"
gtex_expression_dir="$3"
gtex_covariate_dir="$4"
gtex_genotype_dir="$5"
gtex_fusion_weights_data_dir="$6"
gtex_fusion_weights_dir="$7"

# Using Tiffany's paths because I didn't feel like downloading myself
gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
gemma_path="/n/groups/price/tiffany/subpheno/gemma-0.98.1-linux-static"

source ~/.bash_profile
python3 preprocess_gtex_data_for_fusion_weights_analysis.py $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $gtex_fusion_weights_data_dir



gene_summary_file=$gtex_fusion_weights_data_dir$tissue_name"_gene_summary.txt"



tissue_gtex_fusion_weights_dir=$gtex_fusion_weights_dir$tissue_name"/"
module load R/3.5.1

mkdir $tissue_gtex_fusion_weights_dir


cis_window_size="500000"
# Loop through lines (genes) of gene summary file while skipping header
sed 1d $gene_summary_file | while read gene_id chrom_num tss gene_pheno_file covariate_file; do

	echo $gene_id
	# Get stand and end position of the cis window for this gene
	p0="$(($tss - $cis_window_size))"
	p1="$(($tss + $cis_window_size))"

	# Set OUT directory for this gene
	OUT=$tissue_gtex_fusion_weights_dir$tissue_name"_"$gene_id"_1KG_only_fusion_input"

	# Run PLINK to set up gene for fusion
	plink --bfile $gtex_genotype_dir$tissue_name"_GTEx_v8_genotype_EUR_"$chrom_num --pheno $gene_pheno_file --make-bed --out $OUT --keep $gene_pheno_file --chr $chrom_num --from-bp $p0 --to-bp $p1
	
	# Set FINAL_OUT directory for this gene
	FINAL_OUT=$tissue_gtex_fusion_weights_dir$tissue_name"_"$gene_id"_1KG_only_fusion_output"
	#TMP=$fusion_output_dir$data_set_name"_"$gene_id"_1KG_only_fusion_temp_output"
	#Rscript "FUSION_EDITED.compute_weight_sparsity.R" --bfile $OUT --hsq_p 0.01 --tmp $TMP --covar $covariate_file --out $FINAL_OUT --verbose 1 --crossval 0 --save_hsq --PATH_gcta ${fusion_code_dir}"gcta_nr_robust" --models lasso_fixed_lambda
	Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --covar $covariate_file --PATH_gcta $gcta_path --PATH_gemma $gemma_path --verbose 0 --save_hsq --models lasso,top1

	rm $OUT*
done


if false; then
source ~/.bash_profile
fi
if false; then
python3 generate_fusion_pos_file_for_single_tissue.py $gene_summary_file $tissue_gtex_fusion_weights_dir $tissue_name
fi

