#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)



tissue_name="$1"
xt_pc_gene_list_file="$2"
gtex_expression_dir="$3"
gtex_covariate_dir="$4"
gtex_genotype_dir="$5"
gtex_fusion_weights_data_dir="$6"
gtex_fusion_weights_dir="$7"
ukbb_sumstats_hg38_dir="$8"

echo $tissue_name

# Using Tiffany's paths because I didn't feel like downloading myself
gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
gemma_path="/n/groups/price/tiffany/subpheno/gemma-0.98.1-linux-static"


tissue_gtex_fusion_weights_data_dir=$gtex_fusion_weights_data_dir$tissue_name"/"
mkdir $tissue_gtex_fusion_weights_data_dir

source ~/.bash_profile
python3 preprocess_gtex_data_for_fusion_weights_analysis.py $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $tissue_gtex_fusion_weights_data_dir

gene_summary_file=$tissue_gtex_fusion_weights_data_dir$tissue_name"_gene_summary.txt"

for chrom_num in $(seq 1 22); do 
	sh filter_gtex_genotype_data_to_variants_in_ukbb.sh $chrom_num $tissue_name $gtex_genotype_dir $ukbb_sumstats_hg38_dir $tissue_gtex_fusion_weights_data_dir
done

tissue_gtex_fusion_weights_dir=$gtex_fusion_weights_dir$tissue_name"/"
mkdir $tissue_gtex_fusion_weights_dir

module load R/3.5.1




cis_window_size="500000"
# Loop through lines (genes) of gene summary file while skipping header
sed 1d $gene_summary_file | while read gene_id chrom_num tss gene_pheno_file covariate_file; do
	echo $gene_id
	# Get stand and end position of the cis window for this gene
	p0="$(($tss - $cis_window_size))"
	p1="$(($tss + $cis_window_size))"

	# deal with issue if tss is too close to zero
	if (( $p0 < 0 )); then
		p0="0"
	fi


	# Set OUT directory for this gene
	OUT=$tissue_gtex_fusion_weights_dir$tissue_name"_"$gene_id"_1KG_only_fusion_input"

	# Run PLINK to set up gene for fusion
	plink --bfile ${tissue_gtex_fusion_weights_data_dir}${tissue_name}"_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_"${chrom_num} --keep-allele-order --pheno $gene_pheno_file --make-bed --out $OUT --keep $gene_pheno_file --chr $chrom_num --from-bp $p0 --to-bp $p1
	
	# Set FINAL_OUT directory for this gene
	FINAL_OUT=$tissue_gtex_fusion_weights_dir$tissue_name"_"$gene_id"_1KG_only_fusion_output"
	Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --covar $covariate_file --PATH_gcta $gcta_path --PATH_gemma $gemma_path --verbose 0 --save_hsq --models top1,blup,lasso

	rm $OUT*
done


source ~/.bash_profile
python3 generate_fusion_pos_file_for_single_tissue.py $gene_summary_file $tissue_gtex_fusion_weights_dir $tissue_name


