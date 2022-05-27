#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-70:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)



tissue_name="$1"
pseudotissue_name="$2"
xt_pc_gene_list_file="$3"
gtex_expression_dir="$4"
gtex_covariate_dir="$5"
gtex_genotype_dir="$6"
gtex_fusion_weights_data_dir="$7"
gtex_susie_pmces_fusion_weights_dir="$8"
gtex_preprocessed_for_susie_dir="$9"
gtex_susie_results_dir="${10}"

echo $tissue_name"_"$pseudotissue_name

# Using Tiffany's paths because I didn't feel like downloading myself
gcta_path="/n/groups/price/tiffany/subpheno/fusion_twas-master/gcta_nr_robust"
gemma_path="/n/groups/price/tiffany/subpheno/gemma-0.98.1-linux-static"


tissue_gtex_fusion_weights_data_dir=$gtex_fusion_weights_data_dir$tissue_name"/"

gene_summary_file=$tissue_gtex_fusion_weights_data_dir$tissue_name"_gene_summary.txt"


tissue_gtex_fusion_weights_dir=$gtex_susie_pmces_fusion_weights_dir$tissue_name"/"
mkdir $tissue_gtex_fusion_weights_dir



module load R/4.0.1

cis_window_size="500000"



# Loop through lines (genes) of gene summary file while skipping header
sed 1d $gene_summary_file | while read gene_id chrom_num tss gene_pheno_file covariate_file; do
	# Get stand and end position of the cis window for this gene
	p0="$(($tss - $cis_window_size))"
	p1="$(($tss + $cis_window_size))"

	if (( $p0 < 0 )); then
		p0="0"
	fi

	# Set OUT directory for this gene
	OUT=$tissue_gtex_fusion_weights_dir$tissue_name"_"$gene_id"_1KG_only_fusion_input"

	# Run PLINK to set up gene for fusion
	plink --bfile ${tissue_gtex_fusion_weights_data_dir}${tissue_name}"_GTEx_v8_genotype_EUR_overlap_1kg_and_ukbb_"${chrom_num} --keep-allele-order --pheno $gene_pheno_file --make-bed --out $OUT --keep $gene_pheno_file --chr $chrom_num --from-bp $p0 --to-bp $p1
	
	# Set FINAL_OUT directory for this gene
	FINAL_OUT=$tissue_gtex_fusion_weights_dir$tissue_name"_"$gene_id"_1KG_only_fusion_output"
	Rscript FUSION.compute_weights_from_susie_pmces.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --covar $covariate_file --PATH_gcta $gcta_path --PATH_gemma $gemma_path --crossval 0 --susievarnames ${gtex_preprocessed_for_susie_dir}${gene_id}"_variant_ids.txt" --susieres ${gtex_susie_results_dir}${gene_id}"_"${pseudotissue_name}"_no_ambiguous_variants_susie_res.RDS" --hsq_p 1.0 --verbose 0 --save_hsq 

	rm $OUT*
done
