##################
# Input data
##################

# Directory containing summary statistics
ukbb_sumstats_hg19_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_MAF0.001_v3/"

# Hapmap3 rsids
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"

# UKBB genotype
ukbb_hg38_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####2. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory="/n/groups/price/ben/tools/liftOver_x86/"

# File containing gtex tissues to do analysis on and their sample size
gtex_pseudotissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_sample_names/pseudotissue_info.txt"
# File containing gtex tissue info
gtex_tissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_sample_names/tissue_info.txt"

# File containing gtex genes
xt_gene_list_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/cross_tissue_gene_list.txt"

# Directory containing files of eqtl summary stats
eqtl_summary_stats_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_eqtl_summary_stats/"

# GTEx expression dir
gtex_expression_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/"

# GTEx covariate dir
gtex_covariate_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_covariates/"

# GTEx genotype dir
gtex_genotype_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_genotype/"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"

# Genotype data from 1KG
ref_1kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# Gtex tissue colors file
gtex_tissue_colors_file="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/input_data/gtex_tissue_colors.txt"

# UKBB in sample LD
# Generated by Martin (though this file is temporary)
ukbb_in_sample_ld_dir="/n/scratch3/users/j/jz286/imp_geno.gdreg_ld/"
ukbb_in_sample_genotype_dir="/n/scratch3/users/j/jz286/imp_geno/"





##################
# Output data
##################
# Output root directory
output_root="/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/"

# Directory containing hg38 ukbb summary stats
ukbb_sumstats_hg38_dir="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/ukbb_sumstats_hg38/"

# Directory containing misc. items
misc_dir=$output_root"misc/"

# Directory containing UKBB sumstats for genome-wide susie
ukbb_preprocessed_for_genome_wide_susie_dir=$output_root"ukbb_preprocessed_for_genome_wide_susie/"

# Directory containing hg38 GTEx preprocessed for Susie
gtex_gene_set_dir_dir=$output_root"gtex_gene_sets/"

# Directory containing GTEx processed expression
gtex_processed_expression_dir=$output_root"gtex_processed_expression/"

# Directory containing GTEx processed expression
gtex_processed_genotype_dir=$output_root"gtex_processed_genotype/"

# Directory containing pseudotissue GTEx gene model input summary files
gtex_pseudotissue_gene_model_input_dir=$output_root"gtex_pseudotissue_gene_model_input/"

# Directory containing GTEx Susie gene models
gtex_susie_gene_models_dir=$output_root"gtex_susie_gene_models/"

# Directory containing preprocessed TGFM data
preprocessed_tgfm_data_dir=$output_root"preprocessed_tgfm_data/"

# Directory containing number of genes and variants
num_genes_and_variants_dir=$output_root"num_genes_and_variants/"

# Directory containing TGFM heritability estimates
tgfm_heritability_results_dir=$output_root"tgfm_heritability_results/"

# Directory containing visualizations of TGFM h2 estimates
visualize_tgfm_h2_dir=$output_root"visualize_tgfm_heritability_estimates/"

##################
# Analysis
##################

########################################
# Liftover UKBB summary statistics to hg38
########################################
if false; then
sbatch liftover_ukbb_summary_statistics_from_hg19_to_hg38.sh $liftover_directory $ukbb_sumstats_hg19_dir $ukbb_sumstats_hg38_dir
fi

########################################
# Convert hapmap3 rsids to snp ids
########################################
hapmap3_snpid_file=$misc_dir"w_hm3.noMHC.snp_id_list.txt"
if false; then
sh convert_hapmap3_rsids_to_snpids.sh $hapmap3_rsid_file $ukbb_hg38_genotype_dir $hapmap3_snpid_file
fi


########################################
# Preprocess data for UKBB genome-wide Susie Analysis
########################################
if false; then
sh preprocess_data_for_genome_wide_ukbb_susie_analysis.sh $ukbb_sumstats_hg38_dir $gtex_genotype_dir $ref_1kg_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_sumstats_hg19_dir $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir
fi



########################################
# Filter gene list to protein coding genes
########################################
xt_pc_gene_list_file=$gtex_gene_set_dir_dir"cross_tissue_protein_coding_gene_list.txt"
if false; then
python3 filter_gene_list_to_protein_coding_genes.py $xt_gene_list_file $xt_pc_gene_list_file $gene_annotation_file
fi



########################################
# Preprocess GTEx gene expression data for gene model analysis (in each tissue seperately)
########################################
if false; then
sed 1d $gtex_tissue_file | while read tissue_name sample_size pseudotissue_name; do
	sbatch preprocess_gtex_expression_for_gene_model_analysis_in_single_tissue.sh $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_processed_expression_dir $gtex_processed_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir
done
fi

########################################
# Create pseudotissue gene model input summary file
########################################
num_jobs="10"
if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	echo $pseudotissue_name
	sh create_pseudotissue_gene_model_input_summary_file.sh $pseudotissue_name $composit_tissue_string $gtex_processed_expression_dir $gtex_pseudotissue_gene_model_input_dir $num_jobs
done
fi



########################################
# Run GTEx gene model analysis (in each pseudo-tissue seperately)
########################################
num_jobs="10"  # Must match above
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
		sbatch create_susie_gene_model_in_a_single_pseudotissue.sh $pseudotissue_name $composit_tissue_string $gtex_pseudotissue_gene_model_input_dir $gtex_processed_expression_dir $gtex_processed_genotype_dir $gtex_susie_gene_models_dir $job_number
	done
done
fi


########################################
# Organize GTEx gene model results (create pos file)
########################################
if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	sbatch organize_susie_gene_model_results_in_a_single_pseudotissue.sh $pseudotissue_name $gtex_pseudotissue_gene_model_input_dir $gtex_susie_gene_models_dir
done
fi
pseudotissue_name="Adipose_Subcutaneous"
if false; then
sh organize_susie_gene_model_results_in_a_single_pseudotissue.sh $pseudotissue_name $gtex_pseudotissue_gene_model_input_dir $gtex_susie_gene_models_dir
fi



########################################
# Preprocess data for TGFM
########################################
# Number of parallel jobs
num_jobs="200"
# FIle summarizing ukkbb windows
ukkbb_window_summary_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows_and_processed_data.txt"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch preprocess_data_for_tgfm.sh $ukkbb_window_summary_file $hapmap3_snpid_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $job_number $num_jobs
done
fi

########################################
# Get number of genes and numbr of variants used in analysis
########################################
if false; then
python3 get_number_of_genes_and_number_of_variants_used.py $ukkbb_window_summary_file $gtex_susie_gene_models_dir $gtex_pseudotissue_file $num_genes_and_variants_dir
fi


########################################
# LDSC-style genome-wide heritability estimates
########################################
learn_intercept="fixed_intercept"
learn_intercept="learn_intercept"
if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
	sbatch ldsc_style_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $learn_intercept $tgfm_heritability_results_dir"tgfm_ldsc_style_heritability_"$trait_name"_"$learn_intercept"_"
done
fi










########################################
# RSS likelihood heritability estimates
########################################
standardize_expression_boolean="True"
trait_names=( "biochemistry_Cholesterol" "body_WHRadjBMIz" "bp_DIASTOLICadjMEDz" "body_BMIz" "blood_RED_COUNT")
if false; then
for trait_name in "${trait_names[@]}"; do
	sbatch rss_likelihood_parallel_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $tgfm_heritability_results_dir"tgfm_rss_likelihood_parallel_style_heritability_"$trait_name"_standardize_expr_"$standardize_expression_boolean"_"
done
fi
if false; then
trait_name="blood_WHITE_COUNT"
sbatch rss_likelihood_parallel_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $tgfm_heritability_results_dir"tgfm_rss_likelihood_parallel_style_heritability_"$trait_name"_standardize_expr_"$standardize_expression_boolean"_"
fi


########################################
# RSS likelihood heritability estimates
########################################
standardize_expression_boolean="False"
trait_name="blood_WHITE_COUNT"
if false; then
sbatch rss_likelihood_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $tgfm_heritability_results_dir"tgfm_rss_likelihood_style_heritability_"$trait_name"_standardize_expr_"$standardize_expression_boolean"_"
fi

#trait_names=( "biochemistry_Cholesterol" "blood_WHITE_COUNT" "body_WHRadjBMIz" "bp_DIASTOLICadjMEDz" "body_BMIz" "blood_RED_COUNT")
trait_names=( "biochemistry_Cholesterol" "body_WHRadjBMIz" "bp_DIASTOLICadjMEDz" "body_BMIz" "blood_RED_COUNT")
if false; then
for trait_name in "${trait_names[@]}"; do
	sbatch rss_likelihood_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $tgfm_heritability_results_dir"tgfm_rss_likelihood_style_heritability_"$trait_name"_standardize_expr_"$standardize_expression_boolean"_"
done
fi


########################################
# RSS likelihood SVI heritability estimates
########################################
standardize_expression_boolean="False"
trait_name="blood_WHITE_COUNT"
if false; then
sbatch rss_likelihood_svi_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $tgfm_heritability_results_dir"tgfm_rss_likelihood_svi_style_heritability_"$trait_name"_standardize_expr_"$standardize_expression_boolean"_"
fi

#trait_names=( "biochemistry_Cholesterol" "blood_WHITE_COUNT" "body_WHRadjBMIz" "bp_DIASTOLICadjMEDz" "body_BMIz" "blood_RED_COUNT")
trait_names=( "biochemistry_Cholesterol" "body_WHRadjBMIz" "bp_DIASTOLICadjMEDz" "body_BMIz" "blood_RED_COUNT")
if false; then
for trait_name in "${trait_names[@]}"; do
	sbatch rss_likelihood_svi_genome_wide_heritability_estimates_for_a_trait.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $standardize_expression_boolean $tgfm_heritability_results_dir"tgfm_rss_likelihood_svi_style_heritability_"$trait_name"_standardize_expr_"$standardize_expression_boolean"_"
done
fi








if false; then
source ~/.bash_profile
module load R/3.5.1
fi
Rscript visualize_tgfm_heritability_estimates.R $gtex_pseudotissue_file $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" $num_genes_and_variants_dir $tgfm_heritability_results_dir $visualize_tgfm_h2_dir






