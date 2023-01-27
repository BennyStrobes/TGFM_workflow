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

# LDSC baseline LD Dir
ldsc_baseline_ld_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"
ldsc_baseline_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baseline_v1.2/"

# Ldscore regression code
ldsc_code_dir="/n/groups/price/ben/tools/ldsc/"

# Mod-ldscore regression code
curr_dir=`pwd`
mod_ldsc_code_dir=$curr_dir"/modified_sldsc/"

# Summary statistics
full_sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"

# hg38 sldsc weights
sldsc_h38_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# Directory containing quasi-independent ld bllocks
quasi_independent_dir="/n/groups/price/ben/quasi_independent_ld_blocks/"


##################
# Output data
##################
# Output root directory
output_root="/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/"
perm_output_root="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8/"

# Directory containing hg38 ukbb summary stats
ukbb_sumstats_hg38_dir="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/ukbb_sumstats_hg38/"


# Directory containing hg38 ukbb summary stats
quasi_independent_ld_blocks_hg38_dir=$output_root"quasi_independent_ld_blocks_hg38/"

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
gtex_susie_gene_models_dir=$perm_output_root"gtex_susie_gene_models/"

# Directory containing preprocessed TGFM data
preprocessed_tgfm_data_dir=$output_root"preprocessed_tgfm_data/"

# Directory containing preprocessed TGFM-SLDSC data
preprocessed_tgfm_sldsc_data_dir=$output_root"preprocessed_tgfm_sldsc_data/"

# Directory containing number of genes and variants
num_genes_and_variants_dir=$output_root"num_genes_and_variants/"

# Directory containing TGFM heritability estimates
tgfm_heritability_results_dir=$output_root"tgfm_heritability_results/"

# Directory containing sparse-LDSC heritability estimates
sparse_ldsc_heritability_results_dir=$output_root"sparse_ldsc_heritability_results/"

# Directory containing TGFM results
tgfm_results_dir=$output_root"tgfm_results/"

# Directory containing visualizations of TGFM h2 estimates
visualize_tgfm_h2_dir=$output_root"visualize_tgfm_heritability_estimates/"

# Directory containing standard sldsc processed data
standard_sldsc_processed_data_dir=$output_root"standard_sldsc_processed_data/"

# Directory containing standard sldsc results
standard_sldsc_results_dir=$output_root"standard_sldsc_results/"

# Sparse heritability visualization dir
visualize_sparse_h2_dir=$output_root"visualize_sparse_h2/"

tgfm_sldsc_results_dir=$output_root"tgfm_sldsc_results/"

visualize_tgfm_sldsc_dir=$output_root"visualize_tgfm_sldsc/"

visualize_tgfm_dir=$output_root"visualize_tgfm/"


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
# Liftover quasi-independent ld blocks to hg38
########################################
if false; then
sh liftover_quasi_independent_ld_blocks_to_hg38.sh $liftover_directory $quasi_independent_dir $quasi_independent_ld_blocks_hg38_dir
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
	sbatch create_pseudotissue_gene_model_input_summary_file.sh $pseudotissue_name $composit_tissue_string $gtex_processed_expression_dir $gtex_pseudotissue_gene_model_input_dir $num_jobs
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

########################################
# Preprocess data for TGFM-S-LDSC
########################################
if false; then
sbatch preprocess_data_for_tgfm_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_sldsc_data_dir
fi

# Cis heritable genes
gene_type="cis_heritable_gene"
if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	chromosome_group="odd"
	sbatch preprocess_gene_ld_scores_for_tgfm_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $pseudotissue_name $chromosome_group $gtex_susie_gene_models_dir $preprocessed_tgfm_sldsc_data_dir $gene_type
	chromosome_group="even"
	sbatch preprocess_gene_ld_scores_for_tgfm_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $pseudotissue_name $chromosome_group $gtex_susie_gene_models_dir $preprocessed_tgfm_sldsc_data_dir $gene_type
done
fi
if false; then
sbatch organize_gene_ld_scores_for_tgfm_sldsc.sh $gtex_pseudotissue_file $preprocessed_tgfm_sldsc_data_dir $gene_type $gtex_susie_gene_models_dir
fi


# Only components of genes
gene_type="component_gene"
if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	chromosome_group="odd"
	sbatch preprocess_gene_ld_scores_for_tgfm_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $pseudotissue_name $chromosome_group $gtex_susie_gene_models_dir $preprocessed_tgfm_sldsc_data_dir $gene_type
	chromosome_group="even"
	sbatch preprocess_gene_ld_scores_for_tgfm_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $pseudotissue_name $chromosome_group $gtex_susie_gene_models_dir $preprocessed_tgfm_sldsc_data_dir $gene_type
done
fi

if false; then
sbatch organize_gene_ld_scores_for_tgfm_sldsc.sh $gtex_pseudotissue_file $preprocessed_tgfm_sldsc_data_dir $gene_type $gtex_susie_gene_models_dir
fi


########################################
# Run TGFM-S-LDSC
########################################
if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
	sbatch run_tgfm_sldsc.sh $preprocessed_tgfm_sldsc_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $tgfm_sldsc_results_dir $trait_name $mod_ldsc_code_dir $quasi_independent_ld_blocks_hg38_dir
done
fi
if false; then
trait_name="body_WHRadjBMIz"
sh run_tgfm_sldsc.sh $preprocessed_tgfm_sldsc_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $tgfm_sldsc_results_dir $trait_name $mod_ldsc_code_dir $quasi_independent_ld_blocks_hg38_dir
fi
########################################
# Preprocess data for TGFM
########################################
# Number of parallel jobs
num_jobs="200"
#gene_type="component_gene"
gene_type="cis_heritable_gene"
# FIle summarizing ukkbb windows
ukkbb_window_summary_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows_and_processed_data.txt"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch preprocess_data_for_tgfm.sh $ukkbb_window_summary_file $hapmap3_snpid_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $job_number $num_jobs $gene_type $preprocessed_tgfm_sldsc_data_dir
done
fi

if false; then
job_number="0"
sh preprocess_data_for_tgfm.sh $ukkbb_window_summary_file $hapmap3_snpid_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $job_number $num_jobs $gene_type $preprocessed_tgfm_sldsc_data_dir
fi

########################################
# Get number of genes and numbr of variants used in analysis
########################################
if false; then
python3 get_number_of_genes_and_number_of_variants_used.py $ukkbb_window_summary_file $gtex_susie_gene_models_dir $gtex_pseudotissue_file $num_genes_and_variants_dir
fi

########################################
# Run TGFM
########################################
gene_type="cis_heritable_genes"
independent_trait_names_file=$ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_expr_mediated.txt"
num_jobs="5"
if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_expr_mediated.txt" | while read trait_name study_file samp_size h2; do
for job_number in $(seq 0 $(($num_jobs-1))); do
	sbatch run_tgfm.sh $trait_name $ukkbb_window_summary_file $gtex_pseudotissue_file $preprocessed_tgfm_data_dir $tgfm_sldsc_results_dir $samp_size $gene_type $tgfm_results_dir $job_number $num_jobs
done
done
fi

# Organize TGFM results across parallel runs
if false; then
sh organize_tgfm_results_across_parallel_runs.sh $tgfm_results_dir $gene_type $num_jobs $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_expr_mediated.txt"
fi

########################################
# Visualize results
########################################
if false; then
source ~/.bash_profile
module load R/3.5.1
fi
Rscript visualize_tgfm_sldsc_results.R $independent_trait_names_file $tgfm_sldsc_results_dir $preprocessed_tgfm_sldsc_data_dir $gtex_tissue_colors_file $visualize_tgfm_sldsc_dir


if false; then
Rscript visualize_tgfm_results.R $independent_trait_names_file $tgfm_sldsc_results_dir $tgfm_results_dir $preprocessed_tgfm_sldsc_data_dir $gtex_tissue_colors_file $visualize_tgfm_dir
fi






########################################
# Run standard ldsc
########################################
# Summary statistic directory
sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"

# Ldscore regression code
ldsc_code_dir="/n/groups/price/ldsc/ldsc/"

# Ldsc weights
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# LDSC baselineLD annotations (hg19)
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
ldsc_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

if false; then
sh preprocess_data_for_standard_sldsc.sh $ldsc_baseline_ld_hg19_annotation_dir $standard_sldsc_processed_data_dir
fi

if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
	trait_name_full="UKB_460K."${trait_name}
	echo $trait_name_full
	sbatch run_ldsc.sh $trait_name_full $sumstat_dir $ldsc_code_dir $ldsc_baseline_ld_hg19_annotation_dir $ldsc_weights_dir $ldsc_genotype_dir $standard_sldsc_processed_data_dir $standard_sldsc_results_dir 
done
fi




if false; then
source ~/.bash_profile
module load R/3.5.1
Rscript visualize_tgfm_sldsc_estimates.R $gtex_pseudotissue_file $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" $preprocessed_tgfm_sldsc_data_dir $tgfm_sldsc_results_dir $visualize_tgfm_sldsc_dir $tgfm_heritability_results_dir
fi

if false; then
source ~/.bash_profile
module load R/3.5.1

Rscript visualize_tgfm_heritability_estimates.R $gtex_pseudotissue_file $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" $num_genes_and_variants_dir $tgfm_heritability_results_dir $standard_sldsc_results_dir $standard_sldsc_processed_data_dir $ldsc_baseline_ld_hg19_annotation_dir $visualize_tgfm_h2_dir
fi

if false; then
Rscript visualize_sparse_tgfm_heritability_estimates.R $gtex_pseudotissue_file $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" $num_genes_and_variants_dir $tgfm_heritability_results_dir $sparse_ldsc_heritability_results_dir $visualize_sparse_h2_dir
fi




