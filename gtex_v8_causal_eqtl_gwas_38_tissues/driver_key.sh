
##################
# Input data
##################

# Directory containing summary statistics
ukbb_sumstats_hg19_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_MAF0.001_v3/"

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

# UKBB download data
ukbb_download_data="/n/groups/price/UKBiobank/download_500K/"

# UKBB Genotype data
ukbb_genotype_data="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"

# UKBB Phenotype files
ukbb_pheno_file1="/n/groups/price/steven/RareVariants/Final/UKB_new_sumstats/UKB_v3.061518.tab"
ukbb_pheno_file2="/n/groups/price/UKBiobank/app19808mosaic/bloodQC/ukb4777.blood_v2.covars.tab"
ukbb_pheno_file3="/n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"

# Genotype data from 1KG
ref_1kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# GTEx expression dir
gtex_expression_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/"

# GTEx covariate dir
gtex_covariate_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_covariates/"

# GTEx genotype dir
gtex_genotype_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_genotype/"


# EpiMap enhancer-gene links
epimap_input_dir="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/input_data/epimap/per_group_enhancer_gene_links/"

# ABC enhancer gene links
abc_input_dir="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/input_data/abc/"

# Gtex tissue colors file
gtex_tissue_colors_file="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/input_data/gtex_tissue_colors.txt"

# ABC-enhancer, gwas-variant enrichment across biosamples (supp table 6 of Nassar et al)
# wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03446-x/MediaObjects/41586_2021_3446_MOESM8_ESM.txt
abc_enhancer_gwas_variant_enrichment_file="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/input_data/abc/gwas_variant_abc_enhancer_enrichments_across_biosamples_supp_table_6.txt"

# Gene set annotation file
gene_set_annotation_file="/n/groups/price/martin/data_GDREG/sam_gene_annotation/geneset_list_curated.csv"

##################
# Output data
##################
# Output root directory
output_root="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas_38_tissues/"

# Directory containing hg38 ukbb summary stats
ukbb_sumstats_hg38_dir=$output_root"ukbb_sumstats_hg38/"

# Directory containing UKBB sumstats for genome-wide susie
ukbb_preprocessed_for_genome_wide_susie_dir=$output_root"ukbb_preprocessed_for_genome_wide_susie/"

# Directory containing Susie results when applied to genome-wide UKBB
ukbb_genome_wide_susie_results_dir=$output_root"ukbb_genome_wide_susie_results/"

# Directory containing Susie organized results when applied to genome-wide UKBB
ukbb_genome_wide_susie_organized_results_dir=$output_root"ukbb_genome_wide_susie_organized_results/"

# Directory containing visualizations of tissue component mapping
tissue_component_mapping_viz=$output_root"visualize_tissue_mapping_ukbb_genome_wide_susie/"

# Directory containing hg38 GTEx preprocessed for Susie
gtex_gene_set_dir_dir=$output_root"gtex_gene_sets/"

# Directory containing GTEx fusion weights input
gtex_fusion_weights_data_dir=$output_root"gtex_fusion_weights_data/"

# Directory containing GTEx fusion weights
gtex_fusion_weights_dir=$output_root"gtex_fusion_weights/"

# Directory containing organized GTEx fusion weights
pseudotissue_gtex_fusion_weights_dir=$output_root"pseudotissue_gtex_fusion_weights/"

# Directory containing fusion processed data
gtex_fusion_processed_intermediate_data=$output_root"gtex_fusion_processed_intermediate_data/"

# Directory containing hg38 GTEx preprocessed for Susie
gtex_preprocessed_for_susie_dir=$output_root"gtex_preprocessed_for_susie/"

# Directory containing Susie results when applied to GTEx
gtex_susie_results_dir=$output_root"gtex_susie_results/"

# Directory containing GTEx Susie PMCES fusion weights
gtex_susie_pmces_fusion_weights_dir=$output_root"gtex_susie_pmces_fusion_weights/"

# Directory containing organized GTEx Susie PMCES fusion weights
pseudotissue_gtex_susie_pmces_fusion_weights_dir=$output_root"pseudotissue_gtex_susie_pmces_fusion_weights/"

# Directory containing rss multivariate twas data
pseudotissue_gtex_rss_multivariate_twas_data_dir=$output_root"pseudotissue_gtex_rss_multivariate_twas_data/"

# Directory containing rss multivariate twas results
pseudotissue_gtex_rss_multivariate_twas_dir=$output_root"pseudotissue_gtex_rss_multivariate_twas/"

# Directory containing visualizations of multivariate twas
rss_multivariate_twas_visualization_dir=$output_root"visualize_rss_multivariate_twas/"

# Directory containing visualizations of multivariate twas
rss_multivariate_twas_genome_prior_visualization_dir=$output_root"visualize_rss_multivariate_twas_genome_prior/"

# Directory containing visualizations of multivariate twas
rss_joint_multivariate_twas_visualization_dir=$output_root"visualize_rss_multivariate_twas_joint/"

# Directory containing visualizations of fine-mapped components
viz_fine_mapped_components_dir=$output_root"visualize_example_fine_mapped_components/"

# Directory containing processed coloc data
coloc_data_dir=$output_root"coloc_data/"

# Directory containing coloc results
coloc_results_dir=$output_root"coloc_results/"

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
# Run genome-wide susie analysis
########################################
# Preprocess data for UKBB genome-wide Susie Analysis
if false; then
sh preprocess_data_for_genome_wide_ukbb_susie_analysis.sh $ukbb_sumstats_hg38_dir $gtex_genotype_dir $ref_1kg_genotype_dir $ukbb_preprocessed_for_genome_wide_susie_dir
fi

# Run SuSiE on UKBB genome-wide data (windows)
window_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows_and_processed_data.txt"
if false; then
sh run_susie_genome_wide_shell.sh $window_file $ukbb_genome_wide_susie_results_dir
fi

# Organize results of SuSiE run genome-wide 
window_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows_and_processed_data.txt"
if false; then
sbatch organize_susie_genome_wide_results.sh $window_file $ukbb_genome_wide_susie_results_dir $ukbb_genome_wide_susie_organized_results_dir
fi







########################################
# Filter gene list to protein coding genes
########################################
xt_pc_gene_list_file=$gtex_gene_set_dir_dir"cross_tissue_protein_coding_gene_list.txt"
if false; then
python3 filter_gene_list_to_protein_coding_genes.py $xt_gene_list_file $xt_pc_gene_list_file $gene_annotation_file
fi




########################################
# Standard Univariate fusion-TWAS analysis
########################################
# Create fusion weights in GTEx
if false; then
sed 1d $gtex_tissue_file | while read tissue_name sample_size pseudotissue_name; do
	echo $tissue_name
	sbatch create_fusion_weights_in_a_single_tissue.sh $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_fusion_weights_data_dir $gtex_fusion_weights_dir $ukbb_sumstats_hg38_dir
done
fi


# Create fusion weights position file in each PSEUDOTISSUE (NOTE: PREVIOUS STEP WAS DONE USING TISSUES)
if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat_boolean composit_tissue_string; do
	sbatch generate_meta_analyzed_fusion_weights_and_fusion_pos_file_for_single_pseudotissue.sh $pseudotissue_name $composit_tissue_string $gtex_fusion_weights_data_dir $gtex_fusion_weights_dir $pseudotissue_gtex_fusion_weights_dir
done
fi






########
# Some notes on uncongruent variant ids:
## 1. First allele in bim file is usaually regarded as effect allele
## 2. First allele in (Left to right; called 'Allele1') in BOLT-LMM file is effect allele
## 3. GTEx Ids go 'chrX_POS_a1_a2_b38' where a2 is the effect alele. 
## 3_continued. However, in previous step (preprocess_data_for_genome_wide_ukbb_susie_analysis.sh) we used gtex as reference, but only negated UKBB if left to right in UKBB didnt line up with left to right in gtex variant id (which should be negative of true answer 100% of time??????)
## 4. gtex_fusion_weights_data_dir contains list of gtex variants (and their ids) that are also all in UKBB in 1KG (Adipose_Subcutaneous chosen randomly) as all gtex tissues have the same variants
if false; then
sbatch reformat_variant_ids_in_ukbb_and_1kg_to_be_ammendable_with_gtex_for_fusion.sh $gtex_fusion_weights_data_dir"Adipose_Subcutaneous/" $ukbb_sumstats_hg38_dir $ref_1kg_genotype_dir $gtex_fusion_processed_intermediate_data
fi



########################################
# Generate SuSiE posterior mean causal eqtl effect sizes
########################################
# Preprocess GTEx data for susie
if false; then
sh preprocess_gtex_data_for_susie_analysis_wrapper.sh $gtex_pseudotissue_file $xt_pc_gene_list_file $eqtl_summary_stats_dir $ref_1kg_genotype_dir $ukbb_sumstats_hg38_dir $gtex_preprocessed_for_susie_dir
fi




# Run Susie on GTEx data
gene_file=$gtex_preprocessed_for_susie_dir"susie_input_gene_organization_file.txt"
if false; then
sh run_susie_shell.sh $gene_file $gtex_susie_results_dir
fi


# Create Susie Posterior mean causal effect sizes (PMCES) weights in each tissue
if false; then
sed 1d $gtex_tissue_file | while read tissue_name sample_size pseudotissue_name; do
	sbatch create_susie_pmces_fusion_weights_in_a_single_tissue.sh $tissue_name $pseudotissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_fusion_weights_data_dir $gtex_susie_pmces_fusion_weights_dir $gtex_preprocessed_for_susie_dir $gtex_susie_results_dir
done
fi


# Create Susie Posterior mean causal effect sizes (PMCES) position file in each PSEUDOTISSUE (NOTE: PREVIOUS STEP WAS DONE USING TISSUES)
if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat_boolean composit_tissue_string; do
	echo $pseudotissue_name"\t"$composit_tissue_string
	sbatch generate_susie_pmces_fusion_pos_file_for_single_pseudotissue.sh $pseudotissue_name $composit_tissue_string $gtex_fusion_weights_data_dir $gtex_susie_pmces_fusion_weights_dir $pseudotissue_gtex_susie_pmces_fusion_weights_dir
done
fi

########################################
# Overlap SuSiE components with eQTLs using COLOC
########################################
trait_file=$ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_independent.txt"
gtex_susie_input_data_file=$gtex_preprocessed_for_susie_dir"susie_input_gene_organization_file.txt"
if false; then
sh preprocess_coloc_data.sh $gtex_susie_input_data_file $gtex_fusion_processed_intermediate_data $coloc_data_dir $trait_file $coloc_results_dir $gtex_pseudotissue_file
fi



########################################
# Run TGFM 
#######################################
if false; then

gene_version="cis_heritable_genes"
trait_name="blood_WHITE_COUNT"
samp_size="326723"
sh multivariate_rss_twas_shell.sh $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_fusion_weights_dir $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_set_annotation_file $coloc_results_dir
fi

gene_version="cis_heritable_genes"
trait_name="bp_DIASTOLICadjMEDz"
samp_size="310831"
sh multivariate_rss_twas_shell.sh $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_fusion_weights_dir $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_set_annotation_file $coloc_results_dir
if false; then
gene_version="cis_heritable_genes"
trait_name="biochemistry_Cholesterol"
samp_size="321595"
sh multivariate_rss_twas_shell.sh $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_fusion_weights_dir $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_set_annotation_file $coloc_results_dir

gene_version="cis_heritable_genes"
trait_name="body_WHRadjBMIz"
samp_size="336847"
sh multivariate_rss_twas_shell.sh $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_fusion_weights_dir $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_set_annotation_file $coloc_results_dir

fi
trait_file=$ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_independent.txt"
trait_file=$ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_independent_minus_wbc.txt"
echo $trait_file
echo $trait_file
if false; then
sed 1d $trait_file | while read trait_name study_file samp_size bolt_lmm_h2; do
	echo $trait_name"_"$samp_size
	sh multivariate_rss_twas_shell.sh $trait_name $ukbb_genome_wide_susie_organized_results_dir $gtex_pseudotissue_file $pseudotissue_gtex_fusion_weights_dir $pseudotissue_gtex_susie_pmces_fusion_weights_dir $gtex_fusion_processed_intermediate_data $samp_size $pseudotissue_gtex_rss_multivariate_twas_data_dir $pseudotissue_gtex_rss_multivariate_twas_dir $gene_version $gene_set_annotation_file $coloc_results_dir
done
fi







trait_file=$ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_independent.txt"
if false; then
source ~/.bash_profile
module load R/3.5.1
fi

if false; then
Rscript visualize_joint_rss_twas.R $trait_file $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $rss_joint_multivariate_twas_visualization_dir $gtex_tissue_colors_file
fi


if false; then
Rscript visualize_rss_twas.R $trait_file $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $rss_multivariate_twas_visualization_dir
fi


if false; then
Rscript visualize_rss_twas_with_and_without_genome_prior.R $trait_file $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $rss_multivariate_twas_genome_prior_visualization_dir
fi


if false; then
Rscript visualize_rss_twas_thresholds.R $trait_file $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $rss_multivariate_twas_genome_prior_visualization_dir
fi


if false; then
Rscript visualize_fine_mapped_components.R $trait_file $gtex_pseudotissue_file $pseudotissue_gtex_rss_multivariate_twas_dir $viz_fine_mapped_components_dir $gtex_tissue_colors_file
fi






