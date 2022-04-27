
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
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"

# Genotype data from 1KG
ref_1kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# GTEx expression dir
gtex_expression_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/"

# GTEx covariate dir
gtex_covariate_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_covariates/"

# GTEx genotype dir
gtex_genotype_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_genotype/"


# EpiMap enhancer-gene links
epimap_input_dir="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas/input_data/epimap/per_group_enhancer_gene_links/"

# ABC enhancer gene links
abc_input_dir="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas/input_data/abc/"

# ABC-enhancer, gwas-variant enrichment across biosamples (supp table 6 of Nassar et al)
# wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03446-x/MediaObjects/41586_2021_3446_MOESM8_ESM.txt
abc_enhancer_gwas_variant_enrichment_file="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas/input_data/abc/gwas_variant_abc_enhancer_enrichments_across_biosamples_supp_table_6.txt"


##################
# Output data
##################
# Output root directory
output_root="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas/"

# Directory containing hg38 ukbb summary stats
ukbb_sumstats_hg38_dir=$output_root"ukbb_sumstats_hg38/"

# Directory containing UKBB sumstats for genome-wide susie
ukbb_preprocessed_for_genome_wide_susie_dir=$output_root"ukbb_preprocessed_for_genome_wide_susie/"

# Directory containing Susie results when applied to genome-wide UKBB
ukbb_genome_wide_susie_results_dir=$output_root"ukbb_genome_wide_susie_results/"

# Directory containing Susie organized results when applied to genome-wide UKBB
ukbb_genome_wide_susie_organized_results_dir=$output_root"ukbb_genome_wide_susie_organized_results/"

# Directory containing epimap ukbb component overlap
epimap_ukbb_genome_wide_susie_overlap_dir=$output_root"epimap_ukbb_genome_wide_susie_overlap/"

# Directory containing abc ukbb component overlap
abc_ukbb_genome_wide_susie_overlap_dir=$output_root"abc_ukbb_genome_wide_susie_overlap/"

# Directory containing visualizations of tissue component mapping
tissue_component_mapping_viz=$output_root"visualize_tissue_mapping_ukbb_genome_wide_susie/"


# Directory containing hg38 GTEx preprocessed for Susie
gtex_gene_set_dir_dir=$output_root"gtex_gene_sets/"

# Directory containing hg38 GTEx preprocessed for Susie
gtex_preprocessed_for_susie_dir=$output_root"gtex_preprocessed_for_susie/"

# Directory containing Susie results when applied to GTEx
gtex_susie_results_dir=$output_root"gtex_susie_results/"

# Directory containing UKBB summary stats preprocessed for Susie
ukbb_preprocessed_for_susie_dir=$output_root"ukbb_preprocessed_for_susie/"

# Directory containing Susie results when applied to UKBB
ukbb_susie_results_dir=$output_root"ukbb_susie_results/"

# Directory containing causal effect size regression of GTEx onto traits
gtex_onto_trait_causal_effect_size_regression_dir=$output_root"gtex_onto_trait_causal_effect_size_regression/"

# Directory containing visualization of causal effect size regression of GTEx onto traits
gtex_onto_trait_causal_effect_size_regression_viz_dir=$output_root"viz_gtex_onto_trait_causal_effect_size_regression/"

# Directory containing GTEx fusion weights input
gtex_fusion_weights_data_dir=$output_root"gtex_fusion_weights_data/"

# Directory containing GTEx fusion weights
gtex_fusion_weights_dir=$output_root"gtex_fusion_weights/"

# Directory containing GTEx fusion associations
gtex_fusion_associations_dir=$output_root"gtex_fusion_associations/"





##################
# Analysis
##################


# Liftover UKBB summary statistics to hg38
if false; then
sbatch liftover_ukbb_summary_statistics_from_hg19_to_hg38.sh $liftover_directory $ukbb_sumstats_hg19_dir $ukbb_sumstats_hg38_dir
fi



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
sh organize_susie_genome_wide_results.sh $window_file $ukbb_genome_wide_susie_results_dir $ukbb_genome_wide_susie_organized_results_dir
fi




# Epimap ukbb component overlap
global_ukbb_component_file=$ukbb_genome_wide_susie_organized_results_dir"organized_susie_components_files.txt"
if false; then
sh epimap_ukbb_component_overlap.sh $global_ukbb_component_file $epimap_input_dir $epimap_ukbb_genome_wide_susie_overlap_dir $liftover_directory
fi

# ABC UKBB component overlap
global_ukbb_component_file=$ukbb_genome_wide_susie_organized_results_dir"organized_susie_components_files.txt"
if false; then
sbatch abc_ukbb_component_overlap.sh $global_ukbb_component_file $abc_input_dir $abc_ukbb_genome_wide_susie_overlap_dir $liftover_directory
fi










if false; then
module load R/3.5.1
fi
if false; then
Rscript visualize_tissue_mapping_of_ukbb_components.R $tissue_component_mapping_viz $ukbb_genome_wide_susie_organized_results_dir $epimap_ukbb_genome_wide_susie_overlap_dir $abc_ukbb_genome_wide_susie_overlap_dir $abc_enhancer_gwas_variant_enrichment_file
fi







# Filter gene list to protein coding genes
xt_pc_gene_list_file=$gtex_gene_set_dir_dir"cross_tissue_protein_coding_gene_list.txt"
if false; then
python3 filter_gene_list_to_protein_coding_genes.py $xt_gene_list_file $xt_pc_gene_list_file $gene_annotation_file
fi

# Preprocess GTEx data for susie
if false; then
sh preprocess_gtex_data_for_susie_analysis_wrapper.sh $gtex_pseudotissue_file $xt_pc_gene_list_file $eqtl_summary_stats_dir $ref_1kg_genotype_dir $ukbb_sumstats_hg38_dir $gtex_preprocessed_for_susie_dir
fi



# Run Susie on GTEx data
gene_file=$gtex_preprocessed_for_susie_dir"susie_input_gene_organization_file.txt"
if false; then
sh run_susie_shell.sh $gene_file $gtex_susie_results_dir
fi


# Preprocess UKBB data for susie
gtex_gene_file=$gtex_preprocessed_for_susie_dir"susie_input_gene_organization_file.txt"
if false; then
sbatch preprocess_ukbb_data_for_susie_analysis_wrapper.sh $gtex_gene_file $ukbb_sumstats_hg38_dir $ukbb_preprocessed_for_susie_dir
fi



# UKBB Susie dir
ukbb_gene_file=$ukbb_preprocessed_for_susie_dir"susie_input_ukkbb_gene_organization_file.txt"
if false; then
sh run_susie_shell.sh $ukbb_gene_file $ukbb_susie_results_dir
fi


# Regress posterior mean causal gtex effect sizes onto trait posterior mean causal effect sizes
# causal effect size regression
trait_name="blood_WHITE_COUNT"
if false; then
sh run_causal_effect_size_regression.sh $trait_name $gtex_gene_file $gtex_susie_results_dir $ukbb_susie_results_dir $gtex_pseudotissue_file $gtex_onto_trait_causal_effect_size_regression_dir $gtex_onto_trait_causal_effect_size_regression_viz_dir
fi




# Create fusion weights in GTEx
tissue_name="Adipose_Subcutaneous"
if false; then
sh create_fusion_weights_in_a_single_tissue.sh $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_fusion_weights_data_dir $gtex_fusion_weights_dir
fi


if false; then
sed 1d $gtex_tissue_file | while read tissue_name sample_size pseudotissue_name; do
	echo $tissue_name
	sh create_fusion_weights_in_a_single_tissue.sh $tissue_name $xt_pc_gene_list_file $gtex_expression_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_fusion_weights_data_dir $gtex_fusion_weights_dir
done
fi




tissue_name="Adipose_Subcutaneous"
if false; then
sh run_fusion_trait_association.sh $tissue_name $gtex_fusion_weights_dir $gtex_fusion_associations_dir
fi







