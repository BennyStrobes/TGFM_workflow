
##################
# Input data
##################

# Directory containing summary statistics
summary_stat_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_MAF0.001_v3/"

#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####2. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory="/n/groups/price/ben/tools/liftOver_x86/"

# File containing gtex tissues to do analysis on and their sample size
gtex_tissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_sample_names/pseudotissue_info.txt"

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

# Direcotry containing 1KG European genotype data
ref_1kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"









##################
# Output data
##################
# Output root directory
output_root="/n/groups/price/ben/causal_eqtl_gwas/gtex_v8_causal_eqtl_gwas/"

# Directory containing all gtex associations (in gene-level organized format)
processed_gtex_associations_gene_level_dir=$output_root"processed_gtex_associations_gene_level/"

# Directory containing reference genotype data from 1kg EUR. 
# Directory contains seperate file for each gene
reference_genotype_gene_level_dir=$output_root"reference_1kg_eur_genotype_gene_level/"

# Directory processed trait associations (in gene-level organized format)
processed_trait_associations_gene_level_dir=$output_root"processed_trait_associations_gene_level/"
