#############################
# Simulation of:
####### 1. TGFM-SLDSC
####### 2. TGFM fine-mapping
##############################




############################
# Input data
############################

# Directory created by Martin containing UKBB genotype for 334K unrelated European individuals
ukbb_genotype_dir="/n/scratch3/users/j/jz286/imp_geno/"

# Ldsc weights for hg19
ldsc_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/"

# LDSC baseline annotations (hg19)
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"

# LDSC baselineLD annotations (hg19)
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
ldsc_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"



############################
# Output directories
############################
# Output roots (one temporary on scratch and one permanent)
temp_output_root="/n/scratch3/users/b/bes710/causal_eqtl_gwas/simulation/"
perm_output_root="/n/groups/price/ben/causal_eqtl_gwas/simulation/"

# Directory containing processed genotype data
processed_genotype_data_dir=$temp_output_root"processed_genotype/"






############################
# Simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"
