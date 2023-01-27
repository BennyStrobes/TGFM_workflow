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
ldsc_baseline_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline_v1.2/"

# LDSC baselineLD annotations (hg19)
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# LDSC 1KG genotype files (hg19)
kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# Gencode hg19 gene annotation file
gencode_gene_annotation_file="/n/groups/price/ben/gene_annotation_files/gencode.v19.annotation.gtf.gz"



############################
# Output directories
############################
# Output roots (one temporary on scratch and one permanent)
temp_output_root="/n/scratch3/users/b/bes710/causal_eqtl_gwas/simulation/"
perm_output_root="/n/groups/price/ben/causal_eqtl_gwas/simulation/"

# Directory containing processed genotype data
processed_genotype_data_dir=$temp_output_root"processed_genotype/"

# Directory containing simulated gene positions
simulated_gene_position_dir=$temp_output_root"simulated_gene_positions/"






############################
# Simulation parameters
############################
# Number of simulated individuals in GWAS
n_gwas_individuals="100000"

# Chromosome to simulate on 
chrom_num="21"

# cis window arround genes to define eQTLs
cis_window="25000"





############################
# Prepare genotype data for analysis:
## 1. Filter number of individuals in original data
## 2. Filter sites to be those in LDSC annotation file
## 3. Convert to plink bed files
############################
if false; then
sh prepare_ukbb_genotype_data_for_simulation_on_single_chromosome.sh $ukbb_genotype_dir $processed_genotype_data_dir $chrom_num $n_gwas_individuals $ldsc_baseline_hg19_annotation_dir $kg_genotype_dir
fi




############################
# Prepare gene file for simulation:
# Genes are defined by actual tss
# Limit to protein coding genes
# In simulation, I will assume gene is expressed in each tissue
############################
simulated_gene_position_file=${simulated_gene_position_dir}"gene_positions_chr"${chrom_num}".txt"
if false; then
sh prepare_simulated_gene_position_list.sh $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file
fi














