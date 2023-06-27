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
# File containing gtex pseudotissues and their assigned tissue category
gtex_pseudotissue_category_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_sample_names/pseudotissue_categories.txt"
# File containing gtex tissue info
gtex_tissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_sample_names/tissue_info.txt"

# File containing gtex genes
xt_gene_list_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/cross_tissue_gene_list.txt"

# GTEx expression dir
gtex_expression_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_expression/"

# GTEx covariate dir
gtex_covariate_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_covariates/"

# GTEx genotype dir
gtex_genotype_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_genotype/"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"

# Drug target gene list
# Downloaded here: https://www.nature.com/articles/s41588-022-01167-z on 6/9/23
drug_target_gene_list_file="/n/groups/price/ben/gene_annotation_files/scdrs_gold_geneset.txt"


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

# Directory containing epimap data
epimap_data_dir="/n/groups/price/nolan/epimap/"

# Mod-ldscore regression code
curr_dir=`pwd`
mod_ldsc_code_dir=$curr_dir"/modified_sldsc/"

# Summary statistics
full_sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"

# hg38 sldsc weights
sldsc_h38_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

# Directory containing quasi-independent ld bllocks
quasi_independent_dir="/n/groups/price/ben/quasi_independent_ld_blocks/"

# Directory containing pops results
# Downloaded from https://www.finucanelab.org/data on 6/20/23
pops_results_summary_file="/n/groups/price/ben/pops_data/PoPS_FULLRESULTS.txt.gz"

# Input single cell expression data 
input_sc_h5py_file="/n/groups/price/scdata/Perez_Science_2021/lupus.h5ad"

# pseudobulk h5ad file
input_sc_h5py_pseudobulk_file="/n/groups/price/scdata/Perez_Science_2021/lupus-pseudobulk.h5ad"

# Gene annotation file (hg19)
hg19_gene_annotation_file="/n/groups/price/ben/reference_data/gene_annotation_files/gencode.v19.annotation.gff3"

# Gene annotation file (hg38)
hg38_gene_annotation_file="/n/groups/price/ben/reference_data/gene_annotation_files/gencode.v38.annotation.gtf.gz"

# Directory containing genotype data
sc_genotype_data_dir="/n/groups/price/scdata/Perez_Science_2021/genotype/plink_maf10/"

# individual file
sc_individual_info_file="/n/groups/price/scdata/Perez_Science_2021/genotype/indv_info_from_scdata.tsv"





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
ukbb_preprocessed_for_genome_wide_susie_dir=$perm_output_root"ukbb_preprocessed_for_genome_wide_susie/"

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

# Directory containing TGFM iterative prior results
iterative_tgfm_prior_results_dir=$perm_output_root"iterative_tgfm_prior/"

# Directory containing visualizations of TGFM h2 estimates
visualize_tgfm_h2_dir=$output_root"visualize_tgfm_heritability_estimates/"

# Directory containing standard sldsc processed data
standard_sldsc_processed_data_dir=$output_root"standard_sldsc_processed_data/"

# Directory containing standard sldsc results
standard_sldsc_results_dir=$output_root"standard_sldsc_results/"

# Directory containing drug target gene set enrichment analyses
drug_target_gene_set_enrichment_dir=$perm_output_root"drug_target_gene_set_enrichment/"

# Directory containing Epimap enrichment raw data
epimap_enrichment_raw_data_dir=$output_root"raw_epimap_enrichment_data/"

# Directory containing Epimap enrichments
epimap_enrichment_dir=$perm_output_root"epimap_enrichment/"

#Directory containing pops enrichments
pops_enrichment_dir=$perm_output_root"pops_enrichment/"

# Sparse heritability visualization dir
visualize_sparse_h2_dir=$output_root"visualize_sparse_h2/"

tgfm_sldsc_results_dir=$perm_output_root"tgfm_sldsc_results/"

visualize_tgfm_sldsc_dir=$output_root"visualize_tgfm_sldsc/"

visualize_tgfm_dir=$output_root"visualize_tgfm/"

# Directory containing processed sc expression data
processed_sc_expression_dir=$output_root"processed_sc_expression/"

# Processed single cell genotype dir
processed_sc_genotype_dir=$output_root"processed_sc_genotype/"

# Processed single cell pseudobulk expression
sc_pseudobulk_expression_dir=$output_root"sc_pseudobulk_expression/"

# Visualize processed single cell expression dir
visualize_processed_sc_expression_dir=$perm_output_root"visualize_processed_sc_expression/"

# Input data for fusion
sc_fusion_input_dir=$output_root"sc_fusion_input/"

# Directory containing sc PBMC Susie gene models
sc_pbmc_susie_gene_models_dir=$perm_output_root"sc_pbmc_susie_gene_models/"


# Directory containing preprocessed sc TGFM data
preprocessed_sc_tgfm_data_dir=$output_root"preprocessed_sc_tgfm_data/"

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
# NOTE. THIS CURRENTLY HAS A WEIRD HACK IN THE CHROMOSOME ANALYSIS TO DEAL WITH LD CORRECTION AND MISSING DATA
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


# Only components of genes
gene_type="cis_heritable_gene"
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
sh organize_gene_ld_scores_for_tgfm_sldsc.sh $gtex_pseudotissue_file $preprocessed_tgfm_sldsc_data_dir $gene_type $gtex_susie_gene_models_dir
fi


########################################
# Run TGFM-S-LDSC
########################################
if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
	sbatch run_tgfm_sldsc.sh $preprocessed_tgfm_sldsc_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $tgfm_sldsc_results_dir $trait_name $mod_ldsc_code_dir $quasi_independent_ld_blocks_hg38_dir
done
fi





########################################
# Visualize TGLR results
########################################
if false; then
source ~/.bash_profile
module load R/3.5.1
Rscript visualize_tgfm_sldsc_results.R $independent_trait_names_file $tgfm_sldsc_results_dir $preprocessed_tgfm_sldsc_data_dir $gtex_tissue_colors_file $visualize_tgfm_sldsc_dir
fi

########################################
# Preprocess data for TGFM
########################################
# Number of parallel jobs
num_jobs="40"
#gene_type="cis_heritable_gene"
gene_type="component_gene"
# FIle summarizing ukkbb windows
ukkbb_window_summary_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows_and_processed_data.txt"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch preprocess_data_for_tgfm.sh $ukkbb_window_summary_file $hapmap3_snpid_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $job_number $num_jobs $gene_type $preprocessed_tgfm_sldsc_data_dir $tgfm_sldsc_results_dir
done
fi

# Organize preprocessed TGFM results across parallel jobs
if false; then
sh organize_processed_tgfm_input_data.sh $num_jobs $gene_type $preprocessed_tgfm_data_dir
fi



########################################
# Run TGFM
########################################
gene_type="component_gene"
num_jobs="8"

if false; then
trait_name="blood_MONOCYTE_COUNT"
for job_number in $(seq 0 $(($num_jobs-1))); do
	tgfm_input_summary_file=${preprocessed_tgfm_data_dir}${gene_type}"_tgfm_input_data_summary.txt"
	tgfm_output_stem=${tgfm_results_dir}"tgfm_results_"${trait_name}"_"${gene_type}
	sbatch run_tgfm_shell.sh $trait_name $tgfm_input_summary_file $tgfm_output_stem $gtex_pseudotissue_file $job_number $num_jobs
done
fi

if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
	for job_number in $(seq 0 $(($num_jobs-1))); do
		tgfm_input_summary_file=${preprocessed_tgfm_data_dir}${gene_type}"_tgfm_input_data_summary.txt"
		tgfm_output_stem=${tgfm_results_dir}"tgfm_results_"${trait_name}"_"${gene_type}
		sbatch run_tgfm_shell.sh $trait_name $tgfm_input_summary_file $tgfm_output_stem $gtex_pseudotissue_file $job_number $num_jobs
	done
done
fi







########################################
# Compute iterative prior
########################################
tgfm_input_summary_file=${preprocessed_tgfm_data_dir}${gene_type}"_tgfm_input_data_summary.txt"

if false; then
trait_name="blood_RED_COUNT"
	tgfm_output_stem=${tgfm_results_dir}"tgfm_results_"${trait_name}"_"${gene_type}
	sbatch learn_iterative_tgfm_component_prior.sh $trait_name $tgfm_output_stem $gtex_pseudotissue_file ${preprocessed_tgfm_data_dir}${gene_type} $tgfm_input_summary_file $iterative_tgfm_prior_results_dir
fi
if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
	tgfm_output_stem=${tgfm_results_dir}"tgfm_results_"${trait_name}"_"${gene_type}
	sbatch learn_iterative_tgfm_component_prior.sh $trait_name $tgfm_output_stem $gtex_pseudotissue_file ${preprocessed_tgfm_data_dir}${gene_type} $tgfm_input_summary_file $iterative_tgfm_prior_results_dir
done
fi


#################################
# Run TGFM with iterative prior
#################################
gene_type="component_gene"
num_jobs="8"
if false; then
sed 1d $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2.txt" | while read trait_name study_file sample_size h2; do
echo $trait_name

for job_number in $(seq 0 $(($num_jobs-1))); do
	tgfm_input_summary_file=${preprocessed_tgfm_data_dir}${gene_type}"_tgfm_input_data_summary.txt"
	tgfm_output_stem=${tgfm_results_dir}"tgfm_results_"${trait_name}"_"${gene_type}
	sbatch run_tgfm_with_iterative_prior_shell.sh $trait_name $tgfm_input_summary_file $tgfm_output_stem $gtex_pseudotissue_file $job_number $num_jobs
done
done
fi





#################################
# Organize TGFM Results across parallel runs
#################################
if false; then
sbatch organize_tgfm_results_across_parallel_runs.sh $tgfm_results_dir $gene_type $num_jobs $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_readable.txt" $gtex_pseudotissue_file $gtex_pseudotissue_category_file ${preprocessed_tgfm_data_dir}${gene_type} $ukbb_preprocessed_for_genome_wide_susie_dir $tgfm_sldsc_results_dir
fi

#################################
# Run drug target gene set enrichment analysis
#################################
if false; then
sh run_drug_target_gene_set_enrichment_analysis.sh $tgfm_results_dir $gene_type $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_readable.txt" $gene_annotation_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $drug_target_gene_list_file $drug_target_gene_set_enrichment_dir 
fi

#################################
# Run POPs enrichment analysis
#################################
if false; then
sh run_pops_enrichment_analysis.sh $tgfm_results_dir $gene_type $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_readable.txt" $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $pops_results_summary_file $pops_enrichment_dir 
fi

#################################
# Run epimap cell type enrichment analysis
#################################
independent_trait_list_file=$ukbb_sumstats_hg38_dir"ukbb_hg38_independent_sumstat_files_with_samp_size_and_h2_readable.txt"
if false; then
sh run_epimap_enrichment_analysis.sh $tgfm_results_dir $gtex_susie_gene_models_dir $liftover_directory $gtex_pseudotissue_file $independent_trait_list_file $epimap_data_dir $epimap_enrichment_raw_data_dir $epimap_enrichment_dir
fi


#################################
# Run gene set enrichment analysis
#################################
if false; then
sh run_gene_set_enrichment_analysis.sh $tgfm_results_dir $gene_type $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_readable.txt" $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $gene_set_enrichment_dir 
fi
#################################
# Visualize TGFM results
#################################
if false; then
source ~/.bash_profile
module load R/3.5.1
fi
if false; then
Rscript visualize_tgfm_results.R $ukbb_sumstats_hg38_dir"ukbb_hg38_sumstat_files_with_samp_size_and_h2_readable.txt" $tgfm_sldsc_results_dir $tgfm_results_dir $preprocessed_tgfm_sldsc_data_dir $gtex_tissue_colors_file $visualize_tgfm_dir $iterative_tgfm_prior_results_dir $epimap_enrichment_dir
fi










#################################
# Process SC expression
#################################
if false; then
sbatch process_sc_expression.sh $input_sc_h5py_file $input_sc_h5py_pseudobulk_file $processed_sc_expression_dir $sc_individual_info_file $sc_genotype_data_dir $visualize_processed_sc_expression_dir $processed_sc_genotype_dir 
fi

#################################
# Process SC genotype data
#################################
filtered_sample_info_file=$processed_sc_expression_dir"individual_info_european_rna_and_dna.txt"
if false; then
sh process_sc_genotype.sh $processed_sc_genotype_dir $sc_genotype_data_dir $filtered_sample_info_file $liftover_directory $ukbb_preprocessed_for_genome_wide_susie_dir
fi




#################################
# Generate pseudobulk expression
#################################
if false; then
sh generate_pseudobulk_expression.sh $processed_sc_expression_dir $input_sc_h5py_file $processed_sc_genotype_dir $sc_pseudobulk_expression_dir $hg19_gene_annotation_file $gene_annotation_file $xt_gene_list_file
fi


#################################
# Prepare sc expression data for fusion weights analysiss
#################################
pb_cell_type_file=${sc_pseudobulk_expression_dir}"pseudobulk_data_set_summary_filtered.txt"
num_parallel_jobs="10"
echo $pb_cell_type_file
if false; then
sed 1d $pb_cell_type_file | while read data_set_name cluster_method pb_cell_type_name expr_file cov_file num_donors num_genes num_cells_per_indi_file; do
	sh preprocess_sc_data_for_fusion_weights_analysis_in_single_pseudobulk_cell_type.sh $pb_cell_type_name $sc_pseudobulk_expression_dir $sc_fusion_input_dir $num_parallel_jobs $processed_sc_genotype_dir
done
fi


########################################
# Run sc pseudobulk gene model analysis (in each pseudobulk cell type seperately)
########################################
pb_cell_type_name="B"
job_number="0"
if false; then
sed 1d $pb_cell_type_file | while read data_set_name cluster_method pb_cell_type_name expr_file cov_file num_donors num_genes num_cells_per_indi_file; do
	for job_number in $(seq 0 $(($num_parallel_jobs-1))); do
		sbatch create_susie_gene_model_in_a_single_pseudobulk_cell_type.sh $pb_cell_type_name $pb_cell_type_name $sc_fusion_input_dir $sc_pbmc_susie_gene_models_dir $job_number
	done
done
fi


########################################
# Organize sc pseudobulk gene model results (create pos file)
########################################
if false; then
sed 1d $pb_cell_type_file | while read data_set_name cluster_method pb_cell_type_name expr_file cov_file num_donors num_genes num_cells_per_indi_file; do
	sbatch organize_susie_gene_model_results_in_a_single_pseudotissue.sh $pb_cell_type_name $sc_fusion_input_dir${pb_cell_type_name}"/" $sc_pbmc_susie_gene_models_dir
done
fi



########################################
# Preprocess data for SC TGFM
########################################
# Number of parallel jobs
num_jobs="40"
#gene_type="cis_heritable_gene"
gene_type="component_gene"
# FIle summarizing ukkbb windows
ukkbb_window_summary_file=$ukbb_preprocessed_for_genome_wide_susie_dir"genome_wide_susie_windows_and_processed_data.txt"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch preprocess_data_for_sc_tgfm.sh $ukkbb_window_summary_file $hapmap3_snpid_file $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_sc_tgfm_data_dir $job_number $num_jobs $gene_type $preprocessed_tgfm_sldsc_data_dir $pb_cell_type_file $sc_pbmc_susie_gene_models_dir
done
fi





