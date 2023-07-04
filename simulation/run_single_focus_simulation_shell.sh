#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=60GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
simulated_gene_expression_dir="$8"
simulated_learned_gene_models_dir="$9"
simulated_trait_dir="${10}"
simulated_gwas_dir="${11}"
simulated_tgfm_input_data_dir="${12}"
eqtl_sample_size="${13}"
simulated_focus_input_dir="${14}"
simulated_focus_results_dir="${15}"



# Create directory specific to this precise simulation
focus_gene_models_dir=${simulated_focus_input_dir}"simulation_"${simulation_number}"/"
mkdir $focus_gene_models_dir

echo $simulation_name_string
echo $eqtl_sample_size


##############################
# Preprocess gene model data to put into fusion format
source ~/.bash_profile
python3 prepare_gene_models_for_focus.py $simulation_number $chrom_num $simulation_name_string $simulated_gene_expression_dir $simulated_learned_gene_models_dir $eqtl_sample_size $processed_genotype_data_dir $focus_gene_models_dir

module load R/4.0.1
focus_gene_model_python_file=${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_gene_models_python_summary.txt"
focus_gene_model_pos_file=${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_gene_models_python_summary.pos"
Rscript prepare_gene_models_for_focus_pt2.R $simulation_number $chrom_num $simulation_name_string $eqtl_sample_size $focus_gene_models_dir $focus_gene_model_python_file $focus_gene_model_pos_file




##############################
# Now convert gene models to FOCUS gene model input
module load gcc/9.2.0
module load python/3.9.14
module load R/4.1.2
focus import ${focus_gene_model_pos_file} fusion --tissue all --name ${simulation_name_string}"_eqtlss_"${eqtl_sample_size} --assay rnaseq --output ${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_focus_input"




##############################
# Generate gwas data
source ~/.bash_profile
focus_gwas_summary_stat_file=${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_gwas_summary_stats.txt"
global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"
python3 generate_gwas_input_for_focus.py $global_window_file $simulation_number $chrom_num $simulation_name_string $eqtl_sample_size ${simulated_gwas_dir} $processed_genotype_data_dir $n_gwas_individuals $focus_gwas_summary_stat_file

module load gcc/9.2.0
module load python/3.9.14
module load R/4.1.2
focus_cleaned_gwas_summary_stat_stem=${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_cleaned_gwas_summary_stats"
focus munge $focus_gwas_summary_stat_file --output $focus_cleaned_gwas_summary_stat_stem



##############################
# Delete unnecessary files
rm ${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}*"RDat"
rm ${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}*"snps_file.txt"
rm ${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}*"susie_weights.txt"
rm ${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}*"gwas_summary_stats.txt"


##############################
# Run focus fine-mapping
module load gcc/9.2.0
module load python/3.9.14
module load R/4.1.2
sum_stats_file=${focus_cleaned_gwas_summary_stat_stem}".sumstats.gz"
genotype_file=${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}
db_file=${focus_gene_models_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_focus_input.db"
output_stem=${simulated_focus_results_dir}${simulation_name_string}"_eqtlss_"${eqtl_sample_size}"_focus_res"
focus finemap $sum_stats_file $genotype_file $db_file --chr ${chrom_num} --out $output_stem

