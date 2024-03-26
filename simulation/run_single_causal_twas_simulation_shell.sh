#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=50GB                         # Memory total in MiB (for all cores)




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
simulated_causal_twas_input_data_dir="${12}"
simulated_causal_twas_gene_models_dir="${13}"
simulated_causal_twas_results_dir="${14}"
eqtl_sample_size="${15}"
processed_ctwas_genotype_data_dir="${16}"

source ~/.bash_profile

module load gcc/9.2.0
module load R/4.3.1


# Prepare gwas summary statistics for causal-TWAS
##*#*#
input_gwas_file=$simulated_gwas_dir"simulation_realistic_qtl_ss_"$simulation_number"_chrom1_cis_window_100000_ss_100000_ge_h2_075_gt_arch_1_caus_t_qtl_arch_default_merged_gwas_summary_stats.txt"
###
ctwas_gwas_sumstats_file=$simulated_causal_twas_input_data_dir$simulation_name_string"_eqtl_ss_"${eqtl_sample_size}"_ctwas_gwas_sumstats.RDS"
ctwas_gwas_ss_file=$simulated_causal_twas_input_data_dir$simulation_name_string"_eqtl_ss_"${eqtl_sample_size}"_ctwas_gwas_ss.RDS"
Rscript prepare_gwas_summary_statistics_for_causal_twas.R $input_gwas_file $n_gwas_individuals $ctwas_gwas_sumstats_file $ctwas_gwas_ss_file




# Prepare gene models for causal-TWAS
# Make gene model sub-directory (for just this simulation)
new_gene_model_dir=${simulated_causal_twas_gene_models_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}
mkdir $new_gene_model_dir
mkdir $new_gene_model_dir"/multitissue"
##*#*#
input_causal_eqtl_effects_file=$simulated_gene_expression_dir"simulation_realistic_qtl_ss_"$simulation_number"_chrom1_cis_window_100000_ss_100000_ge_h2_075_gt_arch_1_caus_t_qtl_arch_default""_causal_eqtl_effect_summary.txt"
simulated_learned_gene_models_stem=$simulated_learned_gene_models_dir"simulation_"$simulation_number"/simulation_realistic_qtl_ss_"$simulation_number"_chrom1_cis_window_100000_ss_100000_ge_h2_075_gt_arch_1_caus_t_qtl_arch_default_"
######
genotype_bim_file=$processed_genotype_data_dir"simulated_gwas_data_"${chrom_num}".bim"
######
pos_file=$new_gene_model_dir"/multitissue.pos"
tmp_pos_file=$new_gene_model_dir"/multitissue.tmp_pos"
sh prepare_gene_models_for_causal_twas.sh $input_causal_eqtl_effects_file $simulated_learned_gene_models_stem $genotype_bim_file $pos_file $tmp_pos_file $new_gene_model_dir"/multitissue/multitissue" $new_gene_model_dir"/"


Rscript run_ctwas.R $ctwas_gwas_sumstats_file $processed_ctwas_genotype_data_dir $new_gene_model_dir"/multitissue" $simulated_causal_twas_results_dir${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}
