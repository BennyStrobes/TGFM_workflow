#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-9:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=14GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
per_element_heritability="$8"
n_non_mediated_variants_per_gene="$9"
fraction_causal_genes="${10}"
simulated_gene_expression_dir="${11}"
simulated_learned_gene_models_dir="${12}"
simulated_trait_dir="${13}"
simulated_gwas_dir="${14}"
simulated_twas_dir="${15}"


source ~/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date


#######################################################
# Step 1: Simulate gene expression and fit gene models
#######################################################
echo "Simulation Step 1"
python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir


#######################################################
# Step 2: Simulate trait values
#######################################################
echo "Simulation Step 2"
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $per_element_heritability $n_non_mediated_variants_per_gene $fraction_causal_genes $simulated_trait_dir $n_gwas_individuals $processed_genotype_data_dir



#######################################################
# Step 3: Run GWAS on all snps
#######################################################
echo "Simulation Step 3"
python3 run_gwas_on_simulated_trait_at_all_snps.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulated_gwas_dir $simulated_gene_expression_dir


#######################################################
# Step 4: Run TWAS
#######################################################
echo "Simulation Step 4"
python3 run_various_twas_methods.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_trait_dir $simulated_twas_dir









##############
# OLD
##############

if false; then
#######################################################
# Step 1: Simulate gene expression and fit gene models
#######################################################
echo "Simulation Step 1"
python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir

#######################################################
# Step 2: Simulate trait values
#######################################################
echo "Simulation Step 2"
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_trait_dir $n_gwas_individuals

#######################################################
# Step 3: Run GWAS on only hapmap3 snps
#######################################################
echo "Simulation Step 3"
python3 run_gwas_on_simulated_trait_at_only_hapmap3_snps.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $ldsc_weights_dir $simulated_gwas_dir

#######################################################
# Step 4: Generate gene ld-scores
#######################################################
echo "Simulation Step 4"
python3 generate_gene_ld_scores.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_ld_scores_dir

#######################################################
# Step 5: Generate variant annotation-weighted ld-scores
#######################################################
echo "Simulation Step 5"
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
python ${ldsc_code_dir}ldsc.py\
	--l2\
	--bfile ${processed_genotype_data_dir}100G.EUR.QC.filtered.${chrom_num}\
	--ld-wind-cm 1\
	--annot ${processed_genotype_data_dir}baseline.${chrom_num}.annot\
	--out ${simulated_ld_scores_dir}${simulation_name_string}"_baseline".${chrom_num}\
	--print-snps ${simulated_ld_scores_dir}${simulation_name_string}"_regression_snp_ids.txt"


#######################################################
# Step 5b: Generate LD score weights
#######################################################
echo "Simulation Step 5b"
# Filter genotype data to just regression snps in 1KG
source ~/.bash_profile
plink2 --bfile ${processed_genotype_data_dir}"100G.EUR.QC.filtered."${chrom_num} --extract ${simulated_ld_scores_dir}${simulation_name_string}"_regression_snp_ids.txt" --threads 1 --make-bed --keep-allele-order --out ${simulated_ld_scores_dir}${simulation_name_string}"_100G_regression_snps_only."${chrom_num}
# Create regression snp weights
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12
python ${ldsc_code_dir}ldsc.py\
	--l2\
	--bfile ${simulated_ld_scores_dir}${simulation_name_string}"_100G_regression_snps_only."${chrom_num}\
	--ld-wind-cm 1\
	--out ${simulated_ld_scores_dir}${simulation_name_string}"_regression_weights".${chrom_num}\
	--print-snps ${simulated_ld_scores_dir}${simulation_name_string}"_regression_snp_ids.txt"
# Delete uncessary plink file
rm ${simulated_ld_scores_dir}${simulation_name_string}"_100G_regression_snps_only."${chrom_num}*

 
#######################################################
# Step 6: Organize data for tgfm-sldsc
#######################################################
echo "Simulation Step 6"
source ~/.bash_profile
python3 organize_data_for_sldsc.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $processed_genotype_data_dir $simulated_gwas_dir $ldsc_weights_dir $simulated_ld_scores_dir 

#######################################################
# Step 7: Run TGFM-sldsc
#######################################################
echo "Simulation Step 7"
# Use eQTL PMCES
eqtl_sample_size_arr=( "100" "200" "300" "500" "1000" "inf")

for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	trait_file=${simulated_ld_scores_dir}${simulation_name_string}"_ldsc_ready_summary_statistics.txt"
	python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --n-blocks 200 --chisq-max 1000 --ref-ld ${simulated_ld_scores_dir}${simulation_name_string}"_joint_baseline_variant_"${eqtl_sample_size}"_gene_ld_scores" --w-ld ${simulated_ld_scores_dir}${simulation_name_string}"_regression_weights."${chrom_num} --print-delete-vals --print-coefficients --out ${simulated_sldsc_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_sldsc_results"

	# Organize sldsc results
	source ~/.bash_profile
	python3 organize_tgfm_sldsc_results.py ${simulated_sldsc_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_pmces_sldsc_results" ${simulated_ld_scores_dir}${simulation_name_string}"_joint_baseline_variant_"${eqtl_sample_size}"_gene_ld_scores" $n_gwas_individuals
done

# Use susie distr eQTLS
eqtl_sample_size_arr=( "100" "200" "300" "500" "1000")

for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	source /n/groups/price/ben/environments/sldsc/bin/activate
	module load python/2.7.12
	trait_file=${simulated_ld_scores_dir}${simulation_name_string}"_ldsc_ready_summary_statistics.txt"
	python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --n-blocks 200 --chisq-max 1000 --ref-ld ${simulated_ld_scores_dir}${simulation_name_string}"_joint_baseline_variant_"${eqtl_sample_size}"_susie_distr_gene_ld_scores" --w-ld ${simulated_ld_scores_dir}${simulation_name_string}"_regression_weights."${chrom_num} --print-delete-vals --print-coefficients --out ${simulated_sldsc_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_distr_sldsc_results"

	# Organize sldsc results
	source ~/.bash_profile
	python3 organize_tgfm_sldsc_results.py ${simulated_sldsc_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_susie_distr_sldsc_results" ${simulated_ld_scores_dir}${simulation_name_string}"_joint_baseline_variant_"${eqtl_sample_size}"_susie_distr_gene_ld_scores" $n_gwas_individuals
done


#######################################################
# Step 8: Randomly select chromosome windows to run tgfm on
#######################################################
echo "Simulation Step 8"
global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
python3 randomly_select_windows_to_run_tgfm_on_in_simulation.py $global_window_file $simulation_window_list_file $simulation_number




#######################################################
# Step 9: Run GWAS on simulated trait on only snps in TGFM windows.
# Also computes in-sample LD for TGFM windows
#######################################################
echo "Simulation Step 9"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
python3 run_gwas_on_simulated_trait_at_snps_in_tgfm_windows.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulation_window_list_file $simulated_gwas_dir $simulated_tgfm_input_data_dir



#######################################################
# Step 10: Preprocess data for TGFM
#######################################################
echo "Simulation Step 10"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "100" "200" "300" "500" "1000" "inf" )
eqtl_type="susie_pmces"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_sldsc_results_dir $simulated_tgfm_input_data_dir $eqtl_type
done


simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "100" "200" "300" "500" "1000")
eqtl_type="susie_distr"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_sldsc_results_dir $simulated_tgfm_input_data_dir $eqtl_type
done





fi







date















