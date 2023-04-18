#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-25:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=22GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string_orig="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
ldsc_real_data_results_dir="$8"
per_element_heritability="$9"
total_heritability="${10}"
fraction_expression_mediated_heritability="${11}"
num_windows_per_sim="${12}"
simulated_gene_expression_dir="${13}"
simulated_learned_gene_models_base_dir="${14}"
simulated_trait_dir="${15}"
simulated_gwas_dir="${16}"
simulated_tgfm_input_data_dir="${17}"
simulated_tgfm_results_dir="${18}"

source ~/.bash_profile
module load R/4.0.1
echo "Simulation"$simulation_number
date

mkdir ${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_learned_gene_models_base_dir}"simulation_"${simulation_number}"/"

if false; then
n_causal_eqtls_per_gene_arr=( "2" "4" "6" "8" "10" )
n_causal_eqtls_per_gene_arr=( "8" "10" )

for n_causal_eqtls_per_gene in "${n_causal_eqtls_per_gene_arr[@]}"
do


simulation_name_string=${simulation_name_string_orig}"_n_eqtl_per_gene_"${n_causal_eqtls_per_gene}

echo $simulation_name_string

#######################################################
# Step 1: Randomly select chromosome windows to run tgfm on
# This will affect which genes we have to create gene models for
#######################################################
echo "Simulation Step 1"
global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
if false; then
python3 randomly_select_windows_to_run_tgfm_on_in_simulation.py $global_window_file $simulation_window_list_file $simulation_number $num_windows_per_sim
fi

#######################################################
# Step 2: Extract revised gene position file (limited to genes in one of above-specified windows)
#######################################################
simulated_revised_gene_position_file=${simulated_learned_gene_models_dir}${simulation_name_string}"_revised_gene_position.txt"
echo "Simulation Step 2"
if false; then
python3 extract_revised_gene_position_file.py $simulation_window_list_file $simulated_gene_position_file $simulated_revised_gene_position_file
fi



#######################################################
# Step 3: Simulate gene expression and fit gene models
#######################################################
echo "Simulation Step 3"
if false; then
python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_revised_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $n_causal_eqtls_per_gene
fi


#######################################################
# Step 4: Simulate trait values
#######################################################
echo "Simulation Step 4"
if false; then
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_trait_dir $n_gwas_individuals
fi


#######################################################
# Step 5: Run GWAS on simulated trait on only snps in TGFM windows.
# Also computes in-sample LD for TGFM windows
#######################################################
echo "Simulation Step 5"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
if false; then
python3 run_gwas_on_simulated_trait_at_snps_in_tgfm_windows.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulation_window_list_file $simulated_gwas_dir $simulated_tgfm_input_data_dir
fi

#######################################################
# Step 6: Preprocess data for TGFM
#######################################################
echo "Simulation Step 6"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size="inf"
eqtl_type="susie_pmces"
if false; then
python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
fi




#######################################################
# Step 7: Run TGFM (need to test)
#######################################################
echo "Simulation Step 7"
source ~/.bash_profile
module load R/4.0.1

if false; then
eqtl_type="susie_pmces"
eqtl_sample_size="inf"
est_resid_var="False"
init_method="best"
ln_pi_method="uniform"
echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
fi

eqtl_type="susie_pmces"
eqtl_sample_size="inf"
est_resid_var="False"
init_method="susie_inf"
ln_pi_method="uniform"

echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
done
fi


###########################
# TGFM with null permutation for one run


n_causal_eqtls_per_gene="6"
simulation_name_string=${simulation_name_string_orig}"_n_eqtl_per_gene_"${n_causal_eqtls_per_gene}


eqtl_type="susie_pmces"
eqtl_sample_size="inf"
est_resid_var="False"
init_method="finemap"
ln_pi_method="uniform"

echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}







