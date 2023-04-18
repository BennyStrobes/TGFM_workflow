#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-40:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=22GB                         # Memory total in MiB (for all cores)



simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
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










#######################################################
# Step 1: Simulate gene expression 
#######################################################
echo "Simulation Step 1"
python3 simulate_causal_eqtl_effect_sizes.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir



#######################################################
# Step 2: Simulate trait values
#######################################################
echo "Simulation Step 2"
python3 simulate_trait_values.py $simulation_number $chrom_num $cis_window $simulated_gene_expression_dir $simulation_name_string $processed_genotype_data_dir $ldsc_real_data_results_dir $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $simulated_trait_dir $n_gwas_individuals



#######################################################
# Now Looop through number of causal genetic elements
#######################################################
n_causal_genetic_elements_arr=( "0_5" "6_10" "11_15" )


for n_causal_genetic_elements in "${n_causal_genetic_elements_arr[@]}"
do
	echo '############################################'
	echo ${n_causal_genetic_elements}" causal genetic elements"
	echo '############################################'

	#######################################################
	# Step 3: Randomly select chromosome windows to run tgfm on
	# This will affect which genes we have to create gene models for
	#######################################################
	echo "Simulation Step 3"
	global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"
	simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_tgfm_windows.txt"
	python3 randomly_select_windows_to_run_tgfm_on_in_simulation_based_on_n_causal_elements.py $global_window_file $simulation_window_list_file $simulation_number $num_windows_per_sim ${n_causal_genetic_elements} ${simulated_trait_dir} ${simulation_name_string} ${simulated_gene_position_file} ${processed_genotype_data_dir} ${chrom_num}


	#######################################################
	# Step 4: Extract revised gene position file (limited to genes in one of above-specified windows)
	#######################################################
	simulated_revised_gene_position_file=${simulated_learned_gene_models_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_revised_gene_position.txt"
	echo "Simulation Step 4"
	python3 extract_revised_gene_position_file.py $simulation_window_list_file $simulated_gene_position_file $simulated_revised_gene_position_file


	# FIT GENE MODELS
	#######################################################
	# Step 5: Simulate gene expression and fit gene models
	#######################################################
	echo "Simulation Step 5"
	python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_revised_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $n_causal_genetic_elements

	# RUN GWAS ON SIMULATED TRAIT
	# Step 6: Run GWAS on simulated trait on only snps in TGFM windows.
	# Also computes in-sample LD for TGFM windows
	echo "Simulation Step 6"
	simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_tgfm_windows.txt"
	python3 run_gwas_on_simulated_trait_at_snps_in_tgfm_windows.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulation_window_list_file $simulated_gwas_dir $simulated_tgfm_input_data_dir $n_causal_genetic_elements



	#######################################################
	# Step 6: Preprocess data for TGFM
	#######################################################
	echo "Simulation Step 6"
	simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_tgfm_windows.txt"
	annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
	eqtl_sample_size_arr=( "300" "500" "1000" "inf" )
	eqtl_type="susie_pmces"
	echo $eqtl_type
	for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
	do
		python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir $n_causal_genetic_elements
	done

	#######################################################
	# Step 7: Run TGFM
	#######################################################
	echo "Simulation Step 7"
	eqtl_sample_size_arr=( "300" "500" "1000" "inf" )
	ln_pi_method_arr=( "uniform" )
	init_method_arr=( "best" )
	eqtl_type="susie_pmces"
	for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
	do
		for ln_pi_method in "${ln_pi_method_arr[@]}"
		do
			for init_method in "${init_method_arr[@]}"
			do
			est_resid_var="False"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
			done
		done
	done



	#######################################################
	# Step 8: Preprocess data for bootstrapped-TGFM
	#######################################################
	echo "Simulation Step 8"
	simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_tgfm_windows.txt"
	annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
	eqtl_sample_size_arr=( "300" "500" "1000" )
	eqtl_type="susie_distr"
	echo $eqtl_type
	for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
	do
		echo $eqtl_sample_size
		python3 preprocess_data_for_bootstrapped_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir $n_causal_genetic_elements
	done


	#######################################################
	# Step 9: Run bootstrapped-TGFM
	#######################################################
	echo "Simulation Step 9"
	eqtl_sample_size_arr=( "300" "500" "1000" )
	ln_pi_method_arr=( "uniform" )
	init_method_arr=( "best" )
	eqtl_type="susie_distr"
	est_resid_var="False"
	for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
	do
		for ln_pi_method in "${ln_pi_method_arr[@]}"
		do
			for init_method in "${init_method_arr[@]}"
			do
				echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}
				tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_bootstrapped_tgfm_input_data_summary.txt"
				tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_"${n_causal_genetic_elements}"_eqtl_ss_"${eqtl_sample_size}"_bootstrapped_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
				python3 run_bootstrapped_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
			done
		done
	done


done




























# OLD. Now need to simulate first and extract windows in ranges

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
python3 simulate_gene_expression_and_fit_gene_model.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_revised_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir
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
if false; then
echo "Simulation Step 5"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
python3 run_gwas_on_simulated_trait_at_snps_in_tgfm_windows.py $simulation_number $chrom_num $simulation_name_string $processed_genotype_data_dir $simulated_trait_dir $simulation_window_list_file $simulated_gwas_dir $simulated_tgfm_input_data_dir
fi

#######################################################
# Step 6: Preprocess data for TGFM
#######################################################
echo "Simulation Step 6"
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "100" "300" "500" "1000" "inf" )
eqtl_type="susie_pmces"
if false; then
echo $eqtl_type
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done



simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "100" "300" "500" "1000" )
eqtl_type="susie_distr"
echo $eqtl_type
if false; then
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done
fi

simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "100" "300" "500" "1000" )
eqtl_type="fusion_lasso_pmces"
echo $eqtl_type
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done
fi

if false; then

simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "100" "300" "500" "1000" )
eqtl_type="marginal_distr"
echo $eqtl_type
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done
fi


#######################################################
# Step 7: Run TGFM (need to test)
#######################################################
echo "Simulation Step 7"
if false; then
source ~/.bash_profile
module load R/4.0.1
eqtl_sample_size_arr=( "100" "300" "500" "1000" )
ln_pi_method_arr=( "uniform" )
init_method_arr=( "best" )
eqtl_type="susie_distr"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	for ln_pi_method in "${ln_pi_method_arr[@]}"
	do
		for init_method in "${init_method_arr[@]}"
		do
			est_resid_var="False"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}

			est_resid_var="True"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
		done
	done
done



eqtl_sample_size_arr=( "100" "300" "500" "1000" )
ln_pi_method_arr=( "uniform" )
init_method_arr=( "best" )
eqtl_type="fusion_lasso_pmces"

for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	for ln_pi_method in "${ln_pi_method_arr[@]}"
	do
		for init_method in "${init_method_arr[@]}"
		do
			est_resid_var="False"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
			
			est_resid_var="True"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
		done
	done
done





eqtl_sample_size_arr=( "100" "300" "500" "1000" "inf" )
ln_pi_method_arr=( "uniform" )
init_method_arr=( "best" )
eqtl_type="susie_pmces"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	for ln_pi_method in "${ln_pi_method_arr[@]}"
	do
		for init_method in "${init_method_arr[@]}"
		do
			est_resid_var="False"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
			
			est_resid_var="True"
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
		done
	done
done
fi







if false; then
############################
# Fine-map using RSS-impute likelihood
############################

eqtl_type="marginal_distr"
eqtl_sample_size="500"
python3 preprocess_data_for_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir



eqtl_type="marginal_distr"
eqtl_sample_size="500"
est_resid_var="False"
init_method="best"
ln_pi_method="uniform"

echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}

echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_rss_impute_likelihood_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
python3 run_rss_imputation_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
fi








if false; then
eqtl_type="susie_pmces"
eqtl_sample_size="inf"
est_resid_var="False"
init_method="null_perm"
ln_pi_method="uniform"

echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}"_"${est_resid_var}
tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_tgfm_input_data_summary.txt"
tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
python3 run_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
fi


#######################################################
# Step 8: Preprocess data for bootstrapped-TGFM
#######################################################
echo "Simulation Step 8"
if false; then
simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "500" )
eqtl_type="sparse_marginal_distr"
echo $eqtl_type
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_bootstrapped_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done

eqtl_sample_size_arr=( "500" )
ln_pi_method_arr=( "uniform" )
init_method_arr=( "best" )
eqtl_type="sparse_marginal_distr"
est_resid_var="False"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	for ln_pi_method in "${ln_pi_method_arr[@]}"
	do
		for init_method in "${init_method_arr[@]}"
		do
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_bootstrapped_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_bootstrapped_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_bootstrapped_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
		done
	done
done
fi

if false; then
rm ${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_500_sparse_marginal_distr"*"npy"


rm ${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_500_susie_distr"*"npy"
rm ${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_500_marginal_dist"*"npy"

echo "start"

simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "500" )
eqtl_type="marginal_distr"
echo $eqtl_type
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_bootstrapped_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done

eqtl_sample_size_arr=( "500" )
ln_pi_method_arr=( "uniform" )
init_method_arr=( "best" )
eqtl_type="marginal_distr"
est_resid_var="False"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	for ln_pi_method in "${ln_pi_method_arr[@]}"
	do
		for init_method in "${init_method_arr[@]}"
		do
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_bootstrapped_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_bootstrapped_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_bootstrapped_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
		done
	done
done


rm ${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_500_marginal_dist"*"npy"


simulation_window_list_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_tgfm_windows.txt"
annotation_file=${processed_genotype_data_dir}baseline.${chrom_num}.annot
eqtl_sample_size_arr=( "500" )
eqtl_type="susie_distr"
echo $eqtl_type
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	python3 preprocess_data_for_bootstrapped_tgfm.py $simulation_number $chrom_num $simulation_name_string $n_gwas_individuals $eqtl_sample_size $simulation_window_list_file $annotation_file $simulated_gwas_dir $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulated_tgfm_input_data_dir $eqtl_type $processed_genotype_data_dir
done


eqtl_sample_size_arr=( "500" )
ln_pi_method_arr=( "uniform" )
init_method_arr=( "best" )
eqtl_type="susie_distr"
est_resid_var="False"
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do
	for ln_pi_method in "${ln_pi_method_arr[@]}"
	do
		for init_method in "${init_method_arr[@]}"
		do
			echo ${eqtl_sample_size}"_"${ln_pi_method}"_"${init_method}"_"${eqtl_type}
			tgfm_input_file=${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_"${eqtl_type}"_bootstrapped_tgfm_input_data_summary.txt"
			tgfm_output_stem=${simulated_tgfm_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_bootstrapped_"${eqtl_type}"_ln_pi_"${ln_pi_method}"_init_"${init_method}"_resid_var_"${est_resid_var}
			python3 run_bootstrapped_tgfm.py ${tgfm_input_file} ${tgfm_output_stem} ${ln_pi_method} ${init_method} ${est_resid_var}
		done
	done
done

rm ${simulated_tgfm_input_data_dir}${simulation_name_string}"_eqtl_ss_500_susie_distr"*"npy"

fi


