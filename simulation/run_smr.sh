#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)





simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
simulated_gene_expression_dir="$8"
simulated_gwas_dir="$9"
ge_h2="${10}"
eqtl_architecture="${11}"
per_element_heritability="${12}"
total_heritability="${13}"
fraction_expression_mediated_heritability="${14}"
gene_trait_architecture="${15}"
simulated_trait_dir="${16}"
simulated_smr_data_dir="${17}"
simulated_tmp_smr_eqtl_dir="${18}"
simulated_smr_results_dir="${19}"
smr_source_code_dir="${20}"



date

echo $simulation_name_string
echo $eqtl_sample_size


source /home/bes710/.bash_profile

# Make learned gene models output root seperate for each simulatino
mkdir ${simulated_tmp_smr_eqtl_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_tmp_smr_eqtl_dir}"simulation_"${simulation_number}"/"



eqtl_sample_size_arr=( "realistic" "300" "500" "1000")
for eqtl_sample_size in "${eqtl_sample_size_arr[@]}"
do

	echo "eQTL sample size: "${eqtl_sample_size}



	#######################################################
	# Step 1: Get eQTL summary stats for smr
	#######################################################
	module load R/4.0.1
	eqtl_flist_file=${simulated_learned_gene_models_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_gene_list.flist"
	python3 get_eqtl_sumstats_for_smr.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $ge_h2 $eqtl_architecture $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $gene_trait_architecture $eqtl_sample_size $eqtl_flist_file




	#######################################################
	# Step 2: User SMR to make besd file (sparse data structure) storing eqtl data
	#######################################################
	besd_file_stem=${simulated_smr_data_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}".mybesd"
	${smr_source_code_dir}"smr" --eqtl-flist ${eqtl_flist_file} --make-besd --out $besd_file_stem



	#######################################################
	# Step 3: Run SMR
	#######################################################
	# GWAS summary statistics (produced for smr, not jlim)
	gwas_summary_file=${simulated_gwas_dir}${simulation_name_string}"_merged_unstandardized_jlim_gwas_summary_stats.ma"
	# Reference genotype (using genotype from 1000 eqtl individauls)
	ref_geno_plink_stem=${processed_genotype_data_dir}"simulated_eqtl_1000_data_1"
	# SMR output
	smr_output_root=${simulated_smr_results_dir}${simulation_name_string}"_eqtl_ss_"${eqtl_sample_size}"_smr_res"

	# Run SMR
	${smr_source_code_dir}"smr" --bfile ${ref_geno_plink_stem} --gwas-summary ${gwas_summary_file} --beqtl-summary ${besd_file_stem} --out ${smr_output_root}




	#######################################################
	# Step 4: Delete SMR eqtl summary statistics
	#######################################################
	python3 delete_smr_eqtl_summary_statistics.py $eqtl_flist_file $besd_file_stem

done

# TO DO: 
# 2. Allow analysis at realistic QTL sample size
# 3. Add code to delete gwas region files
# 4. Run gwas regions


date