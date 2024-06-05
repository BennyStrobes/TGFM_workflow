#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-55:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=12GB                         # Memory total in MiB (for all cores)




simulation_number="$1"
chrom_num="$2"
cis_window="$3"
n_gwas_individuals="$4"
simulation_name_string="$5"
simulated_gene_position_file="$6"
processed_genotype_data_dir="$7"
simulated_gene_expression_dir="$8"
simulated_gwas_dir="$9"
eqtl_sample_size="${10}"
jlim_source_code_dir="${11}"
processed_jlim_genotype_data_dir="${12}"
simulated_jlim_data_dir="${13}"
simulated_jlim_results_dir="${14}"
simulated_jlim_smr_eqtl_sumstats_dir="${15}"
ge_h2="${16}"
eqtl_architecture="${17}"
per_element_heritability="${18}"
total_heritability="${19}"
fraction_expression_mediated_heritability="${20}"
gene_trait_architecture="${21}"
jlim_ref_geno_dir="${22}"

source /home/bes710/.bash_profile
module load R/4.0.1

echo $simulation_name_string
echo $eqtl_sample_size






# Make learned gene models output root seperate for each simulatino
mkdir ${simulated_jlim_smr_eqtl_sumstats_dir}"simulation_"${simulation_number}
simulated_learned_gene_models_dir=${simulated_jlim_smr_eqtl_sumstats_dir}"simulation_"${simulation_number}"/"


global_window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_3_mb.txt"




#######################################################
# Step 1: Get eQTL Summary stats for jlim and smr
#######################################################
python3 get_eqtl_sumstats_for_jlim_and_smr.py $simulation_number $chrom_num $cis_window $simulated_gene_position_file $simulated_gene_expression_dir $simulated_learned_gene_models_dir $simulation_name_string $processed_genotype_data_dir $ge_h2 $eqtl_architecture $per_element_heritability $total_heritability $fraction_expression_mediated_heritability $gene_trait_architecture $eqtl_sample_size




#######################################################
# Step 2: Run JLIM for each gene-tissue pair
#######################################################
python3 run_jlim_for_each_gene_tissue_pair.py $simulation_number $chrom_num $eqtl_sample_size $simulation_name_string $processed_genotype_data_dir $processed_jlim_genotype_data_dir $jlim_ref_geno_dir $jlim_source_code_dir $simulated_learned_gene_models_dir $simulated_jlim_data_dir $simulated_jlim_results_dir $simulated_gwas_dir






# Source code example for JLIM 2.5.0
if false; then
sh ${jlim_source_code_dir}"run_jlim.sh" --maintr-file ${jlim_source_code_dir}"examples/MS/MS.1.160697074.160933065.txt" \
--sectr-file ${jlim_source_code_dir}"examples/LCL/locus.1.160697074.160933065/LCL.SLAMF7_ENSG00000026751.11.assoc.linear.gz" \
--ref-ld "/n/groups/price/ben/tools/JLIM/refld.1kg.nfe.b37/" \
--index-snp 1:160711804 \
--output-file ${jlim_source_code_dir}"examples/MS-Geuvadis_LCL_SLAMF7.out" \
--sectr-sample-size 278
fi