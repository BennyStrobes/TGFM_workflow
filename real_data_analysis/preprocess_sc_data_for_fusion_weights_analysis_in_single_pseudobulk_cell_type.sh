#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)



pb_cell_type_name="$1"
sc_pseudobulk_expression_dir="$2"
sc_fusion_input_dir="$3"
num_parallel_jobs="$4"
processed_sc_genotype_dir="$5"



cell_type_fusion_input_dir=${sc_fusion_input_dir}$pb_cell_type_name"/"
mkdir $cell_type_fusion_input_dir

# Prepare expression and covariates for fusion
python3 preprocess_sc_data_for_fusion_weights_analysis_in_single_pseudobulk_cell_type.py $pb_cell_type_name $sc_pseudobulk_expression_dir $cell_type_fusion_input_dir $num_parallel_jobs


# Prepare genotypes for fusion
# Limit to samples from exclusively this cell type
plink --bfile $processed_sc_genotype_dir"plink_geno_hg38_ukbb_overlap" --keep $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_fusion_ready_cov" --indiv-sort f $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_fusion_ready_cov" --make-bed --keep-allele-order --out $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_genotype"



# Need to split into seperate chromosomes (this is for computational efficiency when running fusion)
for chrom_num in $(seq 1 22); do 
	plink --bfile $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_genotype" --chr ${chrom_num} --make-bed --keep-allele-order --out $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_genotype_chr_"${chrom_num}
done


# Delete unnecessary genotypes
rm $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_genotype.bed"
rm $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_genotype.bim"
rm $cell_type_fusion_input_dir"cell_types_"${pb_cell_type_name}"_genotype.fam"
