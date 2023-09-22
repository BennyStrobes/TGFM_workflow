#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




processed_sc_expression_dir="$1"
input_h5py_file="$2"
processed_genotype_dir="$3"
pseudobulk_expression_dir="$4"
gene_annotation_file="$5"
hg38_gene_annotation_file="$6"
gtex_gene_list="$7"

echo $hg38_gene_annotation_file
echo $gtex_gene_list

python3 generate_pseudobulk_expression.py $processed_sc_expression_dir $input_h5py_file $processed_genotype_dir $pseudobulk_expression_dir $gene_annotation_file


python3 get_sc_pseudobulk_gene_tss.py $pseudobulk_expression_dir $hg38_gene_annotation_file $gtex_gene_list
