#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2GB                         # Memory total in MiB (for all cores)





chrom_num="$1"
gencode_gene_annotation_file="$2"
simulated_gene_position_file="$3"

source ~/.bash_profile


python3 prepare_simulated_gene_position_list.py $chrom_num $gencode_gene_annotation_file $simulated_gene_position_file

