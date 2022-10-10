#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)


gene_file="$1"
trait_name="$2"
gtex_tissue_file="$3"
coloc_output_dir="$4"
job_number="$5"
total_jobs="$6"

source ~/.bash_profile

echo $trait_name
echo $job_number


python3 run_adaptive_prior_coloc_first_pass.py $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $job_number $total_jobs

