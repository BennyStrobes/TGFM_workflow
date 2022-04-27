#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=40GB                         # Memory total in MiB (for all cores)


window_file="$1"
output_dir="$2"
job_number="$3"
total_jobs="$4"


source ~/.bash_profile

module load R/4.0.1

echo $job_number

Rscript run_susie_genome_wide_in_parallel.R $window_file $output_dir $job_number $total_jobs
