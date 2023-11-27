#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



source ~/.bash_profile

liftover_directory="$1"
ukbb_sumstats_hg19_dir="$2"
ukbb_sumstats_hg38_dir="$3"


python3 liftover_ukbb_summary_statistics_from_hg19_to_hg38.py $liftover_directory $ukbb_sumstats_hg19_dir $ukbb_sumstats_hg38_dir


