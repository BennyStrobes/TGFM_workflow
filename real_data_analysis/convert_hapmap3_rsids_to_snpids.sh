#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)





rsid_file="$1"
ukbb_hg38_genotype_dir="$2"
snpid_file="$3"

source ~/.bash_profile

python3 convert_hapmap3_rsids_to_snpids.py $rsid_file $ukbb_hg38_genotype_dir $snpid_file