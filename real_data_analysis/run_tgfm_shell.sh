#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-25:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=85GB                         # Memory total in MiB (for all cores)

# Most at 60, 4's and 6's at 120

trait_name="$1"
tgfm_input_summary_file="$2"
tgfm_output_stem="$3"
gtex_pseudotissue_file="$4"
job_number="$5"
num_jobs="$6"
ignore_tissues="$7"


# additional parameters
init_method="best"
est_resid_var="False"

source ~/.bash_profile
module load R/4.0.1

echo $trait_name
echo $job_number
date


# UNIFORM PRIOR
echo "Uniform"
########################
# Run TGFM-PMCES
#########################
# Variant-gene-tissue prior from sparse model
ln_pi_method="uniform"
echo $ln_pi_method
new_tgfm_output_stem=${tgfm_output_stem}"_susie_pmces_"${ln_pi_method}
python3 run_tgfm_pmces.py ${trait_name} ${tgfm_input_summary_file} ${new_tgfm_output_stem} ${job_number} ${num_jobs} ${init_method} ${est_resid_var} ${ln_pi_method} ${gtex_pseudotissue_file} ${ignore_tissues}
date
























