#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in




trait_name="$1"
tgfm_stem="$2"
gtex_pseudotissue_file="$3"
processed_tgfm_input_stem="$4"
tgfm_input_summary_file="$5"
iterative_tgfm_prior_results_dir="$6"



source ~/.bash_profile
module load R/4.0.1

echo $trait_name
date



########################
# PMCES (Existing best approach)
#########################
# Variant-gene prior
init_ln_pi_method="uniform"
new_tgfm_stem=${tgfm_stem}"_susie_pmces_"${init_ln_pi_method}
version="pmces"
python3 learn_iterative_tgfm_component_prior_pip_level_bootstrapped.py $trait_name $new_tgfm_stem $version $processed_tgfm_input_stem $gtex_pseudotissue_file $tgfm_input_summary_file $iterative_tgfm_prior_results_dir



# Variant-gene prior
if false; then
init_ln_pi_method="uniform"
new_tgfm_stem=${tgfm_stem}"_susie_pmces_"${init_ln_pi_method}
version="pmces"
python3 learn_iterative_tgfm_component_prior_w_ard_pip_level_bootstrapped.py $trait_name $new_tgfm_stem $version $processed_tgfm_input_stem $gtex_pseudotissue_file $tgfm_input_summary_file $iterative_tgfm_prior_results_dir
fi




