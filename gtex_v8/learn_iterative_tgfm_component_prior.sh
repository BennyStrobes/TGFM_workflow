#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




trait_name="$1"
tgfm_stem="$2"
gtex_pseudotissue_file="$3"
processed_tgfm_input_stem="$4"
tgfm_input_summary_file="$5"



source ~/.bash_profile
module load R/4.0.1

echo $trait_name
date

########################
# PMCES
#########################
# Variant-gene prior
init_ln_pi_method="variant_gene"
new_tgfm_stem=${tgfm_stem}"_susie_pmces_"${init_ln_pi_method}
version="pmces"
python3 learn_iterative_tgfm_component_prior.py $trait_name $new_tgfm_stem $version $processed_tgfm_input_stem $gtex_pseudotissue_file $tgfm_input_summary_file

date

########################
# Sampler
#########################
init_ln_pi_method="variant_gene"
new_tgfm_stem=${tgfm_stem}"_susie_sampler_"${init_ln_pi_method}
version="sampler"
python3 learn_iterative_tgfm_component_prior.py $trait_name $new_tgfm_stem $version $processed_tgfm_input_stem $gtex_pseudotissue_file $tgfm_input_summary_file
date
