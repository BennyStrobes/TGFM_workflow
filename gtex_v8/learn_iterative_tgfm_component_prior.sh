#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)




trait_name="$1"
tgfm_stem="$2"
gtex_pseudotissue_file="$3"
processed_tgfm_input_stem="$4"



source ~/.bash_profile
module load R/4.0.1

echo $trait_name


########################
# PMCES
#########################
if false; then
# Variant-gene prior
init_ln_pi_method="variant_gene"
new_tgfm_stem=${tgfm_stem}"_susie_pmces_"${init_ln_pi_method}
version="pmces"
python3 learn_iterative_tgfm_component_prior.py $trait_name $new_tgfm_stem $version $processed_tgfm_input_stem $gtex_pseudotissue_file
fi

########################
# Sampler
#########################
init_ln_pi_method="variant_gene"
new_tgfm_stem=${tgfm_stem}"_susie_sampler_"${init_ln_pi_method}
version="sampler"
if false; then
python3 learn_iterative_tgfm_component_prior.py $trait_name $new_tgfm_stem $version $processed_tgfm_input_stem $gtex_pseudotissue_file
fi