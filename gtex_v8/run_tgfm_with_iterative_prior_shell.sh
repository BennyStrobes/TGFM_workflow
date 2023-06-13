#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-50:30                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=40GB                         # Memory total in MiB (for all cores)



trait_name="$1"
tgfm_input_summary_file="$2"
tgfm_output_stem="$3"
gtex_pseudotissue_file="$4"
job_number="$5"
num_jobs="$6"


# additional parameters
init_method="best"
est_resid_var="False"

source ~/.bash_profile
module load R/4.0.1

echo $trait_name
echo $job_number

date

echo "CURRENTLY RUNNING HACKY VERSION"


########################
# Run TGFM-SAMPLER
#########################
# Variant-gene prior
ln_pi_method="uniform_pmces_iterative_variant_gene_tissue_pip_level_sampler"
echo $ln_pi_method
new_tgfm_output_stem=${tgfm_output_stem}"_susie_sampler_"${ln_pi_method}
if false; then
python3 run_tgfm_sampler.py ${trait_name} ${tgfm_input_summary_file} ${new_tgfm_output_stem} ${job_number} ${num_jobs} ${init_method} ${est_resid_var} ${ln_pi_method} ${gtex_pseudotissue_file}
fi
python3 run_tgfm_sampler_hacky.py ${trait_name} ${tgfm_input_summary_file} ${new_tgfm_output_stem} ${job_number} ${num_jobs} ${init_method} ${est_resid_var} ${ln_pi_method} ${gtex_pseudotissue_file}


date
if false; then
# Variant-gene prior
ln_pi_method="uniform_pmces_iterative_variant_gene_tissue_pip_level_pmces"
echo $ln_pi_method
new_tgfm_output_stem=${tgfm_output_stem}"_susie_sampler_"${ln_pi_method}
python3 run_tgfm_sampler.py ${trait_name} ${tgfm_input_summary_file} ${new_tgfm_output_stem} ${job_number} ${num_jobs} ${init_method} ${est_resid_var} ${ln_pi_method} ${gtex_pseudotissue_file}
fi






