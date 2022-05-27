#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1


gene_file="$1"
gtex_susie_results_dir="$2"


# Number of parallel jobs
total_jobs="100"



#########NOTE: OUTPUT FILE CURRENTLY SAYS NO_AMBIGUOUS_VARIANTS. SHOUULD BE FIXED
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
	echo $job_number
	sbatch run_susie_in_parallel.sh $gene_file $gtex_susie_results_dir $job_number $total_jobs
done
fi
