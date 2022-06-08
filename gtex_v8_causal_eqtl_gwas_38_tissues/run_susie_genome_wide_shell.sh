#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1


window_file="$1"
susie_results_dir="$2"


# Number of parallel jobs
total_jobs="100"

if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
	echo $job_number
	sbatch run_susie_genome_wide_in_parallel.sh $window_file $susie_results_dir $job_number $total_jobs
done
fi
