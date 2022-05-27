#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)




trait_name="$1"
gtex_tissue_file="$2"
coloc_input_dir="$3"
coloc_output_dir="$4"


# Number of parallel jobs
total_jobs="10"

gene_file=$coloc_input_dir$trait_name"_processed_gene_list.txt"
if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	sbatch run_adaptive_prior_coloc_first_pass.sh $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $job_number $total_jobs
done
fi


pseudocount="100"
if false; then
sbatch run_adaptive_prior_coloc_second_pass.sh $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $pseudocount
fi

if false; then
for chrom_num in $(seq 1 22); do
	echo $chrom_num
	sh organize_coloc_results_for_prs.sh $gtex_tissue_file $trait_name $coloc_input_dir $coloc_output_dir $chrom_num
done
fi
if false; then
python3 print_number_of_coloc_components_per_tissue.py $gtex_tissue_file $trait_name $coloc_output_dir
fi
