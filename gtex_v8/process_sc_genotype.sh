#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --nodes=1



processed_genotype_dir="$1"
genotype_data_dir="$2"
filtered_sample_info_file="$3"
liftover_directory="$4"
ukbb_preprocessed_for_genome_wide_susie_dir="$5"


plink_format_filter_sample_info_file=$processed_genotype_dir"plink_formatted_filtered_sample_info.txt"

python3 extract_ordered_individuals_in_plink_format.py $filtered_sample_info_file $plink_format_filter_sample_info_file


plink2 --pfile $genotype_data_dir"inds_v2_header.maf10" --threads 1 --keep $plink_format_filter_sample_info_file --indiv-sort f $plink_format_filter_sample_info_file --maf .1 --make-pgen --out $processed_genotype_dir"inds_v2_sample_filter"


plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --threads 1 --freq --out $processed_genotype_dir"inds_v2_sample_filter_frq"

plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --indep-pairwise 200 50 0.25 --out $processed_genotype_dir"inds_v2_sample_filter_independent_snps"

plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --extract $processed_genotype_dir"inds_v2_sample_filter_independent_snps.prune.in" --pca 3 --out $processed_genotype_dir"inds_v2_sample_filter_ind_snp_pcs"



plink2 --pfile $processed_genotype_dir"inds_v2_sample_filter" --threads 1 --make-bed --out $processed_genotype_dir"plink_geno"


python3 liftover_sc_genotype_data_to_hg38.py $processed_genotype_dir"plink_geno.bim" $processed_genotype_dir"plink_geno_hg38" $liftover_directory



ukbb_sc_pbmc_formatted_bim=$processed_genotype_dir"ukbb_snps_sc_pbmc_formatted.bim"
python3 convert_ukbb_sumstat_file_to_sc_pbmc_snp_id_format.py $ukbb_preprocessed_for_genome_wide_susie_dir $ukbb_sc_pbmc_formatted_bim


# Hacky fix
mv $processed_genotype_dir"plink_geno.bed" $processed_genotype_dir"plink_geno_hg38.bed" 
mv $processed_genotype_dir"plink_geno.fam" $processed_genotype_dir"plink_geno_hg38.fam"


plink --bfile $processed_genotype_dir"plink_geno_hg38" --allow-extra-chr --extract $ukbb_sc_pbmc_formatted_bim --make-bed --keep-allele-order --out $processed_genotype_dir"plink_geno_hg38_ukbb_overlap"




