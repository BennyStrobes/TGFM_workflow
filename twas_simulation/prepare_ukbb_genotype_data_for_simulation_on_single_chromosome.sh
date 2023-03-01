
#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)





ukbb_genotype_dir="$1"
processed_genotype_data_dir="$2"
chrom_num="$3"
n_gwas_individuals="$4"
ldsc_baseline_hg19_annotation_dir="$5"
kg_genotype_dir="$6"

source ~/.bash_profile

###############################
# Extract list of variants in ldsc baseline analysis
###############################
ldsc_annotation_rs_id_file=${processed_genotype_data_dir}"ldsc_annotation_rsids_chr"${chrom_num}".txt"
python3 extract_list_of_ldsc_annotation_rs_ids.py $ldsc_baseline_hg19_annotation_dir $chrom_num $kg_genotype_dir $ldsc_annotation_rs_id_file


###############################
# Filter UKBB genotype data to only include those variants in ldsc baseline analysis
###############################
plink2 --pfile ${ukbb_genotype_dir}"ukb_imp_chr"${chrom_num}"_v3" --hwe .01 --extract ${ldsc_annotation_rs_id_file} --maf .05 --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}


###############################
# extract lists of Individuals for each data set
###############################
gwas_individual_file=${processed_genotype_data_dir}"gwas_individuals.txt"
eqtl_individual_stem=${processed_genotype_data_dir}"eqtl_individuals_"
python3 extract_lists_of_simulated_individuals.py ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} $n_gwas_individuals $gwas_individual_file $eqtl_individual_stem


###############################
# Filter UKBB genotype data to only include individuals in simulated gwas data set 
###############################
plink --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${gwas_individual_file} --make-bed --fill-missing-a2 --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}

###############################
# Filter UKBB genotype data to only include individuals in eqtl gwas data set 
###############################
eqtl_sample_size="100"
plink --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --fill-missing-a2 --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="200"
plink --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --fill-missing-a2 --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="300"
plink --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --fill-missing-a2 --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="500"
plink --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --fill-missing-a2 --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="1000"
plink --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --fill-missing-a2 --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}


#########################
# Filter 1KG variants to only those in our analysis
##########################
plink --bfile ${kg_genotype_dir}"1000G.EUR.QC."${chrom_num} --extract ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}".bim" --threads 1 --make-bed --keep-allele-order --out ${processed_genotype_data_dir}"100G.EUR.QC.filtered."${chrom_num}


#########################
# Remove unnecessary plink files
##########################
rm ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}*










