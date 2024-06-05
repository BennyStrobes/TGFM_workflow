#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-50:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=200GB                         # Memory total in MiB (for all cores)

# First three parts ran at 200GB



ukbb_genotype_dir="$1"
processed_genotype_data_root_dir="$2"
chrom_num="$3"
n_gwas_individuals="$4"
ldsc_baseline_hg19_annotation_dir="$5"
kg_genotype_dir="$6"
processed_ctwas_genotype_data_dir="$7"
processed_jlim_genotype_data_dir="$8"
simulated_gene_position_file="$9"

source ~/.bash_profile

###############################
# Extract list of variants in ldsc baseline analysis
###############################
ldsc_annotation_rs_id_file=${processed_genotype_data_root_dir}"ldsc_annotation_rsids_chr"${chrom_num}".txt"
if false; then
python3 extract_list_of_ldsc_annotation_rs_ids.py $ldsc_baseline_hg19_annotation_dir $chrom_num $kg_genotype_dir $ldsc_annotation_rs_id_file
fi

###############################
# Make genotype subdirectory for this gwas sampel size
###############################
processed_genotype_data_dir=${processed_genotype_data_root_dir}"gwas_sample_size_"${n_gwas_individuals}"/"
mkdir $processed_genotype_data_dir
echo $processed_genotype_data_dir

###############################
# Filter UKBB genotype data to only include those variants in ldsc baseline analysis
###############################
if false; then
plink2 \
    --bgen /n/groups/price/UKBiobank/download_500K/ukb_imp_chr"${chrom_num}"_v3.bgen ref-unknown\
    --sample /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample\
    --keep /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/unrelated_337K.txt\
    --extract /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/snp_info/snp_list_chr${chrom_num}.MAF_001_INFO_06.txt\
    --rm-dup force-first\
    --maj-ref\
    --geno 0.1\
    --maf 0.001\
    --hwe 1e-50\
    --make-pgen \
    --threads 1\
    --out ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper"

plink2 --pfile ${processed_genotype_data_dir}"ukb_imp_chr"${chrom_num}"_tmper" --hwe .01 --extract ${ldsc_annotation_rs_id_file} --maf .05 --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}
fi


###############################
# extract lists of Individuals for each data set
###############################
gwas_individual_file=${processed_genotype_data_dir}"gwas_individuals.txt"
eqtl_individual_stem=${processed_genotype_data_dir}"eqtl_individuals_"
ref_genotype_individual_file=${processed_genotype_data_dir}"ref_genotype_individuals.txt"
if false; then
python3 extract_lists_of_simulated_individuals.py ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} $n_gwas_individuals $gwas_individual_file $eqtl_individual_stem $ref_genotype_individual_file
fi

if false; then
###############################
# Filter UKBB genotype data to only include individuals in simulated gwas data set 
###############################
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${gwas_individual_file} --make-bed --keep-allele-order --threads 1 --out ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}

###############################
# Filter UKBB genotype data to only include individuals in eqtl gwas data set 
###############################
eqtl_sample_size="100"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="200"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="300"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="500"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}

eqtl_sample_size="1000"
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${eqtl_individual_stem}${eqtl_sample_size}".txt" --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_eqtl_"${eqtl_sample_size}"_data_"${chrom_num}
fi

###############################
# Filter UKBB genotype data to only include individuals in reference genotype data
###############################
if false; then
plink2 --bfile ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num} --keep ${ref_genotype_individual_file} --make-bed --threads 1 --keep-allele-order --out ${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num}
fi


#########################
# Filter 1KG variants to only those in our analysis
##########################
if false; then
plink2 --bfile ${kg_genotype_dir}"1000G.EUR.QC."${chrom_num} --extract ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}".bim" --threads 1 --make-bed --keep-allele-order --out ${processed_genotype_data_dir}"100G.EUR.QC.filtered."${chrom_num}


#########################
# Remove unnecessary plink files
##########################
rm ${processed_genotype_data_dir}"ukb_imp_chr_"${chrom_num}*
fi


#########################
# Extract genomic annotations for our variants
##########################
genomic_annotation_file=${processed_genotype_data_dir}"baseline."${chrom_num}".annot"
if false; then
python3 extract_genomic_annotations_for_simulation_variants.py ${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}".bim" ${processed_genotype_data_dir}"100G.EUR.QC.filtered."${chrom_num}".bim" $ldsc_baseline_hg19_annotation_dir"baseline."${chrom_num}".annot.gz" $genomic_annotation_file
fi

#########################
# Get reference genotype data for JLIM
##########################
#************NEEDTOCHANGE*********************************
tmp_reference_genotype_stem=${processed_genotype_data_dir}"simulated_reference_genotype_data_"${chrom_num}
tmp_reference_genotype_stem=${processed_genotype_data_dir}"simulated_eqtl_500_data_"${chrom_num}
#************NEEDTOCHANGE*********************************
source /home/bes710/.bash_profile
module load R/4.0.1
jlim_window_summary_file=${processed_genotype_data_dir}"jlim_window_summary_file.txt"
sh get_reference_genotype_data_at_each_gene_for_jlim.sh $tmp_reference_genotype_stem $simulated_gene_position_file $processed_jlim_genotype_data_dir $jlim_window_summary_file



source /home/bes710/.bash_profile
#########################
# Create 3MB chromosome windows
##########################
if false; then
window_size_mb="3" # In MB
reference_bim=${processed_genotype_data_dir}"simulated_gwas_data_"${chrom_num}".bim"  # Used to know where end and start on this chromosome
window_file=${processed_genotype_data_dir}"chromosome_"${chrom_num}"_windows_"${window_size_mb}"_mb.txt"
python3 create_chromosome_windows.py $window_size_mb $reference_bim $window_file ${chrom_num}

python3 generate_gwas_in_sample_ld_for_all_windows.py $chrom_num $window_file ${processed_genotype_data_dir}
echo "Done with TGFM windows"



#########################
# Create LD windows for c-TWAS
##########################
module load gcc/9.2.0
module load R/4.3.1
ctwas_regions_file=${processed_ctwas_genotype_data_dir}"ctwas_quasi_independent_regions.bed"
Rscript generate_quasi_independent_ld_windows_for_ctwas.R $chrom_num $processed_ctwas_genotype_data_dir $ctwas_regions_file


source ~/.bash_profile
python3 generate_gwas_in_sample_ld_for_all_ctwas_windows.py $chrom_num $ctwas_regions_file $processed_genotype_data_dir $processed_ctwas_genotype_data_dir
fi


