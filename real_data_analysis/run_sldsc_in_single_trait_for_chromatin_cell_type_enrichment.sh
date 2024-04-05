#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)



trait_name="$1"
ldsc_code_dir="$2"
full_sumstat_dir="$3"
ldsc_baseline_ld_hg19_annotation_dir="$4"
ldsc_cell_type_group_hg19_annotation_dir="$5"
ref_1kg_hg19_genotype_dir="$6"
sldsc_h19_weights_dir="$7"
chromatin_cell_type_group_ldsc_dir="$8"

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12


echo $trait_name

trait_file=$full_sumstat_dir"UKB_460K."${trait_name}".sumstats"

echo "BASELINELD"
# Run baseline_LD
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_hg19_annotation_dir}"baselineLD." --w-ld-chr ${sldsc_h19_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ref_1kg_hg19_genotype_dir}"1000G.EUR.QC." --out ${chromatin_cell_type_group_ldsc_dir}${trait_name}"_sldsc_baselineLD"

for cell_type_group_num in $(seq 1 10); do 
	echo "CELL_TYPE_GROUP"${cell_type_group_num}
	python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_cell_type_group_hg19_annotation_dir}"cell_type_group."${cell_type_group_num}".,"${ldsc_baseline_ld_hg19_annotation_dir}"baselineLD." --w-ld-chr ${sldsc_h19_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ref_1kg_hg19_genotype_dir}"1000G.EUR.QC." --out ${chromatin_cell_type_group_ldsc_dir}${trait_name}"_sldsc_cell_type_group_"${cell_type_group_num}"_baselineLD"
done