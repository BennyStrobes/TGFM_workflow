#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in


tgfm_results_dir="$1"
gene_type="$2"
trait_names_file="3"
gtex_gene_annotation_file="$4"
gtex_susie_gene_models_dir="$5"
preprocessed_tgfm_data_dir="$6"
tgfm_organized_results_dir="$7"
hg38_gene_annotation_file="$8"
hg38_ensamble_trascript_to_entrez_id_file="$9"
biological_pathway_gene_set_file="${10}"
biological_pathway_enrichment_dir="${11}"



ensamble_biological_pathway_gene_set_file=${biological_pathway_enrichment_dir}"biological_pathway_ensamble_id_gene_sets.txt"
if false; then
python3 convert_biological_pathway_gene_set_file_from_entrez_id_to_ensamble_ids.py $hg38_gene_annotation_file $hg38_ensamble_trascript_to_entrez_id_file $biological_pathway_gene_set_file $ensamble_biological_pathway_gene_set_file
fi

if false; then
python3 run_biological_pathway_gene_set_enrichment_analysis.py $tgfm_results_dir $gene_type $trait_names_file $gtex_gene_annotation_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $ensamble_biological_pathway_gene_set_file ${biological_pathway_enrichment_dir} $tgfm_organized_results_dir
fi

if false; then
python3 run_biological_pathway_gene_set_enrichment_analysis2.py $tgfm_results_dir $gene_type $trait_names_file $gtex_gene_annotation_file $gtex_susie_gene_models_dir $preprocessed_tgfm_data_dir $ensamble_biological_pathway_gene_set_file ${biological_pathway_enrichment_dir} $tgfm_organized_results_dir
fi





if false; then
python3 organize_biological_pathway_gene_set_enrichment_analysis_results.py ${biological_pathway_enrichment_dir}
fi