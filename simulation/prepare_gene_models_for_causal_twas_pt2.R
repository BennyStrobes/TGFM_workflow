args = commandArgs(trailingOnly=TRUE)
library(ctwas)
options(warn=1)



create_fusion_rds_object_for_gene_tissue_model <- function(snp_gene_effects, snps, RDat_output_file) {
	# Extract rsids
	rsids <- as.character(snps$V2)

	### 1. Create wgt.matrix
	wgt.matrix <- as.matrix(as.numeric(snp_gene_effects))
	colnames(wgt.matrix) = c("lasso")
	rownames(wgt.matrix) = rsids

	### 2. Create PLACEHOLDER cv.performance matrix (note this is not real data but doesnt matter for downstream purposes)
	cv.performance <- as.matrix(c(1.0, .00001))
	colnames(cv.performance) = c("lasso")
	rownames(cv.performance) = c("rsq", "pval")

	### 3. Create PLACEHOLDER hsq vector (note this is not real data but doesnt matter for downstream purposes)
	hsq <- as.numeric(c(.2, .05))

	### 4. Create PLACEHOLDER hsq.pv vector (note this is not real data but doesnt matter for downstream purposes)
	hsq.pv <- as.numeric(c(.00001))

	### 5. Create PLACEHOLDER N.tot vector (note this is not real data but doesnt matter for downstream purposes)
	N.tot <- as.numeric(c(100))

    # Alter column names of snps
	colnames(snps) <- c("chrom", "id", "cm", "pos", "alt", "ref")

	# Save to output
	save(wgt.matrix, snps , hsq, hsq.pv, N.tot , file = RDat_output_file )
}











#####################
# Command line args
#####################
tmp_pos_file = args[1]
gene_model_output_root = args[2]



# Load in temporary pos file
tmp_pos_df <- read.table(tmp_pos_file, header=TRUE, sep="\t")



# Loop through gene-tissue pairs
n_gt_pairs <- dim(tmp_pos_df)[1]

for (gt_pair_iter in 1:n_gt_pairs) {
	# Extract file-stem for this gene-tissue pair
	wgt_stem <- as.character(tmp_pos_df$WGT)[gt_pair_iter]
	# Extract relevent input files
	snp_gene_effects_input_file <- paste0(gene_model_output_root, wgt_stem, "_snp_gene_effects.txt")
	gene_bim_file <- paste0(gene_model_output_root, wgt_stem, "_snp_bim.txt")

	# Load in data
	snp_gene_effects <- read.table(snp_gene_effects_input_file, header=FALSE)$V1
	gene_bim <- read.table(gene_bim_file, header=FALSE, sep="\t")


	# Now create FUSION RDS obs for this gene
	RDat_output_file <- paste0(gene_model_output_root, wgt_stem, ".wgt.RDat")
	create_fusion_rds_object_for_gene_tissue_model(snp_gene_effects, gene_bim, RDat_output_file)
}