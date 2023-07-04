args = commandArgs(trailingOnly=TRUE)
options(warn=1)







#############################
# Load in command line args
#############################
simulation_number <- args[1]
chrom_num <- args[2]
simulation_name_string <- args[3]
eqtl_sample_size <- args[4]
focus_gene_models_dir <- args[5]
focus_gene_model_python_file <- args[6]
focus_gene_model_pos_file <- args[7]

output_root <- paste0(focus_gene_models_dir, simulation_name_string, "_eqtlss_", eqtl_sample_size)
short_output_root <- paste0(simulation_name_string, "_eqtlss_", eqtl_sample_size)

# Load in input data
input_gene_model_df <- read.table(focus_gene_model_python_file, header=TRUE, sep="\t")

# Get number of genes
n_genes = dim(input_gene_model_df)[1]

sink(focus_gene_model_pos_file)
header_line <- "WGT\tID\tCHR\tP0\tP1\n"
cat(header_line)


# Loop through genes
for (gene_iter in 1:n_genes) {
	gene_name <- input_gene_model_df$gene_tissue[gene_iter]
	chrom_num <- input_gene_model_df$chr[gene_iter]
	tss <- input_gene_model_df$tss[gene_iter]
	susie_weights_file <- input_gene_model_df$wgt.matrix_file[gene_iter]
	snps_file <- input_gene_model_df$snps_file[gene_iter]

	# Load in snps data
	snps <- read.table(snps_file, header=FALSE, sep="\t")

	# Get rs-ids from snps_df
	rs_ids = snps$V2

	# Load in susie weights
	susie_weights <- read.table(susie_weights_file, header=FALSE)
	# Get weight matrix
	wgt.matrix = matrix(0,nrow=length(rs_ids),ncol=1)
	colnames(wgt.matrix) = c("susie")
	rownames(wgt.matrix) = rs_ids
	wgt.matrix[,1] = susie_weights$V1

	# Create placeholder cv matrix
	cv.performance = matrix(0,nrow=2,ncol=1)
	rownames(cv.performance) = c("rsq","pval")
	colnames(cv.performance) = c("susie")
	cv.performance[1,1] = 1.0
	cv.performance[2,1] = 0.0

	# Create placeholder hsq vector
	hsq <- numeric(2)
	hsq[1] = .2
	hsq[2] = .001

	# Create placeholder hsq pv
	hsq.pv <- numeric(0)

	# Save FUSION Rdata
	Rdata_file_name = paste( output_root, "_", gene_name , ".wgt.RDat" , sep='' )
	short_Rdata_file_name = paste0(short_output_root, "_", gene_name , ".wgt.RDat")
	save(wgt.matrix , snps , cv.performance , hsq, hsq.pv , file = Rdata_file_name )

	# Print to output file
	new_line <- paste0(short_Rdata_file_name, "\t", gene_name, "\t", chrom_num, "\t", tss, "\t", tss, "\n")
	cat(new_line)
}
sink()

