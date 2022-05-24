args = commandArgs(trailingOnly=TRUE)
library(hash)
library(susieR)
library('plink2R')
library("susieR")
options(warn=1)


allele.qc = function(a1,a2,ref1,ref2) {
        a1 = toupper(a1)
        a2 = toupper(a2)
        ref1 = toupper(ref1)
        ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="ZZZZZZZ")) # MY Hack to ensure no snps are thrown out b/c of commented out strand ambiguity below b/c not proficient in R
	#snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}

create_gene_tissue_df <- function(tissue_names, gtex_fusion_weights_dir, gtex_fusion_associations_dir, global_chrom_num, trait_name) {
	gene_arr <- c()
	tissue_arr <- c()
	wgt_file_arr <- c()
	z_score_arr <- c()
	pvalue_arr <- c()
	chrom_num_arr <- c()
	tss_arr <- c()
	top_model_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]

		# Tissue pos file (containins all genes with wgt files AND their TSS)
		# TURNS OUT THIS IS NOT NECESSARY
		tissue_pos_file <- paste0(gtex_fusion_weights_dir, tissue_name, "/", tissue_name, ".pos")
		tissue_pos_df <- read.table(tissue_pos_file, header=TRUE)  # Load in file
		tissue_pos_df <- tissue_pos_df[tissue_pos_df$CHR==global_chrom_num,]  # Filter to genes on this chromosome

		# Fusion association results file
		fusion_assoc_results_file <- paste0(gtex_fusion_associations_dir, tissue_name, ".", trait_name, "_", global_chrom_num, ".dat")
		fusion_assoc_results_df <- read.table(fusion_assoc_results_file, header=TRUE, sep="\t")
		ngenes <- dim(fusion_assoc_results_df)[1]


		for (gene_num in 1:ngenes) {
			gene_name <- fusion_assoc_results_df$ID[gene_num]
			twas_Z <- fusion_assoc_results_df$TWAS.Z[gene_num]
			twas_P <- fusion_assoc_results_df$TWAS.P[gene_num]
			p0 <- fusion_assoc_results_df$P0[gene_num]
			wgt_file <- fusion_assoc_results_df$FILE[gene_num]
			top_model <- fusion_assoc_results_df$MODEL[gene_num]

			gene_arr <- c(gene_arr, gene_name)
			tissue_arr <- c(tissue_arr, tissue_name)
			wgt_file_arr <- c(wgt_file_arr, wgt_file)
			z_score_arr <- c(z_score_arr, twas_Z)
			pvalue_arr <- c(pvalue_arr, twas_P)
			tss_arr <- c(tss_arr, p0)
			chrom_num_arr <- c(chrom_num_arr, global_chrom_num)
			top_model_arr <- c(top_model_arr, top_model)
		}
	}
	df <- data.frame(gene=gene_arr, tissue=tissue_arr, wgt_file=wgt_file_arr, twas_nominal_z=z_score_arr, twas_nominal_p=pvalue_arr, chrom_num=chrom_num_arr, gene_tss=tss_arr, top_model=top_model_arr)
	return(df)
}


######################
# Command line args
#######################
global_chrom_num = as.numeric(args[1])
trait_name = args[2]
ukbb_genome_wide_susie_organized_results_dir = args[3]
gtex_tissue_file = args[4]
gtex_fusion_weights_dir = args[5]
kg_genotype_dir = args[6]
gtex_fusion_associations_dir = args[7]
gwas_sample_size = as.numeric(args[8])
gtex_fusion_multivariate_associations_dir = args[9]



# Extract tissue names
tissue_names <- read.table(gtex_tissue_file, header=TRUE, sep="\t")$tissue_name 

# Create data frame of (tissue, gene) on this chromosome that we have wgt files for
fusion_gene_tissue_df <- create_gene_tissue_df(tissue_names, gtex_fusion_weights_dir, gtex_fusion_associations_dir, global_chrom_num, trait_name)


# Extract susie genome wide trait components on this chromosome
organized_susie_trait_component_file <- paste0(ukbb_genome_wide_susie_organized_results_dir, trait_name, "_organized_susie_components.txt")  # File name
organized_susie_trait_component_df <- read.table(organized_susie_trait_component_file, header=TRUE)  # load in file
organized_susie_trait_component_df <- organized_susie_trait_component_df[organized_susie_trait_component_df$chrom_num==global_chrom_num,]  # Filter to chromosomes on the chromosome of interest

# Load in sumstat file (necessary for variant orientation)
sumstat_file <- paste0(kg_genotype_dir, trait_name, "_fusion_processed_sumstats_hg38.txt")
sumstat = read.table(sumstat_file,head=T,as.is=T)

# Load in Ref LD data
ref_ld_chr = paste0(kg_genotype_dir, "1000G.EUR.gtex_formatted.hg38.", global_chrom_num)
genos = read_plink(ref_ld_chr,impute="avg")
#saveRDS(genos, file = "temp_genos.RDS")
#genos <- readRDS("temp_genos.RDS")




#################################
# Some genotype processing (just to match fusion as closely as possible)
#################################
# Match summary data to input, record NA where summary data is missing
m = match( genos$bim[,2] , sumstat$SNP )
sum.missing = is.na(m)
sumstat = sumstat[m,]
sumstat$SNP = genos$bim[,2]
sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]

# QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )

# Flip Z-scores for mismatching alleles
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
	genos$bim = genos$bim[qc$keep,]
	genos$bed = genos$bed[,qc$keep]
	sumstat = sumstat[qc$keep,]
}




twas_nominal_pvalue_table_names <- c()
twas_multivariate_susie_object_names <- c()
nominal_z_files <- c()
pred_expr_corr_files <- c()
gene_tissue_names_files <- c()


# Now begin our loop through Susie genome wide trait components
num_susie_trait_components <- dim(organized_susie_trait_component_df)[1]

for (trait_component_iter in 1:num_susie_trait_components) {
	# Extract relvent fields for this trait component
	trait_component_lead_variant <- organized_susie_trait_component_df$lead_variant_name[trait_component_iter]  # Name of lead variant for trait component
	trait_component_window_id <- organized_susie_trait_component_df$window_id[trait_component_iter]  # Name of window id defining trait component
	trait_component_lead_variant_pos <- organized_susie_trait_component_df$lead_variant_position[trait_component_iter]  # Position of lead variant for trait component
	component_name <- paste0(trait_name, "_", trait_component_lead_variant, "_", trait_component_window_id)
	print(component_name)

	# Extract gene-tissue pairs in cis window of this trait component
	cis_window_gene_indices <- abs(fusion_gene_tissue_df$gene_tss - trait_component_lead_variant_pos) <= 500000.0
	trait_component_fusion_gene_tissue_df <- fusion_gene_tissue_df[cis_window_gene_indices,]


	# Now loop through (gene, tissue pairs in cis window of trait component)
	num_gene_tissue_pairs <- dim(trait_component_fusion_gene_tissue_df)[1]

	# Initialize output data
	pred_expr_mat <- matrix(0, 489, num_gene_tissue_pairs)  #489 as we have 489 genotyped individuals
	gene_tissue_pair_names <- c()
	gene_tissue_pair_z_scores <- c()

	if (num_gene_tissue_pairs > 0) {
	# Write table to output
	nominal_pvalue_table_file_name <- paste0(gtex_fusion_multivariate_associations_dir, component_name, "_twas_nominal_association_table.tsv")
	write.table(trait_component_fusion_gene_tissue_df,nominal_pvalue_table_file_name, quote=FALSE, sep="\t",row.names = FALSE)
	twas_nominal_pvalue_table_names <- c(twas_nominal_pvalue_table_names, nominal_pvalue_table_file_name)


	for (gene_tissue_iter in 1:num_gene_tissue_pairs) {

		# Load in info for this gene-tissue pair
		gene_id <- trait_component_fusion_gene_tissue_df$gene[gene_tissue_iter]
		tissue_id <- trait_component_fusion_gene_tissue_df$tissue[gene_tissue_iter]
		wgt_file <- trait_component_fusion_gene_tissue_df$wgt_file[gene_tissue_iter]
		top_model <- trait_component_fusion_gene_tissue_df$top_model[gene_tissue_iter]
		zscore <- trait_component_fusion_gene_tissue_df$twas_nominal_z[gene_tissue_iter]

		#top_model="blup" ###TEMPORARY###
		###########################
		# CP AND PASTE FUSION_assoc_test.R code to reconstruct LD Matrix
		###########################
		
		# Load weights
		load(wgt_file)
		# Remove NAs (these should not be here)
		wgt.matrix[is.na(wgt.matrix)] = 0
	
		# Match up the SNPs and weights
		m = match( snps[,2] , genos$bim[,2] )
		m.keep = !is.na(m)
		snps = snps[m.keep,]
		wgt.matrix = wgt.matrix[m.keep,,drop=F]
		cur.genos = scale(genos$bed[,m[m.keep]])
		cur.bim = genos$bim[m[m.keep],]
		# Flip WEIGHTS for mismatching alleles
		qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
		wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]


		#cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)	

		mod.best = which(colnames(wgt.matrix)==top_model)
		if (length(mod.best) != 1) {
			print("ASSUMPTION EROROR")
		}
		if (top_model == "top1") {
			wgt.matrix[ -which.max(wgt.matrix[,mod.best]^2)  , mod.best] = 0
		}


		# Get predicted expression in reference panel
		best_weights = (wgt.matrix[, mod.best])

		# Quick errror check
		if (sum(rownames(wgt.matrix) != colnames(cur.genos)) > 0) {
			print("FATAL ASSUMPTION EROROROR")
		}

		#saveRDS(cur.genos, "cur.genos.RDS")
		#saveRDS(best_weights, "best_weights.RDS")

		pred_expr = cur.genos %*% best_weights

		# Add to data to keep track
		pred_expr_mat[, gene_tissue_iter] = pred_expr
		gene_tissue_pair_names <- c(gene_tissue_pair_names, paste0(gene_id, ":", tissue_id))
		gene_tissue_pair_z_scores <- c(gene_tissue_pair_z_scores, zscore)

		#####################
		# This remaining info is no longer necessary (just shows we recapitulate the z scores from fusion-twas)


		# Match up the SNPs and the summary stats ###TEMPORARY###
		#m = match(cur.bim[,2] , sumstat$SNP)

		#cur.Z = sumstat$Z[m]
		#cur.miss = is.na(cur.Z)

		#cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
		#cur.impz = cur.wgt %*% cur.Z[!cur.miss]
		#cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
		#cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

		#all.r2pred = rep(1,length(cur.Z))
		#all.r2pred[ cur.miss ] = cur.r2pred

		#cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
		#cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]

		#cur.twas = cur.twasz / sqrt(cur.twasr2pred)

		#gene_tissue_pair_z_scores <- c(gene_tissue_pair_z_scores, cur.twas) ###TEMPORARY###
		#print(trait_component_fusion_gene_tissue_df[gene_tissue_iter,])
		#print(cur.twas)


	}

	cor_pred_expr = cor(pred_expr_mat)
	#blood_WHITE_COUNT_chr21_39095717_C_T_21:37522576:40522576
	
	# Save gene-tissue names
	twas_gene_tissue_name_file <- paste0(gtex_fusion_multivariate_associations_dir, component_name, "_twas_gene_tissue_names.RDS")
	saveRDS(gene_tissue_pair_names, twas_gene_tissue_name_file)
	gene_tissue_names_files <- c(gene_tissue_names_files, twas_gene_tissue_name_file)	

	# Save z scores
	twas_z_score_file <- paste0(gtex_fusion_multivariate_associations_dir, component_name, "_twas_z_scores.RDS")
	saveRDS(gene_tissue_pair_z_scores, twas_z_score_file)
	nominal_z_files <- c(nominal_z_files, twas_z_score_file)	

	# Save correlation of predicted expression
	twas_expr_corr_file <- paste0(gtex_fusion_multivariate_associations_dir, component_name, "_twas_pred_expr_corr.RDS")
	saveRDS(cor_pred_expr, twas_expr_corr_file)
	pred_expr_corr_files <- c(pred_expr_corr_files, twas_expr_corr_file)	


	tryCatch(
	{
	print(gene_tissue_pair_names)
	print(gene_tissue_pair_z_scores)
	fitted_rss <- susie_rss(z=gene_tissue_pair_z_scores, R=cor_pred_expr,n=gwas_sample_size,  L = 20, max_iter = 3000)
	
	# Save Susie RDS
	susie_rds_file <- paste0(gtex_fusion_multivariate_associations_dir, component_name, "_multivariate_twas_susie.RDS")
	saveRDS(fitted_rss, susie_rds_file)
	twas_multivariate_susie_object_names <- c(twas_multivariate_susie_object_names, susie_rds_file)

	print(summary(fitted_rss))

    },
    error = function(e) {
    	twas_multivariate_susie_object_names <- c(twas_multivariate_susie_object_names, "NA")
    }
    )



	} else {
		twas_nominal_pvalue_table_names <- c(twas_nominal_pvalue_table_names, "NA")
		twas_multivariate_susie_object_names <- c(twas_multivariate_susie_object_names, "NA")
		nominal_z_files <- c(nominal_z_files, "NA")
		pred_expr_corr_files <- c(pred_expr_corr_files, "NA")
		gene_tissue_names_files <- c(gene_tissue_names_files, "NA")	

	}
}


organized_susie_trait_component_df$twas_nominal_association_table <- twas_nominal_pvalue_table_names
organized_susie_trait_component_df$twas_multivariate_twas_obj <- twas_multivariate_susie_object_names
organized_susie_trait_component_df$nominal_z_score_obj <- nominal_z_files
organized_susie_trait_component_df$pred_expression_correlation_object <- pred_expr_corr_files
organized_susie_trait_component_df$gene_tissue_names <- gene_tissue_names_files



output_file <- paste0(gtex_fusion_multivariate_associations_dir, trait_name, "_", global_chrom_num, "_component_multivariate_twas_overlaps.txt")

write.table(organized_susie_trait_component_df,output_file, quote=FALSE, sep="\t",row.names = FALSE)






