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

create_gene_tissue_df <- function(tissue_names, pseudotissue_gtex_susie_pmces_fusion_weights_dir, global_chrom_num, trait_name) {
	gene_arr <- c()
	tissue_arr <- c()
	wgt_file_arr <- c()
	chrom_num_arr <- c()
	tss_arr <- c()
	top_model_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]

		# Tissue pos file (containins all genes with wgt files AND their TSS)
		tissue_pos_file <- paste0(pseudotissue_gtex_susie_pmces_fusion_weights_dir, tissue_name, ".pos")
		tissue_pos_df <- read.table(tissue_pos_file, header=TRUE)  # Load in file
		tissue_pos_df <- tissue_pos_df[tissue_pos_df$CHR==global_chrom_num,]  # Filter to genes on this chromosome

		# Loop through genes
		ngenes <- dim(tissue_pos_df)[1]
		for (gene_num in 1:ngenes) {
			gene_arr <- c(gene_arr, tissue_pos_df$ID[gene_num])
			tissue_arr <- c(tissue_arr, tissue_name)
			wgt_file_arr <- c(wgt_file_arr, tissue_pos_df$WGT[gene_num])
			chrom_num_arr <- c(chrom_num_arr, global_chrom_num)
			tss_arr <- c(tss_arr, tissue_pos_df$P0[gene_num])
			top_model_arr <- c(top_model_arr, "susie_pmces")
		}
	}
	df <- data.frame(gene=gene_arr, tissue=tissue_arr, wgt_file=wgt_file_arr, chrom_num=chrom_num_arr, gene_tss=tss_arr, top_model=top_model_arr)
	return(df)
}

randomly_sample_susie_weights <- function(susie_mu, susie_alpha, susie_sdev) {
	num_components <- dim(susie_mu)[1]
	num_variants <- dim(susie_mu)[2]
	# Initialize sampled weights 
	weights <- rep(0, num_variants)

	for (component_num in 1:num_components) {
		# First sample variant from categorical distribution
		component_variant = sample(1:num_variants, 1, replace=TRUE, prob=susie_alpha[component_num,])

		variant_beta_expected_value = susie_mu[component_num, component_variant]
		variant_beta_sdev = susie_sdev[component_num, component_variant]
		component_beta = rnorm(n=1, mean=variant_beta_expected_value, sd=variant_beta_sdev)
		weights[component_variant] = weights[component_variant] + component_beta
	}
	return(weights)
}


global_chrom_num=args[1]
trait_name=args[2]
ukbb_genome_wide_susie_organized_results_dir=args[3]
gtex_pseudotissue_file=args[4]
pseudotissue_gtex_susie_pmces_fusion_weights_dir=args[5]
gtex_fusion_processed_intermediate_data=args[6]
gwas_sample_size=as.numeric(args[7])
pseudotissue_gtex_stochastic_multivariate_twas_dir=args[8] # OUTPUTDIR

num_sim=100

# Extract tissue names
tissue_names <- read.table(gtex_pseudotissue_file, header=TRUE, sep="\t")$pseudotissue_name


# Create data frame of (tissue, gene) on this chromosome that we have wgt files for
fusion_gene_tissue_df <- create_gene_tissue_df(tissue_names, pseudotissue_gtex_susie_pmces_fusion_weights_dir, global_chrom_num, trait_name)


# Extract susie genome wide trait components on this chromosome
organized_susie_trait_component_file <- paste0(ukbb_genome_wide_susie_organized_results_dir, trait_name, "_organized_susie_components.txt")  # File name
organized_susie_trait_component_df <- read.table(organized_susie_trait_component_file, header=TRUE)  # load in file
organized_susie_trait_component_df <- organized_susie_trait_component_df[organized_susie_trait_component_df$chrom_num==global_chrom_num,]  # Filter to chromosomes on the chromosome of interest

# Load in sumstat file (necessary for variant orientation)
sumstat_file <- paste0(gtex_fusion_processed_intermediate_data, trait_name, "_fusion_processed_sumstats_hg38.txt")
sumstat = read.table(sumstat_file,head=T,as.is=T)

# Load in Ref LD data
ref_ld_chr = paste0(gtex_fusion_processed_intermediate_data, "1000G.EUR.gtex_formatted.hg38.", global_chrom_num)
genos = read_plink(ref_ld_chr,impute="avg")
#saveRDS(genos, file = "temp_genos.RDS")




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





trait_component_genes_table_names <- c()
gene_tissue_names_files <- c()
nominal_z_files <- c()
stochastic_nominal_z_files <- c()
pip_files <- c()
stochastic_pip_files <- c()


# Now begin our loop through Susie genome wide trait components
num_susie_trait_components <- dim(organized_susie_trait_component_df)[1]
print(num_susie_trait_components)
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

	# Initialize output data (for expected value analysis)
	gene_tissue_pair_names <- c()
	pred_expr_mat <- matrix(0, 489, num_gene_tissue_pairs)  #489 as we have 489 genotyped individuals
	gene_tissue_pair_z_scores <- c()
	
	# Initialize output data (for stochastic analysis)
	stochastic_pred_expr_mat_list <- list()
	stochastic_z_scores_list <- list()
	stochastic_pip_list <- list()
	for (sim_num in 1:num_sim) {
		stochastic_pred_expr_mat_list[[sim_num]] = matrix(0, 489, num_gene_tissue_pairs)
		stochastic_z_scores_list[[sim_num]] = rep(0, num_gene_tissue_pairs)
		stochastic_pip_list[[sim_num]] = rep(NA, num_gene_tissue_pairs)
	}


	if (num_gene_tissue_pairs > 0) {
	# Write table to output
	trait_component_genes_table_file_name <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_component_genes.tsv")
	write.table(trait_component_fusion_gene_tissue_df, trait_component_genes_table_file_name, quote=FALSE, sep="\t",row.names = FALSE)
	trait_component_genes_table_names <- c(trait_component_genes_table_names, trait_component_genes_table_file_name)


	for (gene_tissue_iter in 1:num_gene_tissue_pairs) {
		print(gene_tissue_iter)

		# Load in info for this gene-tissue pair
		gene_id <- trait_component_fusion_gene_tissue_df$gene[gene_tissue_iter]
		tissue_id <- trait_component_fusion_gene_tissue_df$tissue[gene_tissue_iter]
		wgt_file <- trait_component_fusion_gene_tissue_df$wgt_file[gene_tissue_iter]
		top_model <- trait_component_fusion_gene_tissue_df$top_model[gene_tissue_iter]

		###########################
		# CP AND PASTE FUSION_assoc_test.R code to reconstruct LD Matrix
		###########################
		
		# Load weights
		load(wgt_file)
		# Remove NAs (these should not be here)
		wgt.matrix[is.na(wgt.matrix)] = 0
		susie_sdev = sqrt(susie_mu2 - (susie_mu^2))
		# Also have susie_mu and susie_alpha
	
		# Match up the SNPs and weights
		m = match( snps[,2] , genos$bim[,2] )
		m.keep = !is.na(m)
		snps = snps[m.keep,]
		wgt.matrix = wgt.matrix[m.keep,,drop=F]
		susie_mu = susie_mu[,m.keep]
		susie_alpha = susie_alpha[,m.keep]
		susie_sdev = susie_sdev[,m.keep]
		cur.genos = scale(genos$bed[,m[m.keep]])
		cur.bim = genos$bim[m[m.keep],]
		# Flip WEIGHTS for mismatching alleles
		qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
		wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
		susie_mu[,qc$flip] = -1*susie_mu[,qc$flip]



		cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)	

		mod.best = which(colnames(wgt.matrix)==top_model)
		if (length(mod.best) != 1) {
			print("ASSUMPTION EROROR")
		}


		# Get predicted expression in reference panel
		best_weights = (wgt.matrix[, mod.best])
		gene_tissue_pair_names <- c(gene_tissue_pair_names, paste0(gene_id, ":", tissue_id))
		pred_expr = cur.genos %*% best_weights
		# Add to data to keep track
		pred_expr_mat[, gene_tissue_iter] = pred_expr

		# Match up the SNPs and the summary stats (Only needs to be done once for all stochastic iterations)
		m = match(cur.bim[,2] , sumstat$SNP)

		cur.Z = sumstat$Z[m]
		cur.miss = is.na(cur.Z)

		cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
		cur.impz = cur.wgt %*% cur.Z[!cur.miss]
		cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
		cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
		#all.r2pred = rep(1,length(cur.Z))
		#all.r2pred[ cur.miss ] = cur.r2pred

		# Needs to be done for each stochastic iterations
		cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
		cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]
		twas_z_score = cur.twasz / sqrt(cur.twasr2pred)

		# Save data
		gene_tissue_pair_z_scores <- c(gene_tissue_pair_z_scores, twas_z_score)

		# Draw random samples from effect size distributions
		for (sim_num in 1:num_sim) {
			simulated_weights <- randomly_sample_susie_weights(susie_mu, susie_alpha, susie_sdev)

			# Compute twas z-score for stochastically simulated weights
			cur.twasz = simulated_weights %*% cur.Z
			cur.twasr2pred = simulated_weights %*% cur.LD %*% simulated_weights
			twas_z_score = cur.twasz / sqrt(cur.twasr2pred)

			stochastic_z_scores_list[[sim_num]][gene_tissue_iter] = twas_z_score
			stochastic_pred_expr_mat_list[[sim_num]][, gene_tissue_iter] = cur.genos %*% simulated_weights
		}
	}


	# Save gene-tissue names
	twas_gene_tissue_name_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_twas_gene_tissue_names.RDS")
	saveRDS(gene_tissue_pair_names, twas_gene_tissue_name_file)
	gene_tissue_names_files <- c(gene_tissue_names_files, twas_gene_tissue_name_file)	
	
	# Save z scores
	twas_z_score_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_twas_z_scores.RDS")
	saveRDS(gene_tissue_pair_z_scores, twas_z_score_file)
	nominal_z_files <- c(nominal_z_files, twas_z_score_file)

	# Save stochastic z scores
	stochastic_z_score_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_stochastic_twas_z_scores.RDS")
	saveRDS(stochastic_z_scores_list, stochastic_z_score_file)
	stochastic_nominal_z_files <- c(stochastic_nominal_z_files, stochastic_z_score_file)	

	# Save stochastic z scores
	stochastic_pred_expr_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_stochastic_pred_expr.RDS")
	saveRDS(stochastic_pred_expr_mat_list, stochastic_pred_expr_file)


	tryCatch(
	{
		fitted_rss <- susie_rss(z=gene_tissue_pair_z_scores, R=cor(pred_expr_mat),n=gwas_sample_size, L = 20, max_iter = 3000)
		# Save PIPS
		pip_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_twas_pips.RDS")
		saveRDS(fitted_rss$pip, pip_file)
		pip_files <- c(pip_files, pip_file)



    },
    error = function(e) {
    	print("SuSiE did not converge")
    	pip_files <- c(pip_files, "NA")
    }
    )

    for (sim_num in 1:num_sim) {
    	tryCatch(
		{
    		fitted_rss_sim <- susie_rss(z=stochastic_z_scores_list[[sim_num]], R=cor(stochastic_pred_expr_mat_list[[sim_num]]),n=gwas_sample_size, L = 20, max_iter = 3000)
    		
    		stochastic_pip_list[[sim_num]] = fitted_rss_sim$pip
    	},
    	error = function(e) {
    		print("SuSiE did not converge")
    	}
    	)
    }

	# Save stochastic pips
	stochastic_pip_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, component_name, "_stochastic_twas_pips.RDS")
	saveRDS(stochastic_pip_list, stochastic_pip_file)
	stochastic_pip_files <- c(stochastic_pip_files, stochastic_pip_file)	


	} else {
		trait_component_genes_table_names <- c(trait_component_genes_table_names, "NA")
		gene_tissue_names_files <- c(gene_tissue_names_files, "NA")
		nominal_z_files <- c(nominal_z_files, "NA")
		stochastic_nominal_z_files <- c(stochastic_nominal_z_files, "NA")
		pip_files <- c(pip_files, "NA")
		stochastic_pip_files <- c(stochastic_pip_files, "NA")	
	}
}


organized_susie_trait_component_df$trait_component_genes_table <- trait_component_genes_table_names
organized_susie_trait_component_df$gene_tissue_names_file <- gene_tissue_names_files
organized_susie_trait_component_df$nominal_z_score_obj <- nominal_z_files
organized_susie_trait_component_df$stochastic_nominal_z_score_obj <- stochastic_nominal_z_files
organized_susie_trait_component_df$pip_obj <- pip_files
organized_susie_trait_component_df$stochastic_pip_obj <- stochastic_pip_files



output_file <- paste0(pseudotissue_gtex_stochastic_multivariate_twas_dir, trait_name, "_", global_chrom_num, "_component_stochastic_multivariate_twas_overlaps.txt")

write.table(organized_susie_trait_component_df,output_file, quote=FALSE, sep="\t",row.names = FALSE)

