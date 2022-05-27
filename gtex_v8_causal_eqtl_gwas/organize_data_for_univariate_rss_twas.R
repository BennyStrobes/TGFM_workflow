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
	#snp[["keep"]] = !((a1=="ZZZZZZZ")) # MY Hack to ensure no snps are thrown out b/c of commented out strand ambiguity below b/c not proficient in R
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
	snp[["flip"]] = (a1 == ref2 & a2 == ref1)
	return(snp)
}

create_gene_tissue_df <- function(tissue_names, pseudotissue_gtex_susie_pmces_fusion_weights_dir, global_chrom_num, trait_name, gene_version) {
	gene_arr <- c()
	tissue_arr <- c()
	wgt_file_arr <- c()
	chrom_num_arr <- c()
	tss_arr <- c()
	top_model_arr <- c()

	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]

		# Tissue pos file (containins all genes with wgt files AND their TSS)
		tissue_pos_file <- paste0(pseudotissue_gtex_susie_pmces_fusion_weights_dir, tissue_name, "_", gene_version, ".pos")
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
pseudotissue_gtex_rss_multivariate_twas_data_dir=args[8] # OUTPUTDIR
gene_version=args[9]




# Extract tissue names
tissue_names <- read.table(gtex_pseudotissue_file, header=TRUE, sep="\t")$pseudotissue_name


# Create data frame of (tissue, gene) on this chromosome that we have wgt files for
fusion_gene_tissue_df <- create_gene_tissue_df(tissue_names, pseudotissue_gtex_susie_pmces_fusion_weights_dir, global_chrom_num, trait_name, gene_version)


# Extract susie genome wide trait components on this chromosome
organized_susie_trait_component_file <- paste0(ukbb_genome_wide_susie_organized_results_dir, trait_name, "_organized_susie_filtered_components.txt")  # File name
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
sumstat$beta[ qc$flip ] = -1 * sumstat$beta[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
	genos$bim = genos$bim[qc$keep,]
	genos$bed = genos$bed[,qc$keep]
	sumstat = sumstat[qc$keep,]
}


component_names <- c()
num_genes <- c()
trait_component_genes_table_names <- c()


# Now begin our loop through Susie genome wide trait components
num_susie_trait_components <- dim(organized_susie_trait_component_df)[1]
print(num_susie_trait_components)
for (trait_component_iter in 1:num_susie_trait_components) {
	print(trait_component_iter)
	# Extract relvent fields for this trait component
	trait_component_lead_variant <- organized_susie_trait_component_df$lead_variant_name[trait_component_iter]  # Name of lead variant for trait component
	trait_component_window_id <- organized_susie_trait_component_df$window_id[trait_component_iter]  # Name of window id defining trait component
	trait_component_lead_variant_pos <- organized_susie_trait_component_df$lead_variant_position[trait_component_iter]  # Position of lead variant for trait component
	component_name <- paste0(trait_name, "_", trait_component_lead_variant, "_", trait_component_window_id)

	# Extract gene-tissue pairs in cis window of this trait component
	cis_window_gene_indices <- abs(fusion_gene_tissue_df$gene_tss - trait_component_lead_variant_pos) <= 250000.0
	trait_component_fusion_gene_tissue_df <- fusion_gene_tissue_df[cis_window_gene_indices,]


	# Now loop through (gene, tissue pairs in cis window of trait component)
	num_gene_tissue_pairs <- dim(trait_component_fusion_gene_tissue_df)[1]

	# Initialize output data (for expected value analysis)
	gene_tissue_pair_names <- c()
	geno_file_names <- c()
	bim_file_names <- c()
	eqtl_susie_mu_file_names <- c()
	eqtl_susie_alpha_file_names <- c()
	eqtl_susie_mu_sd_file_names <- c()	
	gwas_beta_file_names <- c()
	gwas_beta_se_file_names <- c()

	if (num_gene_tissue_pairs > 0) {
	# Write table to output


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


		mod.best = which(colnames(wgt.matrix)==top_model)
		if (length(mod.best) != 1) {
			print("ASSUMPTION EROROR")
		}


		# Get predicted expression in reference panel
		best_weights = (wgt.matrix[, mod.best])

		gene_tissue_pair_names <- c(gene_tissue_pair_names, paste0(gene_id, ":", tissue_id))


		# Match up the SNPs and the summary stats (Only needs to be done once for all stochastic iterations)
		m = match(cur.bim[,2] , sumstat$SNP)

		cur.Z = sumstat$Z[m]
		cur.beta = sumstat$beta[m]
		cur.beta_se = sumstat$beta_se[m]



		##########################
		# Start saving data
		##########################
		# Save geno file
		geno_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version, "_", gene_id, "_", tissue_id, "_reference_genotype.tsv")
		write.table(cur.genos, geno_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		geno_file_names <- c(geno_file_names, geno_file_name)

		# Save BIM file
		bim_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version, "_", gene_id, "_", tissue_id, "_reference_genotype_bim.tsv")
		write.table(cur.bim, bim_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		bim_file_names <- c(bim_file_names, bim_file_name)

		# Save eQTL SuSiE mu file
		eqtl_susie_mu_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version,"_", gene_id, "_", tissue_id, "_eqtl_susie_mu.tsv")
		write.table(susie_mu, eqtl_susie_mu_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		eqtl_susie_mu_file_names <- c(eqtl_susie_mu_file_names, eqtl_susie_mu_file_name)

		# Save eQTL SuSiE alpha file
		eqtl_susie_alpha_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version,"_", gene_id, "_", tissue_id, "_eqtl_susie_alpha.tsv")
		write.table(susie_alpha, eqtl_susie_alpha_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		eqtl_susie_alpha_file_names <- c(eqtl_susie_alpha_file_names, eqtl_susie_alpha_file_name)

		# Save eQTL SuSiE mu sd file
		eqtl_susie_mu_sd_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version,"_", gene_id, "_", tissue_id, "_eqtl_susie_mu_sd.tsv")
		write.table(susie_sdev, eqtl_susie_mu_sd_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		eqtl_susie_mu_sd_file_names <- c(eqtl_susie_mu_sd_file_names, eqtl_susie_mu_sd_file_name)

		# Save beta file
		gwas_beta_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version,"_", gene_id, "_", tissue_id, "_gwas_betas.tsv")
		write.table(cur.beta, gwas_beta_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		gwas_beta_file_names <- c(gwas_beta_file_names, gwas_beta_file_name)

		# Save beta se file
		gwas_beta_se_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name,"_", gene_version,"_", gene_id, "_", tissue_id, "_gwas_betas_se.tsv")
		write.table(cur.beta_se, gwas_beta_se_file_name, quote=FALSE, sep="\t",row.names = FALSE, col.names=FALSE)
		gwas_beta_se_file_names <- c(gwas_beta_se_file_names, gwas_beta_se_file_name)

	}

	trait_component_fusion_gene_tissue_df$geno_file_name = as.character(geno_file_names)
	trait_component_fusion_gene_tissue_df$bim_file_name = as.character(bim_file_names)
	trait_component_fusion_gene_tissue_df$eqtl_susie_mu_file_name = as.character(eqtl_susie_mu_file_names)
	trait_component_fusion_gene_tissue_df$eqtl_susie_alpha_file_name = as.character(eqtl_susie_alpha_file_names)
	trait_component_fusion_gene_tissue_df$eqtl_susie_mu_sd_file_name = as.character(eqtl_susie_mu_sd_file_names)
	trait_component_fusion_gene_tissue_df$gwas_beta_file_name = as.character(gwas_beta_file_names)
	trait_component_fusion_gene_tissue_df$gwas_beta_se_file_name = as.character(gwas_beta_se_file_names)

	trait_component_genes_table_file_name <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, component_name, "_", gene_version, "_component_gene_info.tsv")
	write.table(trait_component_fusion_gene_tissue_df, trait_component_genes_table_file_name, quote=FALSE, sep="\t",row.names = FALSE)
	

	component_names <- c(component_names, component_name)
	num_genes <- c(num_genes, num_gene_tissue_pairs)
	trait_component_genes_table_names <- c(trait_component_genes_table_names, trait_component_genes_table_file_name)


	} else {
		trait_component_genes_table_names <- c(trait_component_genes_table_names, "NA")
		component_names <- c(component_names, component_name)
		num_genes <- c(num_genes, num_gene_tissue_pairs)
	}
}


organized_susie_trait_component_df$trait_component_genes_table <- trait_component_genes_table_names
organized_susie_trait_component_df$num_cis_genes <- num_genes


output_file <- paste0(pseudotissue_gtex_rss_multivariate_twas_data_dir, trait_name, "_", gene_version, "_", global_chrom_num, "_component_rss_multivariate_twas_data.txt")
write.table(organized_susie_trait_component_df,output_file, quote=FALSE, sep="\t",row.names = FALSE)




