args = commandArgs(trailingOnly=TRUE)
library("optparse")
library('plink2R')
library('susieR')
library('methods')
library(hash)
options(warn=1)

# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
	if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
	else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
	return( eff.wgt )
}


tissue_specific_analysis <- function(composit_tissue, covariate_file, OUT, gcta_path) {
	tissue_specific_stem <- paste0(OUT, "_", composit_tissue)

	# Load in fam file 
	fam_file <- paste0(tissue_specific_stem, ".fam")
	fam = read.table(fam_file,as.is=T)

	# Extract Pheno file from fam file 
	pheno = fam[,c(1,2,6)]
	raw.pheno.file = paste0(tissue_specific_stem, ".pheno")
	write.table(pheno,quote=F,row.names=F,col.names=F,file=raw.pheno.file)

	# Load in covariates
	covar = ( read.table(covariate_file,as.is=T,head=F) )
	
	# Match up sample data between pheno and covariates
	m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,]

	# Regress out covariates on phenotype
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
	pheno[,3] = scale(reg$resid)

	# Save pheno to output
	pheno.file = paste0(tissue_specific_stem, ".pheno.resid")
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)

	# recode to the intersection of samples and new phenotype (and use residual expression)
	temp_tissue_specific_stem <- paste0(tissue_specific_stem, "_tmper")
	geno.file <- temp_tissue_specific_stem
	arg = paste( "plink --allow-no-sex --bfile ",tissue_specific_stem," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --keep-allele-order --out ",geno.file,sep='')
	system(arg, ignore.stdout=T)


	#########################
	# Heritability analysis
	########################
	# 1. generate GRM
	arg = paste( "plink --allow-no-sex --bfile ",temp_tissue_specific_stem," --make-grm-bin --out ",temp_tissue_specific_stem,sep='')
	system(arg, ignore.stdout=T)
	# 2. estimate heritability
	arg = paste( gcta_path ," --grm ",temp_tissue_specific_stem," --pheno ",raw.pheno.file," --qcovar ",covariate_file," --out ",temp_tissue_specific_stem," --reml --reml-no-constrain --reml-lrt 1",sep='')
	system(arg, ignore.stdout=T)
	# 3. evaluate LRT and V(G)/Vp
	if ( !file.exists( paste(temp_tissue_specific_stem,".hsq",sep='') ) ) {
		cat("GCTA heritability analysis output file not exist, likely GCTA could not converge, skipping gene\n",file=stderr())
		hsq = NA
		hsq.pv = NA
	} else {
		hsq.file = read.table(file=paste(temp_tissue_specific_stem,".hsq",sep=''),as.is=T,fill=T)
		hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
		hsq.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2]))
	}

	# read in genotypes
	genos = read_plink(geno.file,impute="avg")

	# Scale genotypes
	mafs = apply(genos$bed,2,mean)/2
	sds = apply(genos$bed,2,sd)
	# important : genotypes are standardized and scaled here:
	genos$bed_raw = genos$bed
	genos$bed = scale(genos$bed)
	pheno = genos$fam[,c(1,2,6)]
	pheno[,3] = scale(pheno[,3])




	# Get marginal effects of each snp on on pheno
	sample_size = dim(pheno)[1]
	nvar = dim(genos$bed)[2]
	betas = rep(0.0, nvar)
	beta_ses = rep(0.0, nvar)
	for (var_index in 1:nvar) {
		fit <- lm(pheno[,3] ~ genos$bed[,var_index] +0)
		beta = summary(fit)$coefficients[1]
		beta_se = summary(fit)$coefficients[2]
		betas[var_index] = beta
		beta_ses[var_index] = beta_se
	}

	my_list <- list("beta" = betas, "beta_se" = beta_ses, "genotype" = genos$bed_raw, "heritability"=hsq, "heritability_p"=hsq.pv, "bim"=genos$bim)
	return(my_list)
}


fixed_meta_analysis <- function(beta_vec, beta_se_vec) {
	variances = beta_se_vec**2
	wt = 1.0/variances


	summ = ((wt%*%beta_vec)[1,1])/sum(wt)
	varsum = sum(wt*wt*variances)/(sum(wt)**2)
	se = sqrt(varsum)
	
	my_list <- list("meta_beta"=summ, "meta_beta_se"=se)
	return(my_list)
}

pseudotissue_name <- args[1]
composit_tissue_string <- args[2]
gene_id <- args[3]
chrom_num <- args[4]
composit_covariate_file <- args[5]
OUT <- args[6]
FINAL_OUT <- args[7]
gcta_path <- args[8]
gemma_path <- args[9]



# Get array of composit tissues from composit tissue string
composit_tissues <- strsplit(composit_tissue_string, ",")[[1]]
num_composit_tissues = length(composit_tissues)

# Get array of covariate files from composit covariate file string
covariate_files <- strsplit(composit_covariate_file, ",")[[1]]


if (num_composit_tissues == 1) {
	# Extract info from tissue 1
	composit_tissue1 <- composit_tissues[1]
	covariate_file1 <- covariate_files[1]
	ts_obj1 <- tissue_specific_analysis(composit_tissue1, covariate_file1, OUT, gcta_path)

	# Variant Bim
	variant_bim = ts_obj1$bim

	# Combine with betas and standard errors with meta-analysis
	num_var = length(ts_obj1$beta)
	meta_beta <- ts_obj1$beta
	meta_beta_se <- ts_obj1$beta_se

	# Compute concatenated LD
	meta_ld = cor(ts_obj1$genotype)

	# Meta-analyzed sample size
	meta_sample_size = dim(ts_obj1$genotype)[1]

	# Names of variants
	variant_names = colnames(ts_obj1$genotype)

	# Meta heritabilities
	combined_heritabilities = c(ts_obj1$heritability[1])
	combined_heritabilities_p = c(ts_obj1$heritability_p[1])

}



if (num_composit_tissues == 2) {
	# Extract info from tissue 1
	composit_tissue1 <- composit_tissues[1]
	covariate_file1 <- covariate_files[1]
	ts_obj1 <- tissue_specific_analysis(composit_tissue1, covariate_file1, OUT, gcta_path)

	# Extract info from tissue 2
	composit_tissue2 <- composit_tissues[2]
	covariate_file2 <- covariate_files[2]
	ts_obj2 <- tissue_specific_analysis(composit_tissue2, covariate_file2, OUT, gcta_path)

	# Variant Bim
	variant_bim = ts_obj1$bim

	# Combine with betas and standard errors with meta-analysis
	num_var = length(ts_obj1$beta)
	meta_beta <- rep(0, num_var)
	meta_beta_se <- rep(0, num_var)

	for (var_num in 1:num_var) {
		beta_vec <- c(ts_obj1$beta[var_num], ts_obj2$beta[var_num])
		beta_se_vec <- c(ts_obj1$beta_se[var_num], ts_obj2$beta_se[var_num])
		meta_analysis_obj = fixed_meta_analysis(beta_vec, beta_se_vec)

		meta_beta[var_num] = meta_analysis_obj$meta_beta
		meta_beta_se[var_num] = meta_analysis_obj$meta_beta_se
	}

	# Compute concatenated LD
	combined_geno = rbind(ts_obj1$genotype, ts_obj2$genotype)
	meta_ld = cor(combined_geno)

	# Meta-analyzed sample size
	meta_sample_size = dim(combined_geno)[1]

	# Names of variants
	variant_names = colnames(combined_geno)

	# Meta heritabilities
	combined_heritabilities = c(ts_obj1$heritability[1], ts_obj2$heritability[1])
	combined_heritabilities_p = c(ts_obj1$heritability_p[1], ts_obj2$heritability_p[1])

}


if (num_composit_tissues == 3) {
	# Extract info from tissue 1
	composit_tissue1 <- composit_tissues[1]
	covariate_file1 <- covariate_files[1]
	ts_obj1 <- tissue_specific_analysis(composit_tissue1, covariate_file1, OUT, gcta_path)

	# Extract info from tissue 2
	composit_tissue2 <- composit_tissues[2]
	covariate_file2 <- covariate_files[2]
	ts_obj2 <- tissue_specific_analysis(composit_tissue2, covariate_file2, OUT, gcta_path)

	# Extract info from tissue 3
	composit_tissue3 <- composit_tissues[3]
	covariate_file3 <- covariate_files[3]
	ts_obj3 <- tissue_specific_analysis(composit_tissue3, covariate_file3, OUT, gcta_path)

	# Variant Bim
	variant_bim = ts_obj1$bim

	# Combine with betas and standard errors with meta-analysis
	num_var = length(ts_obj1$beta)
	meta_beta <- rep(0, num_var)
	meta_beta_se <- rep(0, num_var)

	for (var_num in 1:num_var) {
		beta_vec <- c(ts_obj1$beta[var_num], ts_obj2$beta[var_num], ts_obj3$beta[var_num])
		beta_se_vec <- c(ts_obj1$beta_se[var_num], ts_obj2$beta_se[var_num], ts_obj3$beta_se[var_num])
		meta_analysis_obj = fixed_meta_analysis(beta_vec, beta_se_vec)

		meta_beta[var_num] = meta_analysis_obj$meta_beta
		meta_beta_se[var_num] = meta_analysis_obj$meta_beta_se
	}

	# Compute concatenated LD
	combined_geno = rbind(ts_obj1$genotype, ts_obj2$genotype, ts_obj3$genotype)
	meta_ld = cor(combined_geno)

	# Meta-analyzed sample size
	meta_sample_size = dim(combined_geno)[1]

	# Names of variants
	variant_names = colnames(combined_geno)

	# Meta heritabilities
	combined_heritabilities = c(ts_obj1$heritability[1], ts_obj2$heritability[1], ts_obj3$heritability[1])
	combined_heritabilities_p = c(ts_obj1$heritability_p[1], ts_obj2$heritability_p[1],ts_obj3$heritability_p[1])

}






# Run SuSiE
susie_ss = susie_rss(bhat=meta_beta, shat=meta_beta_se,R=meta_ld, n=meta_sample_size)
print(summary(susie_ss))
# Extract relevent fields
susie_mu=data.frame(susie_ss$mu)
susie_mu2=data.frame(susie_ss$mu2)
susie_alpha=data.frame(susie_ss$alpha)
susie_V =susie_ss$V
susie_converged=susie_ss$converged
susie_cs=susie_ss$sets$cs_index
if (is.null(susie_cs)) {
	component_bool=FALSE
} else {
	component_bool=TRUE
}

print(component_bool)


# Save results to output
save(variant_names, variant_bim, combined_heritabilities, combined_heritabilities_p, meta_sample_size, susie_mu, susie_mu2,susie_alpha, susie_V, susie_converged, component_bool, susie_cs, file = paste( FINAL_OUT , ".wgt.RDat" , sep='' ) )


