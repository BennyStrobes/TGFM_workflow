args = commandArgs(trailingOnly=TRUE)
library(ctwas)
library(tools)
options(warn=1)




######################
# Command line args
######################
ctwas_gwas_sumstats_file = args[1]
LD_dir = args[2]
fusion_gene_model_root = args[3]
output_root = args[4]


# Load in z-scores
z_snp <- readRDS(ctwas_gwas_sumstats_file)


res <- impute_expr_z(z_snp = z_snp,
                     weight = fusion_gene_model_root,
                     ld_R_dir = LD_dir,
                     outputdir = output_root,
                     outname = "testing",
                     harmonize_z = F,
                     harmonize_wgt = F,
                     strand_ambig_action_z = "none",
                     recover_strand_ambig_wgt = F)


z_gene <- res$z_gene
z_snp <- res$z_snp
save(z_gene, file = paste0(output_root, "/", "testing", "_z_gene.Rd"))
save(z_snp, file = paste0(output_root, "/", "testing", "_z_snp.Rd"))


#### TEMP
load(paste0(output_root, "/", "testing", "_z_gene.Rd"))
load(paste0(output_root, "/", "testing", "_z_snp.Rd"))
##########


# Load in previously generated imputed genes
ld_exprvarfs <- paste0(output_root, "/", "testing", "_chr", 1:22, ".exprvar")

# Get name of genomic regions file
genomic_regions_file = paste0(LD_dir, "ctwas_quasi_independent_regions.bed")


######
# Run ctwas
thin <- 0.1
max_snp_region <- 20000
ncore <- 1
#We pass these arguments the ctwas_rss and run the full analysis. As specified, using 6 cores and with 56GB of RAM available, this step took approximately 7 hours.


# estimating parameters
ctwas_rss(z_gene = z_gene,
          z_snp = z_snp,
          ld_exprvarfs = ld_exprvarfs,
          ld_R_dir = LD_dir,
          ld_regions_custom = genomic_regions_file,
          outputdir = output_root,
          outname = "testing",
          thin = thin,
          max_snp_region = max_snp_region,
          ncore = ncore)