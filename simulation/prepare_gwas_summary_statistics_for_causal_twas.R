args = commandArgs(trailingOnly=TRUE)
options(warn=1)









input_gwas_file = args[1]
n_gwas_individuals = as.numeric(args[2])
ctwas_gwas_sumstats_file = args[3]
ctwas_gwas_ss_file = args[4]


# Load in gwas data
gwas_input <- read.table(input_gwas_file,header=TRUE)
# extract vector of z-scores
z_vec = gwas_input$BETA/sqrt(gwas_input$BETA_VAR)

# Rearange into new data-frame
z_df <- data.frame(id=as.character(gwas_input$SNP), A1=gwas_input$A1, A2=gwas_input$A2, z=z_vec)



# save the formatted z-scores and GWAS sample size
saveRDS(z_df, file=ctwas_gwas_sumstats_file)
saveRDS(n_gwas_individuals, file=ctwas_gwas_ss_file)



#z_snp <- readRDS(ctwas_gwas_sumstats_file)
#gwas_n <- readRDS(ctwas_gwas_ss_file)

