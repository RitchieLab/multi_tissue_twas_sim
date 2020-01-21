## analyze_twas_simulation.R
## Binglan Li 10/27/2018
## simulate necessary data and analyze TWAS power and type I error


############################################
## Menu
## 1. Food Preparation
## 1.1. Load Necessary Libraries and Scripts
## 1.2. Define Parameters
## 2. Appetizers
## 2.1. Data Simulation
## 2.2. eQTL Detection
## 3. Entree/Main Course
## 3.1. Run single-tissue TWAS
## 3.2. Run integrative TWAS
## 3.3. Evaluate Power and Type I Error Rate of TWAS Results
## 4. Dessert
## 4.1. None. Sorry this is a healthy (aka anti-sweet) restaurant.
############################################


############################################
## 1. Food Preparation
############################################

## 1.1. Load Necessary Libraries and Scripts
### load libraries and functions
suppressMessages(library(optparse))
suppressMessages(library(mvtnorm))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(pls))
suppressMessages(library(GBJ))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
if(!is.loaded("wrapper")){
  #system("R CMD SHLIB optim.c")
  dyn.load("optim.so")
}
source("twas_simulation_util.R")

## 1.2. Define Parameters
# read external 
opt_list <- list(
  make_option("--training", type="integer", default=500, help="eQTL discovery dataset sample size (default = %default)"),
  make_option("--testing", type="integer", default=1000, help="TWAS sample size (default = %default)"),
  make_option("--snps", type="integer", default=60, help="Number of SNPs in a given gene (default = %default)"),
  make_option("--eqtls", type="integer", default=30, help="Number of eQTLs among SNPs in a given gene (default = %default)"),
  make_option("--maf", type="character", default="0.01,0.5", help="Range of minor allele frequency (default = %default)"),
  make_option("--pct-mt-eqtls", type="numeric", default=0.8, help="Percentage of multi-tissue eQTLs among total eQTLs (default = %default)"),
  make_option("--genes", type="integer", default=100, help="Number of genes in a round of simulation (default = %default)"),
  make_option("--expr-tis", type="integer", default=1, help="Number of gene expressing tissues (default = %default)"),
  make_option("--total-tis", type="integer", default=10, help="Number of total tissues (default = %default)"),
  make_option("--h2-ge", type="numeric", default=0.3, help="Heritability of gene expression levels, (default = %default)"),
  make_option("--cor-tissues", type="character", default="0", help="Correlation of gene expression levels among tissues, like '0' or '-1,1' (default = %default)"),
  make_option("--r2-et", type="numeric", default=0.01, help="Variance of simulated traits explained by gene expression levels, (default = %default)"),
  make_option("--output-dir", type="character", action="store", default=paste0(getwd(), "/simulated_data/"), help="Output directory of simulation results"),
  make_option("--output-prefix", type="character", action="store", default="twas_sim", help="Prefix of each output file"),
  make_option("--core", type="numeric", default=1, help="Number of cores to run parallel tasks (default = %default)"),
  make_option("--random-seed", type="numeric", default=1, help="Random seed number (default = %default)")
)
opts <- parse_args(OptionParser(option_list=opt_list))
# parse values
n_train <- opts$training # eQTL discovery dataset sample size
n_test <- opts$testing # TWAS dataset sample size
n_snps <- opts$snps
n_eqtls <- opts$eqtls
maf_range <- as.numeric(strsplit(opts$maf, ",")[[1]])
pct_mt_eqtls <- opts$`pct-mt-eqtls` #pct_mt_eqtls_list <- c(0, 0.2, 0.4, 0.6, 0.8)
n_genes <- opts$genes
n_tissues <- opts$`expr-tis` #n_tissue_list <- c(1, 2, 5, 10)
total_tissues <- opts$`total-tis`
h2_ge <- opts$`h2-ge`
cor_tissues <- as.numeric(strsplit(opts$`cor-tissues`, ",")[[1]])
r2_et <- opts$`r2-et` #r2_et_list <- c(0.00001, 0.0005, 0.005, 0.01)
# set random seed
random_seed <- n_tissues*100 + opts$`random-seed`
set.seed(random_seed, kind = "L'Ecuyer-CMRG")

output_dir <- opts$`output-dir`
output_prefix <- opts$`output-prefix`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

numCores <- opts$core
print(paste0("numCores=", numCores))

# checkpoints
# no checks at the moment. High-quality food is ready to turn into delicious dishes.

# attain parameters from input variables
# adjust the number of multi-tissue eQTLs according to the number of gene expressing tissues
if(n_tissues == 1){ pct_mt_eqtls <- 0 } # if single tissue, there is no multi-tissue eqtls
if(length(output_prefix) > 0){output_prefix <- paste0(output_prefix, "_")} # reformat prefix if given
if(length(maf_range) == 1){maf_range <- c(maf_range, maf_range)}
n_mt_eqtls <- as.integer(round(n_eqtls * pct_mt_eqtls, digits = 0))
train_ind <- seq(1, n_train) # index for training individuals
test_ind <- seq(n_train + 1, n_train + n_test) # index for testing individuals
tissue_names <- paste0("Tissue", seq(1, n_tissues)) # tissue names
causal_tis <- tissue_names[1] # assuming only one tissue is causal

# print parameters
print(paste0("#tissues=",n_tissues,", pct_mt_eqtls=", pct_mt_eqtls, ", h2_gene_trait=", r2_et, ", cor_tissues=", cor_tissues, ", iteration=", output_prefix))

############################################
## 2. Appetizers
############################################

## 2.1. Data Simulation

### 2.1.1. define paths and parameters
# i.e. path = sample_size/tissue/snp_composition/heritability/
sim_dataset_dir <- paste0("train", n_train, "_test", n_test, 
                          "/tis", n_tissues,
			  "/snps", n_snps, "_eqtls", n_eqtls, "_mteqtls", n_mt_eqtls, "_cor_tissues", cor_tissues,   
                          "/h2_ge_", h2_ge, "_r2_et_", format(r2_et, scientific = FALSE),"/")
dir.create(sim_dataset_dir, recursive = TRUE, showWarnings = FALSE)
setwd(sim_dataset_dir)

### set up dishes
geno_dir <- "geno/"
dir.create(geno_dir, recursive = TRUE, showWarnings = FALSE)
expr_dir <- "gene/"
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
pheno_dir <- "pheno/"
dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)
assoc_dir <- "twas_assoc/"
dir.create(assoc_dir, recursive = TRUE, showWarnings = FALSE)

# remove gene specificity annotation files if it exists from previous simulation
file.remove(paste0(expr_dir, output_prefix, "gene_info.txt"))

# record starting time & system setups
bgt <- Sys.time()
registerDoParallel(numCores)

# set up empty variables to store results
twas_summary <- data.frame()
twas_summary_elnet <- data.frame()
twas_summary_glasso <- data.frame()
gene_tsi <- data.frame()

# parallel run n_genes permutations  
twas_summary <- foreach(i=1:n_genes, .combine = 'rbind') %dopar% {
 
  ### 2.1.2. Genotypic data
  # simulate SNP info file
  snp_info <- simulate_geno(n_snps, n_eqtls, n_mt_eqtls, min_maf = maf_range[1], max_maf = maf_range[2])
  
  # simulate allele dosage file
  ## ref allele dosage
  geno <- simulate_dosage(n_train + n_test, snp_info$MAF)
  rownames(geno) <- c(paste0("train_ind", seq(1,n_train)), paste0("test_ind", seq(1,n_test)))
  colnames(geno) <- snp_info$SNP
  ## standardize genotypes
  geno <- scale(geno)/sqrt(ncol(geno))
 
  # simulate gene expression files across tissues
  # multi-tissue eQTLs
  snp_list <- grep("mt", snp_info$SNP, value = TRUE)
  mt_eqtl_weights <- simulate_eqtls(snp_list, n_tissues, cor_tis = cor_tissues)
  mt_eqtl_signal <- simulate_quan_trait(geno[,snp_list], mt_eqtl_weights, h2_ge*pct_mt_eqtls)
  # single-tissue eQTLs
  snp_list <- grep("st", snp_info$SNP, value = TRUE)
  st_eqtl_weights <- simulate_eqtls(snp_list, n_tissues)
  st_eqtl_signal <- simulate_quan_trait(geno[,snp_list], st_eqtl_weights, h2_ge*(1-pct_mt_eqtls))
  # signal
  expr_signal <- mt_eqtl_signal + st_eqtl_signal
  colnames(expr_signal) <- tissue_names
  # noise
  expr_error <- simulate_eqtls(rownames(geno), n_tissues, cor_tissues) # simulate a N*T expression matrix following N(0, Sigma)
  expr_error <- expr_error/sd(expr_error)*sqrt(1-h2_ge) 
  # gene expression matrix
  expr <- expr_signal + expr_error
  colnames(expr) <- tissue_names

  # simulate phenptypic traits from test dataset
  # for simplicity, assuming a gene is related and only related to a trait in one tissue
  ## simulate signal according to Sudha's summary statistics based TWAS
  trait_signal <- expr_signal[, causal_tis]/sd(expr_signal[, causal_tis]) * sqrt(r2_et)
 
  trait_error <- rnorm(n_train + n_test)
  trait_error <- trait_error/sd(trait_error) * sqrt(1-r2_et)
  trait <- trait_signal + trait_error

  
  ## 2.2. eQTL Detection
  ### Detecting eQTLs from test datasets with 1) elastic net (implemented in PrediXcan) and 2) sparse group lasso (implemented in UTMOST)#
  geno_train <- geno[train_ind, ]
  expr_train <- as.matrix(expr[train_ind, ])
  n_folds <- 5
  ### 2.2.1. elastic net (implemented in PrediXcan)
  elnet_weights <- matrix()
  elnet_weights <- detect_eqtls_elnet(geno_train, expr_train, n_folds)
  colnames(elnet_weights) <- tissue_names
  
  ### 2.2.1. sparse group lasso (implemented in UTMOST)#
  glasso_weights <- matrix()
  glasso_weights <- detect_eqtls_glasso(geno_train, expr_train, n_folds)
  colnames(glasso_weights) <- tissue_names
  rownames(glasso_weights) <- snp_info$SNP
  
  
  ############################################
  ## 3. Entree/Main Course
  ############################################
  
  ## 3.1. Run single-tissue TWAS
  geno_test <- geno[test_ind,]
  trait_test <- trait[test_ind]
  
  ### 3.1.1. Predict gene expression levels from estimated eQTL weights
  elnet_pred_expr <- geno_test %*% elnet_weights
  glasso_pred_expr <- geno_test %*% glasso_weights
  
  ### 3.1.2. Run linear regression on trait and predicted gene expression levels
  ### PrediXcan
  assoc_tmp <- run_lm(trait_test, elnet_pred_expr)
  elnet_tmp <- cbind(data.frame(method = "elnet", gene = rep(paste0("gene", i), n_tissues), trait = rep(paste0("trait_gene",i), n_tissues)), assoc_tmp)
  ### UTMOST
  assoc_tmp <- run_lm(trait_test, glasso_pred_expr)
  glasso_tmp <- cbind(data.frame(method = "glasso", gene = rep(paste0("gene", i), n_tissues), trait = rep(paste0("trait_gene",i), n_tissues)), assoc_tmp)


  ## 3.2. Run integrative TWAS
  # run integrative when there are multiple gene expressing tissues

  ### 3.2.1. Principal Component Regression
  ### MulTiXcan = Elastic Net + Principal Component Regression
  #### for principal component regression, response variable must be centered
  #### selection of the number of principal components
  assoc_tmp <- run_pcr(trait_test, elnet_pred_expr)
  pcr_tmp <- cbind(data.frame(method = "elnet_pcr", gene = paste0("gene", i), trait = paste0("trait_gene",i)), assoc_tmp)
  
  ### 3.2.2. GBJ test
  ### UTMOST
  #### extract the list of eQTLs
  #### calculate covariance matrix of tissues with respect to a certain gene as indicated by the UTMOST paper
  #### cov(tissues) = t(eqtl weights) * cov(genotypes) * eqtl weights
  # assoc_tmp <- run_gbj(trait_test, glasso_pred_expr, geno, glasso_weights, use_snp_cov = TRUE)
  assoc_tmp <- run_gbj(trait_test, glasso_pred_expr, use_snp_cov = FALSE)
  gbj_tmp <- cbind(data.frame(method = "glasso_gbj", gene = paste0("gene", i), trait = paste0("trait_gene",i)), assoc_tmp)

  # combine all test results  
  twas_summary <- rbind(elnet_tmp, glasso_tmp, pcr_tmp, gbj_tmp)

  ############################################
  ## Sides
  ############################################
  # save genotypic data and snp info 
  # save info file
  write.table(snp_info, file = paste0(geno_dir,output_prefix, "snp_info_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  #write.table(geno[train_ind,], file = paste0(geno_dir, output_prefix, "geno_train_gene",i,".txt"), quote = FALSE, sep = " ", col.names = TRUE, row.names = FALSE)
  #write.table(geno[test_ind], file = paste0(geno_dir, output_prefix, "geno_test_gene",i,".txt"), quote = FALSE, sep = " ", col.names = TRUE, row.names = FALSE)

  # save expression and real eQTL info 
  #write.table(rbind(st_eqtl_weights, mt_eqtl_weights), file = paste0(true_eqtl_dir, output_prefix, "true_eqtl_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  #write.table(expr[train_ind,], file = paste0(expr_dir, output_prefix, "expr_train_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  #write.table(expr[test_ind,], file = paste0(expr_dir, output_prefix, "expr_test_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  
  # save trait info
  # complete trait info
  #write.table(trait, file = paste0(pheno_dir, output_prefix, "trait_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

  # save predicted eQTL info
  #write.table(elnet_weights, file = paste0(pred_eqtl_dir, output_prefix, "pred_eqtl_elnet_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  #write.table(glasso_weights, file = paste0(pred_eqtl_dir, output_prefix, "pred_eqtl_glasso_gene",i,".txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

  # save tissue specificity information
  if(file.exists(paste0(expr_dir, output_prefix, "gene_info.txt"))){
    write.table(gene_tsi, paste0(expr_dir, output_prefix, "gene_info.txt"), quote = F, sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
  }else{
    write.table(gene_tsi, paste0(expr_dir, output_prefix, "gene_info.txt"), quote = F, sep = "\t", col.names = TRUE, row.names = FALSE)
  }

  # return value(s) of parallel jobs
  return(twas_summary)  
}
# print running time
edt <- Sys.time()
print(edt-bgt)


parameter_setup <- paste(names(opts),unlist(opts), sep = "=", collapse = ";")
write.table(parameter_setup, paste0(output_prefix, "simulation_parameters.txt"), quote = F, sep = "\t", col.names = FALSE, row.names = FALSE, append = FALSE)

############################################
## Sides
############################################
# save twas assocation results
# sort association results before writing into files

twas_method_list <- unique(twas_summary$method)
for(twas_method in twas_method_list){
  twas_results <- twas_summary[twas_summary$method == twas_method, ]
  twas_results <- twas_results[order(twas_results$raw_p_value), ]
  write.table(twas_results, paste0(assoc_dir, output_prefix, twas_method, "_twas_assoc.txt"), quote = F, sep = "\t", col.names = TRUE, row.names = FALSE)
}


## 3.3. Evaluate Power and/or Type I Error Rate of TWAS Results
### probably in another script
### for your own health, please consider taking small meals



