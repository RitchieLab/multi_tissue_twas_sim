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
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(pls))
suppressMessages(library(GBJ))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
#suppressMessages(library(hudson))
if(!is.loaded("wrapper")){
  #system("R CMD SHLIB optim.c")
  dyn.load("~/group/datasets/UTMOST/CTIMP/optim.so")
  #dyn.load("~/Box/Personal_Document/twas_simulation/analysis/simulation_script/CTIMP/optim.so")
}
source("twas_simulation_util_sum_stat.R")

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
  make_option("--expr-tis", type="integer", default=10, help="Number of gene expressing tissues (default = %default)"),
  make_option("--total-tis", type="integer", default=10, help="Number of total tissues (default = %default)"),
  make_option("--h2-ge", type="numeric", default=0.3, help="Heritability of eQTLs and gene expression levels, (default = %default)"),
  make_option("--cor-tissues", type="character", default="0", help="Correlation among tissues, like '0' or '-1,1' (default = %default)"),
  make_option("--r2-et", type="numeric", default=0.01, help="Heritability between gene expression levels and trait, (default = %default)"),
  make_option("--r2-cov", type="numeric", default=0.05, help="Heritability between gene expression levels and trait, (default = %default)"),
  make_option("--output-dir", type="character", action="store", default="~/group/personal/binglan/twas_simulation/results/test_indlvl_vs_sumstat/", help="Output directory of simulation results"),
  make_option("--output-prefix", type="character", action="store", default="", help="Prefix of each output file"),
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
r2_cov <- opts$`r2-cov`
# set random seed
random_seed <- n_tissues*100 + opts$`random-seed`
set.seed(random_seed, kind = "L'Ecuyer-CMRG")

output_dir <- opts$`output-dir`
output_prefix <- opts$`output-prefix`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)


# checkpoints
# no checks at the moment. High-quality food is ready to turn into delicious dishes.

# attain parameters from input variables
# adjust the number of multi-tissue eQTLs according to the number of gene expressing tissues
if(n_tissues == 1){ pct_mt_eqtls <- 0 } # if single tissue, there is no multi-tissue eqtls
if(length(output_prefix) > 1){output_prefix <- paste0(output_prefix, "_")} # reformat prefix if given
if(length(maf_range) == 1){maf_range <- c(maf_range, maf_range)}
n_mt_eqtls <- as.integer(round(n_eqtls * pct_mt_eqtls, digits = 0))
train_ind <- seq(1, n_train) # index for training individuals
test_ind <- seq(n_train + 1, n_train + n_test) # index for testing individuals
tissue_names <- paste0("Tissue", seq(1, n_tissues)) # tissue names
causal_tis <- tissue_names[1] # assuming only one tissue is causal
#trait_es <- rnorm(n = 1) # effect size for all gene-trait relationships in one iteration of simulation


############################################
## 2. Appetizers
############################################

## 2.1. Data Simulation
# record starting time & system setups
bgt <- Sys.time()

# set up empty variables to store results
all_assoc_indlvl <- data.frame()
all_assoc_sumstat <- data.frame()

# parallel run n_genes permutations  
for(i in 1:n_genes){
#for(i in 1:5){
  
  ### 2.1.2. Genotypic data
  # simulate SNP info file
  snp_info <- simulate_geno(n_snps, n_eqtls, n_mt_eqtls, min_maf = maf_range[1], max_maf = maf_range[2])
  
  # simulate allele dosage file
  ## ref allele dosage
  geno <- simulate_dosage(n_train + n_test, snp_info$MAF)
  rownames(geno) <- c(paste0("train_ind", seq(1,n_train)), paste0("test_ind", seq(1,n_test)))
  colnames(geno) <- snp_info$SNP
  #geno <- scale(geno)/sqrt(ncol(geno))
 
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
  ## simulate trait signal with a defined effect size
  #trait_signal <- expr[, tissue_names[1]] * trait_es
  #trait_signal <- trait_signal/sd(trait_signal) * sqrt(r2_et)
  
  ## simulate signal according to Sudha's summary statistics based TWAS
  # which should be the E?
  #trait_signal <- expr[, tissue_names[1]]/sd(expr[, tissue_names[1]]) * sqrt(r2_et)
  trait_signal <- expr_signal[, causal_tis]/sd(expr_signal[, causal_tis]) * sqrt(r2_et)
  
  trait_cov <- rnorm(n_train + n_test)
  trait_cov <- trait_cov/sd(trait_cov) * sqrt(r2_cov)
 
  trait_error <- rnorm(n_train + n_test)
  trait_error <- trait_error/sd(trait_error) * sqrt(1-r2_et-r2_cov)
  trait <- trait_cov + trait_signal + trait_error

  
  ## 2.2. eQTL Detection
  ### Detecting eQTLs from test datasets with 1) elastic net (implemented in PrediXcan) and 2) sparse group lasso (implemented in UTMOST)#
  geno_train <- geno[train_ind, ]
  expr_train <- as.matrix(expr[train_ind, ])
  n_folds <- 5
  
  ### 2.2.1. sparse group lasso (implemented in UTMOST)#
  glasso_weights <- matrix()
  glasso_weights <- detect_eqtls_glasso(geno_train, expr_train, n_folds)
  colnames(glasso_weights) <- tissue_names
  rownames(glasso_weights) <- snp_info$SNP

  
  ############################################
  ## 3. Entree/Main Course
  ############################################
  
  ## 3.1. Run UTMOST
  geno_test <- geno[test_ind,]
  trait_test <- trait[test_ind]
  cov_test <- trait_cov[test_ind]
  
  
  ## 3.2 Run individual level TWAS
  ### Predict gene expression levels from estimated eQTL weights
  glasso_pred_expr <- geno_test %*% glasso_weights
 
  ### Run linear regression on trait and predicted gene expression levels
  ### UTMOST
  tmp <- calc_score_stats(null_model = glm(trait_test ~ cov_test), factor_matrix = glasso_pred_expr, link_function = "linear")
  gbj_result <- GBJ(test_stats = tmp$test_stats, cor_mat = tmp$cor_mat)
  gbj_tmp_ind_lvl <- cbind(SNP = paste0("gene", i), CHR = 'Sim', POS = i, pvalue = gbj_result$GBJ_pvalue)
  all_assoc_indlvl <- rbind(all_assoc_indlvl, gbj_tmp_ind_lvl)
  
  
  ## 3.2. Run summary stat TWAS 
  # run integrative when there are multiple gene expressing tissues
  #### calculate covariance matrix of tissues with respect to a certain gene as indicated by the UTMOST paper
  #### for summary stat, cov(tissues) = t(eqtl weights) * cov(genotypes) * eqtl weights, followed by normalization
  gwas <- run_lm(trait_test, geno_test, cov_test)
  gene_lvl_z <- vector()
  for(tissue_i in 1:n_tissues){
    cov <- sqrt(var(geno_test)/var(glasso_pred_expr[,tissue_i]))
    cov[upper.tri(cov)] <- 0
    cov[lower.tri(cov)] <- 0
    gene_lvl_z[tissue_i] <- t(glasso_weights[,tissue_i]) %*% cov %*% gwas$z_score
  }
  # the following is equivalent to, if not the same as, "assoc_tmp_sum_stat <- run_gbj(twas_z = glasso_tmp$z_score, genotypes = geno, eqtl_weights = glasso_weights, gwas_sum_stat = TRUE)"
  assoc_tmp_sum_stat <- run_gbj(twas_z = gene_lvl_z, genotypes = geno, eqtl_weights = glasso_weights, gwas_sum_stat = TRUE)
  gbj_tmp_sum_stat <- cbind(SNP = paste0("gene", i), CHR = 'Sim', POS = i, pvalue = assoc_tmp_sum_stat$raw_p_value)
  all_assoc_sumstat <- rbind(all_assoc_sumstat, gbj_tmp_sum_stat)
}

# print running time
edt <- Sys.time()
print(edt-bgt)

# store output
write.table(all_assoc_indlvl, paste0(output_prefix, "utmost_sim_ind_lvl.txt"), quote = F, sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)
write.table(all_assoc_sumstat, paste0(output_prefix, "utmost_sim_sum_stat.txt"), quote = F, sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE)

# draw hudson plot to compare results
#Bonferroni_sig_threshold <- 0.05
#gmirror(top = all_assoc_indlvl, bottom = all_assoc_sumstat, 
#        tline = Bonferroni_sig_threshold, bline = Bonferroni_sig_threshold, 
#        toptitle = "Simulated Individual-level TWAS (covariates explained 5% trait variance)", 
#        bottomtitle = "Simulated Summary-based TWAS (covariates explained 5% trait variance)", 
#        highlight_p = Bonferroni_sig_threshold,
#        annotate_p = Bonferroni_sig_threshold,
#        file = paste0(output_prefix, "utmost_sim_indlvl_vs_sumstat_twas.png"))


