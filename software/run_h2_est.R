## run_h2_est.R
## Binglan Li 10/27/2018
## estimate trait heritability of simulated expression regulated traits using either GCTA or regress() in R

############################################
## Menu
## 1. Food Preparation
## 1.1. Load Necessary Libraries and Scripts
## 1.2. Define Parameters
## 2. Appetizers
## 2.1. Data Simulation
## 2.2. eQTL Detection
## 3. Entree/Main Course
## 3.1. Calculate trait h2
############################################


############################################
## 1. Food Preparation
############################################

## 1.1. Load Necessary Libraries and Scripts
### load libraries and functions
if(!require(optparse, quietly=T)){install.packages('optparse',repos="http://cran.us.r-project.org"); require(optparse, quietly=T);}
if(!require(mvtnorm, quietly=T)){install.packages('mvtnorm',repos="http://cran.us.r-project.org"); require(mvtnorm, quietly=T);}
suppressMessages(library(regress))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(pls))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
if(!is.loaded("wrapper")){
  #system("R CMD SHLIB optim.c")
  dyn.load("~/group/datasets/UTMOST/CTIMP/optim.so")
  #dyn.load("optim.so")
}
source("twas_simulation_util.R")

## 1.2. Define Parameters
# read external 
opt_list <- list(
  make_option("--method", type="character", default="regress", help="heritability estimation method (default = %default)"),
  make_option("--sample-size", type="integer", default=5000, help="sample size (default = %default)"),
  make_option("--snps", type="integer", default=60, help="Number of SNPs in a given gene (default = %default)"),
  make_option("--eqtls", type="integer", default=30, help="Number of eQTLs among SNPs in a given gene (default = %default)"),
  make_option("--maf", type="character", default="0.01,0.5", help="Range of minor allele frequency (default = %default)"),
  make_option("--pct-mt-eqtls", type="numeric", default=0.8, help="Percentage of multi-tissue eQTLs among total eQTLs (default = %default)"),
  make_option("--genes", type="integer", default=100, help="Number of genes in a round of simulation (default = %default)"),
  make_option("--phenotype", type="integer", default=1, help="Number of phenotypes in a round of simulation (default = %default)"),
  make_option("--expr-tis", type="integer", default=1, help="Number of gene expressing tissues (default = %default)"),
  make_option("--total-tis", type="integer", default=1, help="Number of total tissues (default = %default)"),
  make_option("--h2-ge", type="numeric", default=0.3, help="Heritability of eQTLs and gene expression levels, (default = %default)"),
  make_option("--cor-tissues", type="character", default="0", help="Correlation among tissues, like '0' or '-1,1' (default = %default)"),
  make_option("--r2-et", type="numeric", default=0.01, help="Heritability between gene expression levels and trait, (default = %default)"),
  make_option("--output-dir", type="character", action="store", default=paste0(getwd(), "/simulated_data/"), help="Output directory of simulation results"),
  make_option("--core", type="numeric", default=1, help="Number of cores to run parallel tasks (default = %default)"),
  make_option("--index-gene", type="numeric", default=1, help="Index of genes (default = %default)")
)
opts <- parse_args(OptionParser(option_list=opt_list))
# parse values
h2_method <- opts$method
n_test <- opts$`sample-size` # trait h2 estimation sample size
n_snps <- opts$snps
n_eqtls <- opts$eqtls
maf_range <- as.numeric(strsplit(opts$maf, ",")[[1]])
pct_mt_eqtls <- opts$`pct-mt-eqtls` #pct_mt_eqtls_list <- c(0, 0.2, 0.4, 0.6, 0.8)
n_genes <- opts$genes
n_trait <- opts$phenotype
n_tissues <- opts$`expr-tis` #n_tissue_list <- c(1, 2, 5, 10)
total_tissues <- opts$`total-tis`
h2_ge <- opts$`h2-ge`
cor_tissues <- as.numeric(strsplit(opts$`cor-tissues`, ",")[[1]])
r2_et <- opts$`r2-et` #r2_et_list <- c(0.00001, 0.0005, 0.005, 0.01)
# set random seed
index_gene <- opts$`index-gene`
random_seed <- n_tissues*(n_genes + n_trait) + index_gene
set.seed(random_seed, kind = "L'Ecuyer-CMRG")

output_dir <- opts$`output-dir`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

numCores <- opts$core
print(paste0("numCores=", numCores))

# checkpoints
# no checks at the moment. High-quality food is ready to turn into delicious dishes.

# attain parameters from input variables
# adjust the number of multi-tissue eQTLs according to the number of gene expressing tissues
if(n_tissues == 1){ pct_mt_eqtls <- 0 } # if single tissue, there is no multi-tissue eqtls
if(length(maf_range) == 1){maf_range <- c(maf_range, maf_range)}
n_mt_eqtls <- as.integer(round(n_eqtls * pct_mt_eqtls, digits = 0))
test_ind <- seq(1, n_test) # index for testing individuals
tissue_names <- paste0("Tissue", seq(1, n_tissues)) # tissue names
causal_tis <- tissue_names[1] # assuming only one tissue is causal

############################################
## 2. Appetizers
############################################

## 2.1. Data Simulation

### 2.1.1. define paths and parameters
# i.e. path = sample_size/tissue/snp_composition/heritability/
sim_dataset_dir <- paste0("sample_size", n_test, 
                          "/tis", n_tissues,
			  "/snps", n_snps, "_eqtls", n_eqtls, "_maf", maf_range[1], "_", maf_range[2],   
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

# remove gene specificity annotation files if it exists from previous simulation
if(file.exists(paste0(expr_dir, "gene_info.txt"))){
  file.remove(paste0(expr_dir, "gene_info.txt"))
}

# record starting time & system setups
bgt <- Sys.time()
registerDoParallel(numCores)


### 2.1.2. Genotypic data
# simulate SNP info file
snp_info <- simulate_geno(n_snps, n_eqtls, n_mt_eqtls, min_maf = maf_range[1], max_maf = maf_range[2])
  
# simulate allele dosage file
## ref allele dosage
geno <- simulate_dosage(n_test, snp_info$MAF)
rownames(geno) <- paste0("test_ind", seq(1,n_test))
colnames(geno) <- snp_info$SNP
eqtl_index <- grep("eqtl", snp_info$SNP, ignore.case = T)
## generate plink ped and map file for gcta
if(h2_method == "gcta"){
  geno_ped <- generate_ped(geno, snp_info)
  geno_map <- data.frame(chr = 1, id = snp_info$SNP, distance = 0, pos = seq(1, dim(snp_info)[1])*100)
}
 
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
  
# save files for trait h2 estimation for gcta
if(h2_method == "gcta"){
  # save genotypic data and snp info 
  # save info file
  write.table(geno_ped, file = paste0(geno_dir, "gene",index_gene,".ped"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(geno_map, file = paste0(geno_dir, "gene",index_gene,".map"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

  # generate GRM for genotype file  
  geno_file_prefix <- paste0(geno_dir, "gene",index_gene)
  system(paste0("plink --file ", geno_file_prefix, " --make-bed --out ", geno_file_prefix))
  system(paste0("gcta64 --bfile ", geno_file_prefix, " --make-grm --out ", geno_file_prefix, "_grm"))
}

# simulate different traits from the same genotype data
## to obtain the heritability variance driven by simulated phenotype
foreach(i=1:n_trait) %dopar% {
  # for simplicity, assuming a gene is related and only related to a trait in one tissue
  ## simulate signal according to Sudha's summary statistics based TWAS
  trait_signal <- expr_signal[, causal_tis]/sd(expr_signal[, causal_tis]) * sqrt(r2_et)

  trait_error <- rnorm(n_test)
  trait_error <- trait_error/sd(trait_error) * sqrt(1-r2_et)
  trait <- trait_signal + trait_error

  ############################################
  ## 3. Entree/Main Course
  ############################################
    
  # estimate heritability using GCTA REML
  if(h2_method == "gcta"){
    # 3.1.1. save files for trait h2 estimation
    # save trait info
    # for trait heritability estimation
    trait_file_prefix <- paste0("trait_gene",index_gene, "_pheno", i, ".txt")
    trait_h2 <- data.frame(fid = paste0("test_ind", seq(1,n_test)), iid = paste0("test_ind", seq(1,n_test)), trait = trait)
    write.table(trait_h2, file = paste0(pheno_dir, trait_file_prefix), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)  

    # 3.1.2. calculate trait h2 using external gcta
    ## trait heritability estimation
    trait_file <- paste0(pheno_dir, trait_file_prefix)
    system(paste0("gcta64 --reml --grm ", geno_file_prefix, "_grm --pheno ", trait_file, " --grm-adj 0 --out ", geno_file_prefix, "_pheno", i, "_reml_h2")) 

    ## remove intermediate files
    system(paste0("rm ", trait_file))
  }

  # estimate heritability using R regress() function 
  if(h2_method == "regress"){    
    # estimate heritability using R regress function
    geno <- scale(geno)/sqrt(ncol(geno))
    GRM <- tcrossprod(geno)
    fmAlt <- regress(trait ~ 1, ~GRM)
    trait_h2 <- fmAlt$sigma[1]/sum(fmAlt$sigma)
    h2_file_prefix <- paste0(geno_dir, "gene",index_gene, "_pheno", i, "_regress_h2.txt")
    write.table(trait_h2, h2_file_prefix, quote = F, sep = "\t", col.names = FALSE, row.names = FALSE)
  }
}

# print running time
edt <- Sys.time()
print(edt-bgt)
#system("rm geno/*.b* geno/*.f* geno/*grm* geno/*map geno/*ped geno/*nosex geno/*log pheno/*")




