## twas_simulation_util.R
## Binglan Li 11/06/2018
## necessary functions for twas simulations

############################################
## Grocery List
## 1. simulate_dosage
## 2. simulate_geno
## 3. calculate_maf
## 4. simulate_eqtls
## 5. simulate_quan_trait
## 6. generate_fold_ids
## 7. detect_eqtls_elnet
## 8. elastic_net_mse
## 9. minmax_lambda
## 10. glasso
## 11. glasso_no_early_stopping
## 12. multi_mse
## 13. detect_eqtls_glasso
## 14. run_lm
## 15. count_sig_assoc
## 16. calculate_tau
## 17. run_pcr
## 18. run_gbj
## 19. generate_ped
############################################


############################################
## 1. simulate_dosage
############################################

simulate_dosage <- function(n, maf){
  # simulate ref/major allele dosages
  # MAF stands for minor allele frequency
  dosage <- matrix(nrow = n, ncol = length(maf))
  for(i in 1:length(maf)){
    prob_1 = maf[i]^2 # prob of homogeneous alternative alleles
    prob_2 = 2 * maf[i] * (1-maf[i]) # prob of heterogyzous alleles
    prob_3 = (1-maf[i])^2
    dosage[,i] <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(prob_1, prob_2, prob_3))
  }
  return(dosage)
}

############################################
## 2. simulate_geno
############################################

simulate_geno <- function(n_snps, n_eqtls = 0, n_mt_eqtls = 0, min_maf = 0.01, max_maf = 0.5){
  # simulate SNP info files containing SNP IDs, minor allele frequency, ref and alt alleles
  stopifnot(n_snps != 0, min_maf <= max_maf)
  if(n_eqtls == 0 ){
    snp_info <- data.frame(SNP = paste0("SNP", seq(1,n_snps)), stringsAsFactors = FALSE)
  }else if(n_mt_eqtls == 0){
    snp_info <- data.frame(
      SNP = c(paste0("st_eQTL", seq(1,n_eqtls)), 
              paste0("SNP", seq(1,n_snps - n_eqtls))),
      stringsAsFactors = FALSE
    )
  }else{
    snp_info <- data.frame(
      SNP = c(paste0("mt_eQTL", seq(1,n_mt_eqtls)), 
              paste0("st_eQTL", seq(1,n_eqtls - n_mt_eqtls)), 
              paste0("SNP", seq(1,n_snps - n_eqtls))),
      stringsAsFactors = FALSE)    
  }
  snp_info$MAF <- runif(n_snps, min = min_maf, max = max_maf)
  snp_info$REF <- sample(c("A", "T", "C", "G"), size = n_snps, replace = TRUE)
  snp_info$ALT <- apply(snp_info, 1, function(x){sample(setdiff(c("A", "T", "C", "G"),x["REF"]), 1)})
  return(snp_info)
}

############################################
## 3. calculate_maf
############################################

calculate_maf <- function(dosage_matrix){
  N <- nrow(dosage_matrix)
  maf <- 1 - colMeans(dosage_matrix)/2
  return(maf)
}

############################################
## 4. simulate_eqtls
############################################

simulate_eqtls <- function(snp_list, n_tis = 1, cor_tis = 0){
  # simulate gene expression data based on provided genotypic data
  stopifnot(n_tis>0)
  
  n_snps <-  length(snp_list)
  if(n_snps == 0){
    snp_weights <- matrix(0, ncol = n_tis, nrow = n_snps)
    colnames(snp_weights) <- paste0("Tissue", seq(1, n_tis))
    return(snp_weights)
    stop
  }
  # simulate weights
  if(n_tis == 1) {
    snp_weights <- as.matrix(rnorm(n_snps, mean = 0, sd = 1))
  }else{
    # simulate weight in multiple genes
    # set sigma matrix
    if(length(cor_tis == 1)){
      # cor(tissue1, tissue2) = cov(tissue1, tissue2)/(var(tissue1) * var(tissue2))
      tis_sigma <- matrix(cor_tis, nrow = n_tis, ncol = n_tis)
      diag(tis_sigma) <- 1
    }else if(length(cor_tis == 2)){
      tis_sigma <- matrix(runif(n_tis*n_tis, min = cor_tis[1], max = cor_tis[2]), nrow = n_tis, ncol = n_tis)
      diag(tis_sigma) <- 1
    }else{
      print("Parameter out of subscription.")
    }
    snp_weights <- rmvnorm(length(snp_list), mean = rep(0, n_tis), sigma = tis_sigma)
  }
  rownames(snp_weights) <- snp_list
  colnames(snp_weights) <- paste0("Tissue", seq(1, n_tis))
  return(snp_weights)
}

############################################
## 5. simulate_quan_trait
############################################

simulate_quan_trait <- function(x, betas, h2){
  # simulate quantitative trait of x*betas with mean = 0 and var = h2
  if(length(betas) == 0){
    quan_trait <- matrix(0, nrow = nrow(x), ncol = ncol(betas))
    return(quan_trait)
  }else{
    quan_trait <- x %*% betas
    quan_trait <- apply(quan_trait, 2, function(x){x/sd(x)*sqrt(h2)})# for standardized genotypes
    #quan_trait <- scale(quan_trait)*sqrt(h2) # for non-standardized genotypes
    return(quan_trait)
  }
}

############################################
## 6. generate_fold_ids
############################################

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

############################################
## 7. detect_eqtls_elnet
############################################

detect_eqtls_elnet <- function(geno, gene_expr, n_folds = 5){
  # identify eQTLs using elastic net
  # expression data are expected to be a N by T matrix, where N represents sample size of observations, T represents tissue numbers
  
  #bgt <- Sys.time()
  # check and set up parameters
  stopifnot(nrow(geno) == nrow(gene_expr))
  n <- nrow(gene_expr)
  n_tis <- ncol(gene_expr)
  elnet_weights <- matrix(nrow = ncol(geno), ncol = n_tis)
  
  # perform elastic net
  for(tis_ind in 1:n_tis){
    cv_fold_ids <- generate_fold_ids(n, n_folds)
    elnet_fit <- tryCatch(cv.glmnet(geno, gene_expr[,tis_ind], nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE, standardize = FALSE),
                          error = function(cond) {message('Error'); message(geterrmessage()); list()})
    best_lam_ind <- which.min(elnet_fit$cvm)
    #pred_eqtl_ind <- which(elnet_fit$glmnet.fit$beta[,best_lam_ind] != 0)
    #elnet_weights[[tis_ind]] <- elnet_fit$glmnet.fit$beta[pred_eqtl_ind, best_lam_ind]
    elnet_weights[, tis_ind] <- elnet_fit$glmnet.fit$beta[, best_lam_ind]
  }
  rownames(elnet_weights) <- colnames(geno)

  
  # print running time
  #edt <- Sys.time()
  #print(edt-bgt)
  
  return(elnet_weights)  
}

############################################
## 8. elastic_net_mse
############################################

elastic_net_mse <- function(lst, X_tune, Y_tune, X_test, Y_test){
  ## adapted from CTIMP (https://github.com/yiminghu/CTIMP), copyright belongs to Dr. Yuming Hu
  ## evaluate the performance of elastic net on each response
  ## args
  ## lst: a list of glmnet object (fitted elastic net model for each response)
  ## X_tune: a matrix of covariate matrices corresponding for each response (for tuning lambda)
  ## Y_tune: a matrix of response vectors (for tuning lambda)
  ## X_test: a matrix of covariate matrices corresponding for each response (for testing performance)
  ## Y_test: a matrix of response vectors (for testing performance)
  ## value
  ## lam: best performing lambda (on (X_tune,Y_tune)) for each response
  ## mse: list of matrices with each element being a matrix of predicted vs observed response
  ## est: estimated effect sizes for each response (B matrix)
  P <- length(lst) # number of tissues
  M <- ncol(X_tune) # number of SNPs
  lam_V <- rep(0, P)
  test_res <- list()
  test_beta <- matrix(0, M, P)
  for(t in 1:P){
    ncv <- length(lst[[t]]$lambda)
    tmp_mse <- rep(0, ncv)
    for(k in 1:ncv){
      tmp_mse[k] <- mean((Y_tune[,t] - X_tune%*%lst[[t]]$glmnet.fit$beta[,k])^2)
    }
    ss <- which.min(tmp_mse)
    test_beta[,t] <- lst[[t]]$glmnet.fit$beta[,ss]
    lam_V[t] <- lst[[t]]$lambda[ss]
    predicted <- X_test%*%lst[[t]]$glmnet.fit$beta[,ss]
    test_res[[t]] <- cbind(Y_test[,t], predicted)
  }
  list(lam = lam_V, mse = test_res, est = test_beta)
}

############################################
## 9. minmax_lambda
############################################

minmax_lambda <- function(lst){
  ## adapted from CTIMP (https://github.com/yiminghu/CTIMP), copyright belongs to Dr. Yuming Hu
  ## get the minimum and maximum of lambda searched in cross-validation of an elastic net model
  ## args
  ## lst: an object returned by glmnet
  ## value
  ## min_lam: smallest lambda searched in glmnet cross-validation
  ## max_lam: largest lambda searched in glmnet cross-validation
  max_lam <- max(unlist(lapply(lst, function(x){max(x$lambda)})))
  min_lam <- min(unlist(lapply(lst, function(x){min(x$lambda)})))
  c(min_lam, max_lam)
}

############################################
## 10. glasso
############################################

glasso <- function(X, Y, X1, Y1, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3){
  ## adapted from CTIMP (https://github.com/yiminghu/CTIMP), copyright belongs to Dr. Yuming Hu
  
  M <- nrow(XY) # number of SNPs
  P <- ncol(Y) # number of tissues
  NN <- rep(nrow(X), P)  # sample size in each tissue
  old_objV1 = rep(0,P)
  for(t in 1:P){
    old_objV1[t] <- 1/2*mean((Y[,t]-X%*%theta[,t])^2)
  }
  old_objV2 = rep(0,P)
  for(t in 1:P){
    old_objV2[t] <- 1/2*mean((Y1[,t]-X1%*%theta[,t])^2)
  }
  beta_j_lasso <- rep(0, P)
  tmp_XYj <- 0
  if(!is.loaded("wrapper")){
    dyn.load("optim.so")
  }
  for(i in 1:maxiter){
    # adjust estimated theta (betas)
    res <- .Call("wrapper", rep(list(XX), P), XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
    new_objV1 <- new_objV2 <- rep(0,P)
    for(t in 1:P){
      new_objV1[t] <- 1/2*mean((Y[,t]-X%*%theta[,t])^2)
    }
    for(t in 1:P){
      new_objV2[t] <- 1/2*mean((Y1[,t]-X1%*%theta[,t])^2)
    }
    #cat("Testing error: ", new_objV2, '\n')
    if(mean(new_objV2) > mean(old_objV2)|mean(new_objV1) > mean(old_objV1)){
      break
    }else{
      old_objV2 = new_objV2
    }
    if(max(abs(new_objV1-old_objV1)) < eps){
      break
    }else{
      old_objV1 <- new_objV1
    }
  }
  list(est = theta, avg_tune_err = mean(new_objV2), tune_err=new_objV2)
}

############################################
## 11. glasso_no_early_stopping
############################################

glasso_no_early_stopping <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3){
  ## adapted from CTIMP (https://github.com/yiminghu/CTIMP), copyright belongs to Dr. Yuming Hu
  
  M <- nrow(XY) # number of SNPs
  P <- ncol(Y) # number of tissues
  NN <- rep(nrow(X), P)  # sample size in each tissue
  old_objV1 <- rep(0,P)
  for(t in 1:P){
    old_objV1[t] <- 1/2*mean((Y[,t]-X%*%theta[,t])^2)
  }
  #cat("Training error: ", mean(old_objV1), '\n')
  beta_j_lasso <- rep(0, P)
  tmp_XYj <- 0
  if(!is.loaded("wrapper")){
    dyn.load("optim.so")
  }
  for(i in 1:maxiter){
    res <- .Call("wrapper", rep(list(XX), P), XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
    new_objV1 <- rep(0,P)
    for(t in 1:P){
      new_objV1[t] <- 1/2*mean((Y[,t]-X%*%theta[,t])^2)
    }
    #cat("Training error: ", mean(new_objV1), '\n')
    if(max(abs(new_objV1-old_objV1)) < eps|mean(new_objV1) > mean(old_objV1)){
      break
    }else{
      old_objV1 <- new_objV1
    }
  }
  list(est = theta, avg_train_err = mean(new_objV1), train_err = new_objV1)
}

############################################
## 12. multi_mse
############################################

multi_mse <- function(theta_est, X_test, Y_test){
  ## adapted from CTIMP (https://github.com/yiminghu/CTIMP), copyright belongs to Dr. Yuming Hu
  
  answer <- list()
  P <- ncol(theta_est)
  for(t in 1:P){
    predicted <- X_test%*%theta_est[,t]
    answer[[t]] <- cbind(Y_test[,t], predicted)
  }
  answer
}

############################################
## 13. detect_eqtls_glasso
############################################

detect_eqtls_glasso <- function(geno, gene_expr, n_folds = 5){
  # identify eQTLs using sparse group lasso
  # adapted from UTMOST CTIMP
  
  #bgt <- Sys.time()
  
  stopifnot(nrow(geno) == nrow(gene_expr))
  n <- nrow(geno)
  n_tis <- ncol(gene_expr)
  n_snps <- ncol(geno)
  ntune <- 100
  
  # tuning parameters from single tissue-based elastic net
  single_res_test <- list()
  single_lam <- matrix(0, n_folds, n_tis)
  single_theta_est <- list()
  # tunning parameters for integrative tissue-based group lasso
  multi_res_test <- list()
  multi_lam <- matrix(0, n_folds, 2)
  multi_theta_est <- list()
  
  res_tune <- list()
  rec_lamv <- matrix(0, n_folds, ntune)
  
  # centralize dosage
  geno_centered <- apply(geno, 2, function(x){x - mean(x)})

  # obtain diagnal of X transpose X
  XX <- t(geno_centered)%*%geno_centered/n_train
  Xnorm <- diag(XX)
  
  cv_fold_ids <- generate_fold_ids(n, n_folds)
  for(f in 1:n_folds){
    test_id <- which(cv_fold_ids == f)
    tune_id <- which(cv_fold_ids == f%%n_folds + 1)
    train_id <- intersect(which(cv_fold_ids != f), which(cv_fold_ids != f%%n_folds + 1))
    
    X_test <- geno_centered[test_id,]
    Y_test <- as.matrix(gene_expr[test_id,])
    X_tune <- geno_centered[tune_id,]
    Y_tune <- as.matrix(gene_expr[tune_id,])
    X_train <- geno_centered[train_id,]
    Y_train <- as.matrix(gene_expr[train_id,])
    
    ## model training
    ## train elastic net and used average lambda as tuning parameters
    single_initial_est <- matrix(0, nrow = ncol(X_train), ncol = n_tis)
    single_summary <- list()
    for(t in 1:n_tis){
      el_model <- cv.glmnet(X_train, Y_train[,t], alpha = 0.5, nfolds = 5, standardize = FALSE)
      single_summary[[t]] <- el_model
      single_initial_est[,t] <- el_model$glmnet.fit$beta[,which.min(el_model$cvm)]
    }
    ## performance of Elastic net on tuning and testing data with various tuning parameters
    els_output <- elastic_net_mse(single_summary, X_tune, Y_tune, X_test, Y_test)
    single_res_test[[f]] <- els_output$mse
    single_lam[f,] <- els_output$lam
    single_theta_est[[f]] <- els_output$est
    remove(els_output)
    
    ## use elastic net ests row norm as weights
    lam_range <- minmax_lambda(single_summary)
    sig_norm <- apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
    sig_norm[sig_norm==0] <- rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
    sig_norm <- sig_norm/sum(sig_norm)
    weights2 <- 1/sig_norm; weights2 <- weights2/sum(weights2);
    
    tis_norm <- apply(single_initial_est, 2, function(x){sum(abs(x))})
    tis_norm[tis_norm==0] <- rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
    tis_norm <- tis_norm/sum(tis_norm)
    weights1 <- 1/tis_norm; weights1 <- weights1/sum(weights1);
    lam_V <- seq(lam_range[1], lam_range[2], length.out = ntune)
    
    initial_numeric <- as.numeric(single_initial_est)
    remove(single_summary); remove(single_initial_est);
    
    ## preparation
    XY <- t(X_train)%*%Y_train/nrow(X_train)
    XX_train <- t(X_train)%*%X_train/nrow(X_train)
    spsz <- rep(n, n_tis)  # record sample size in each tissues
    res_tune[[f]] <- array(-1, dim=c(ntune, ntune, n_tis))
    rec_lamv[f,] <- lam_V
    for(lam1 in 1:ntune){
      for(lam2 in 1:ntune){
        single_est <- matrix(initial_numeric, n_snps, n_tis)
        ans <- glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est)
        if(sum(ans$est!=0)>0){
          res_tune[[f]][lam1,lam2, ] <- ans$tune_err
          remove(single_est); remove(ans);
        }else{
          remove(single_est); remove(ans);
          break
        }			
      }
    }
    avg_tune_res <- apply(res_tune[[f]], c(1,2), mean)
    best.lam <- which(avg_tune_res == min(avg_tune_res[avg_tune_res>=0]), arr.ind = TRUE)[1,]
    single_est <- matrix(initial_numeric, n_snps, n_tis)
    ans <- glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_est)
    multi_res_test[[f]] <- multi_mse(ans$est, X_test, Y_test)
    multi_lam[f,] <- lam_V[best.lam]
    multi_theta_est[[f]] <- ans$est
    remove(single_est); remove(ans);
  }
  ## generate an estimate with whole data
  # initial values 
  single_initial_est <- matrix(0, n_snps, n_tis)
  for(t in 1:n_tis){
    el_model <- cv.glmnet(geno_centered, gene_expr[,t], alpha = 0.5, nfolds = 5, standardize = FALSE)
    single_initial_est[,t] = el_model$glmnet.fit$beta[,which.min(el_model$cvm)]
  }
  
  sig_norm <- apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
  sig_norm[sig_norm==0] <- rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
  sig_norm <- sig_norm/sum(sig_norm)
  weights2 <- 1/sig_norm; weights2 <- weights2/sum(weights2);
  
  tis_norm <- apply(single_initial_est, 2, function(x){sum(abs(x))})
  tis_norm[tis_norm==0] <- rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
  tis_norm <- tis_norm/sum(tis_norm)
  weights1 <- 1/tis_norm; weights1 <- weights1/sum(weights1);
  
  spsz <- rep(n, n_tis)  # record sample size in each tissues
  initial_numeric <- as.numeric(single_initial_est)
  remove(single_initial_est)
  XY <- t(geno_centered)%*%gene_expr/nrow(geno_centered)
  tmp_res <- rep(0, n_folds)
  for(f in 1:n_folds){
    ans <- glasso_no_early_stopping(X=geno_centered, Y=gene_expr, XX=XX, XY=XY, Xnorm=Xnorm, lambda1=multi_lam[f,1]/spsz, lambda2=multi_lam[f,2], theta=matrix(initial_numeric,n_snps,n_tis))
    tmp_res[f] <- ans$avg_train_err
  }
  final.lam <- multi_lam[which.min(tmp_res),]
  ans <- glasso_no_early_stopping(X=geno_centered, Y=gene_expr, XX=XX, XY=XY, Xnorm=Xnorm, lambda1=final.lam[1]/spsz, lambda2=final.lam[2], theta=matrix(initial_numeric,n_snps,n_tis))
  
  # print running time
  #edt <- Sys.time()
  #print(edt-bgt)
  
  # return estimated weights
  return(ans$est)
}

############################################
## 14. run_lm
############################################

run_lm <- function(trait, gene_expr){
  # run linear regression on a simulated trait and predicted gene expression levels
  # return a data frame of information about tissue, effect size, and p-value
  n_tis <- ncol(gene_expr)
  tis_names <- colnames(gene_expr)
  assoc_summary <- data.frame()
  for(i in 1:n_tis){
    ge <- gene_expr[,i]
    if(all(ge == 0)){
      lm_result <- data.frame(tissue = tis_names[i], effect_size = NA, raw_p_value = NA)
    }else{
      lm_summary <- summary(lm(trait ~ ge))
      lm_result <- data.frame(
        tissue = tis_names[i],
        effect_size = lm_summary$coefficients["ge",1], 
        raw_p_value = lm_summary$coefficients["ge",4])
    }
    assoc_summary <- rbind(assoc_summary, lm_result)
  }
  return(assoc_summary)
}

############################################
## 15. count_sig_assoc
############################################

count_sig_assoc <- function(p_values, alpha = 0.05, p_adjust = "none"){
  # count the number of significant p-values among the provided
  n_pvalue <- length(p_values)
  p_values <- p.adjust(p_values, method = p_adjust)
  sig_percent <- length(which(p_values <= alpha))/n_pvalue
  return(sig_percent)
}

############################################
## 16. calculate_tau
############################################

calculate_tsi <- function(expr_lvls = NA, genotypes = NA, eqtl_weights = NA, total_tis = 10, expr_tis = 1, method = "tau", use_geno_cov = FALSE){
  ############################################
  ## use all gene expression levels
  ## not appropriate in this simulation
  ############################################
  #  expr_tis <- ncol(expr_lvls)
  #  noexpr_tis <- total_tis - expr_tis
  #  n_sample <- nrow(expr_lvls)
  #  expr_max <- max(expr_lvls)
  #  expr_lvls <- expr_lvls - min(expr_lvls)
  #  # calculate tau
  #  gene_tau <- (sum(1-expr_lvls/expr_max) + noexpr_tis*n_sample)/(total_tis*n_sample - 1)
  
  if(method == "tau"){
    ############################################
    ## use average expression levels of genes in each tissue
    ## not appropriate in this simulation
    ## because all genes are expected to have mean 0 according to statistical models
    ############################################
    stopifnot(length(expr_lvls) > 0)
    
    expr_tis <- ncol(expr_lvls)
    if(expr_tis == 1){
      gene_tau <- 1
      return(gene_tau)
    }else{
      noexpr_tis <- total_tis - expr_tis
      expr_means <- colMeans(expr_lvls)
      expr_means <- expr_means - min(expr_means)
      expr_max <- max(expr_means)
      gene_tau <- (sum(1-expr_means/expr_max) + noexpr_tis)/(total_tis - 1)
      return(gene_tau)
    }
  }else if(method == "expr_var" && use_geno_cov == FALSE){
    ############################################
    ## use similarities of a gene's expression levels in multiple tissues
    ## 0 = tissue specific
    ## 1 = ubiquitously expressed gene
    ############################################
    stopifnot(length(expr_lvls) > 0)
    
    expr_tis <- ncol(expr_lvls)
    if(expr_tis == 1){
      tis_spec_score <- 1
      return(tis_spec_score)
    }else{
      no_expr_pairs <- total_tis*(total_tis + 1)/2 - expr_tis*(expr_tis + 1)/2
      cor_expr <- cor(expr_lvls)[upper.tri(matrix(ncol = expr_tis, nrow = expr_tis), diag = T)]
      tis_spec_score <- (sum(1 - cor_expr) + no_expr_pairs)/(total_tis*(total_tis + 1)/2)
      return(tis_spec_score)
    }
  }else if(method == "expr_var" && use_geno_cov == TRUE){
    ############################################
    ## use similarities of a gene's expression levels in multiple tissues
    ## use detected eQTLs to evaluate expression similarities among tissues
    ## 0 = tissue specific
    ## 1 = ubiquitously expressed gene
    ############################################
    stopifnot(length(genotypes) > 0, length(eqtl_weights) > 0)
    
    expr_tis <- ncol(eqtl_weights)
    if(expr_tis == 1){
      tis_spec_score <- 1
      return(tis_spec_score)
    }else{
      # originally designed by UTMOST for summary-level statistics
      eqtl_list <- rownames(eqtl_weights)[!apply(eqtl_weights, 1, function(x){all(x == 0)})]
      cov_snps <- cov(genotypes[, eqtl_list])
      eqtl_weights <- eqtl_weights[eqtl_list, ]
      cor_expr <- t(eqtl_weights) %*% cov_snps %*% eqtl_weights
      # normalize covariance
      for(i in 1:n_tissues){
        if(cor_expr[i,i] != 0){
          cor_expr[i,] <- cor_expr[i,]/sqrt(cor_expr[i,i])
          cor_expr[,i] <- cor_expr[,i]/cor_expr[i,i]
        }
      }
      # calculate tissue specificity using expression variation among tissues
      no_expr_pairs <- total_tis*(total_tis + 1)/2 - expr_tis*(expr_tis + 1)/2
      cor_expr <- cor_expr[upper.tri(cor_expr, diag = T)]
      tis_spec_score <- (sum(1 - cor_expr) + no_expr_pairs)/(total_tis*(total_tis + 1)/2)
      return(tis_spec_score)
    }
  }else{
    print("please specify method as 'tau' or 'expr_var'.")
  }
}

############################################
## 17. run_pcr
############################################
run_pcr <- function(trait, gene_expr){
  # run principal component regression
  # return a data frame of information about tissue, effect size, and p-value
  assoc_summary <- data.frame()
  
  pred_expr_eigenvalue <- eigen(t(gene_expr)%*%gene_expr)$values
  n_comp <- length(max(pred_expr_eigenvalue)/pred_expr_eigenvalue < 30)
  pred_expr_scores <- prcomp(gene_expr, center = T, scale. = T)$x[, 1:n_comp] # keep only significant components
  pcr_result <- summary(aov(trait - mean(trait) ~ pred_expr_scores))
  assoc_summary <- data.frame(
    tissue = "integrative",
    effect_size = NA,
    raw_p_value = pcr_result[[1]][["Pr(>F)"]][1])

  return(assoc_summary)
}

############################################
## 18. run_gbj
############################################

run_gbj <- function(trait, gene_expr, genotypes, eqtl_weights, use_snp_cov = TRUE, link_function = "linear"){
  # run GBJ test
  # return a data frame of information about tissue, effect size, and p-value
  if(use_snp_cov){
    # run GBJ test with tissue covriance matrix derived from SNP and eQTL covariance matrix
    stopifnot(length(genotypes) > 0, length(eqtl_weights) > 0)
    # originally designed by UTMOST for summary-level statistics
    eqtl_list <- rownames(eqtl_weights)[!apply(eqtl_weights, 1, function(x){all(x == 0)})]
    cov_snps <- cov(genotypes[, eqtl_list])
    eqtl_weights <- eqtl_weights[eqtl_list, ]
    cov_pred_tissues <- t(eqtl_weights) %*% cov_snps %*% eqtl_weights
    # normalize covariance
    for(i in 1:n_tissues){
      if(cov_pred_tissues[i,i] != 0){
        cov_pred_tissues[i,] <- cov_pred_tissues[i,]/sqrt(cov_pred_tissues[i,i])
        cov_pred_tissues[,i] <- cov_pred_tissues[,i]/cov_pred_tissues[i,i]
      }
    }
    # run GBJ test
    tmp <- calc_score_stats(null_model = glm(trait ~ 1), factor_matrix = gene_expr, link_function = link_function)
    gbj_result <- GBJ(test_stats = tmp$test_stats, cor_mat = cov_pred_tissues)
  }else{
    # run GBJ with covariance matrix derived from predicted gene expression levels
    tmp <- calc_score_stats(null_model = glm(trait ~ 1), factor_matrix = gene_expr, link_function = link_function)
    gbj_result <- GBJ(test_stats = tmp$test_stats, cor_mat = tmp$cor_mat)
  }

  assoc_summary <- data.frame(
    tissue = "integrative",
    effect_size = NA,
    raw_p_value = gbj_result$GBJ_pvalue
  )
  
  return(assoc_summary)
}

############################################
## 19. generate_ped
############################################
generate_ped <- function(genotypes, snp_info){
  # generate ped format data from allele dosages
  geno_ped <- apply(genotypes, 1, function(ind_geno){
    tmp <- ""
    for(single_snp in names(ind_geno)){
      ref_geno <- paste(rep(snp_info[snp_info$SNP == single_snp, "REF"], times = as.integer(ind_geno[single_snp])), sep = "", collapse = "")
      alt_geno <- paste(rep(snp_info[snp_info$SNP == single_snp, "ALT"], times = 2 - as.integer(ind_geno[single_snp])), sep = "", collapse = "")
      tmp <- paste0(tmp, ref_geno, alt_geno, sep = "")
    }
    tmp <- paste(unlist(strsplit(tmp, "")), collapse = " ")
    return(tmp)
  })
  geno_ped <- data.frame(fid = rownames(genotypes), iid = rownames(genotypes), pid = 0, mid = 0, sex = "unknown", phenotype = -9, geno_ped)
  return(geno_ped)
}
