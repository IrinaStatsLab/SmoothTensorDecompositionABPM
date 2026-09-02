library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(doRNG)
library(rTensor)
library(MASS)
library(refund)
library(SmoothHOOI)
library(MFPCA)

source("./fpca_run.R") # script for running univariate FPCA
source("./cp_run.R")   # script for running CP decomposition with imputation and smoothness
source("./mfpca_run.R") # script for running MFPCA
source("./parafac2_run.R") # script for running PARAFAC2 (SCA-PF2 and SCA-IND)
source("./sca_ecp_run.R") # script for running SCA-ECP

## Should run this on HPC 

nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

load("./synthetic_raw.Rda") 

N <- 100
percent <- 0.2
noise_level <- 1
D2 <- SecDiffMat(24)
rank_grid <- as.matrix(expand.grid(r1<-seq(2,6,by=1), r2<-c(2,3)))
lambda_seq <- seq(0,50,by=1)

## p = 50
set.seed(345678)

full_p50 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=50, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data)
  
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data)
  
  ## FPCA, allowing the number of principal components to vary
  fpca_cf_res <- fpca_run(Mmiss, center=FALSE)
  
  fpca_cf_loss <- loss(tnsr = fpca_cf_res$est, true_tnsr=Msmooth@data)
  
  ## FPCA, allowing the number of principal components to vary
  fpca_ct_res <- fpca_run(Mmiss, center=TRUE)
  
  fpca_ct_loss <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data)
  
  ## CP with imputation and smoothness and orthogonality
  cp_ortsmo_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  cp_ortsmo_loss <- loss(tnsr = cp_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP with imputation and smoothness, no orthogonality 
  cp_smooth_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))
  
  cp_smooth_loss <- loss(tnsr = cp_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## CP2 with imputation and smoothness and orthogonality
  cp2_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(2), const=c("ortsmo", "uncons", "uncons"))
  
  cp2_ortsmo_loss <- loss(tnsr = cp2_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP2 with imputation and smoothness, no orthogonality 
  cp2_smooth_res <- cp_run(Mmiss@data, R_seq=c(2), const=c("smooth", "uncons", "uncons"))
  
  cp2_smooth_loss <- loss(tnsr = cp2_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP6 with imputation and smoothness, no orthogonality 
  cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=seq(2,6,by=1))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA-2
  mfpca2_res <- mfpca_run(Mmiss, M_seq=c(2))
  
  mfpca2_loss <- loss(tnsr = mfpca2_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA-3
  mfpca3_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca3_loss <- loss(tnsr = mfpca3_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA-6
  mfpca6_res <- mfpca_run(Mmiss, M_seq=c(6))
  
  mfpca6_loss <- loss(tnsr = mfpca6_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-flexible with imputation and smoothness, no orthogonality 
  pf_smooth_res <- parafac2_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))
  
  pf_smooth_loss <- loss(tnsr = pf_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-2 with imputation and smoothness, no orthogonality 
  pf2_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(2), const=c("smooth", "uncons", "uncons"))
  
  pf2_smooth_loss <- loss(tnsr = pf2_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
  pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  pf3_smooth_loss <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-flexible with imputation and smoothness and orthogonality 
  pf_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  pf_ortsmo_loss <- loss(tnsr = pf_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-2 with imputation and smoothness and orthogonality
  pf2_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(2), const=c("ortsmo", "uncons", "uncons"))
  
  pf2_ortsmo_loss <- loss(tnsr = pf2_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-3 with imputation and smoothness and orthogonality
  pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  pf3_ortsmo_loss <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## SCA-ECP-2 with imputation and smoothness and orthogonality
  ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  ecp_ortsmo_loss <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## SCA-ECP-2 with imputation and smoothness and orthogonality
  ecp2_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(2), const=c("ortsmo", "uncons", "uncons"))
  
  ecp2_ortsmo_loss <- loss(tnsr = ecp2_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## SCA-ECP with imputation and smoothness and orthogonality
  ecp3_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  ecp3_ortsmo_loss <- loss(tnsr = ecp3_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  output <- list("oracle_rank"=oracle_opt$opt_para[1],"oracle_r2"=oracle_opt$opt_para[2], "oracle_loss"=oracle_loss,
                 "kcv_rank"=kcv_opt$opt_para[1],"kcv_r2"=kcv_opt$opt_para[2], "kcv_loss"=kcv_loss, 
                 "fpca_cf_rank"=fpca_cf_res$min_col, "fpca_cf_loss"=fpca_cf_loss, "fpca_ct_rank"=fpca_ct_res$min_col, "fpca_ct_loss"=fpca_ct_loss,
                 "cp_ortsmo_rank"=cp_ortsmo_res$best_R, "cp_ortsmo_loss"=cp_ortsmo_loss, 
                 "cp2_ortsmo_loss"=cp2_ortsmo_loss,"cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp6_ortsmo_loss"=cp6_ortsmo_loss,
                 "cp_smooth_rank"=cp_smooth_res$best_R, "cp_smooth_loss"=cp_smooth_loss, 
                 "cp2_smooth_loss"=cp2_smooth_loss,"cp3_smooth_loss"=cp3_smooth_loss, "cp6_smooth_loss"=cp6_smooth_loss,
                 "pf_smooth_loss"=pf_smooth_loss, "pf_smooth_rank"=pf_smooth_res$best_R,
                 "pf2_smooth_loss"=pf2_smooth_loss, "pf3_smooth_loss"=pf3_smooth_loss,
                 "pf_ortsmo_loss"=pf_ortsmo_loss, "pf_ortsmo_rank"=pf_ortsmo_res$best_R,
                 "pf2_ortsmo_loss"=pf2_ortsmo_loss, "pf3_ortsmo_loss"=pf3_ortsmo_loss,
                 "ecp_ortsmo_loss"=ecp_ortsmo_loss, "ecp_ortsmo_rank"=ecp_ortsmo_res$best_R,
                 "ecp2_ortsmo_loss"=ecp2_ortsmo_loss, "ecp3_ortsmo_loss"=ecp3_ortsmo_loss,
                 "mfpca_rank"=mfpca_res$best_M, "mfpca_true_rank"=mfpca_res$true_M, "mfpca_loss"=mfpca_loss, 
                 "mfpca2_loss"=mfpca2_loss, "mfpca2_true_rank"=mfpca2_res$true_M,
                 "mfpca3_loss"=mfpca3_loss, "mfpca3_true_rank"=mfpca3_res$true_M,
                 "mfpca6_loss"=mfpca6_loss, "mfpca6_true_rank"=mfpca6_res$true_M)
  
  output
}

save(full_p50, file="./full_new/full_p50.Rda")

## p = 500
set.seed(567890)

full_p500 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=500, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  oracle_opt <- oracle_memeff(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                              init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data)
  
  kcv_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                        L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data)
  
  ## FPCA, allowing the number of principal components to vary
  fpca_cf_res <- fpca_run(Mmiss, center=FALSE)
  
  fpca_cf_loss <- loss(tnsr = fpca_cf_res$est, true_tnsr=Msmooth@data)
  
  ## FPCA, allowing the number of principal components to vary
  fpca_ct_res <- fpca_run(Mmiss, center=TRUE)
  
  fpca_ct_loss <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data)
  
  ## CP with imputation and smoothness and orthogonality
  cp_ortsmo_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  cp_ortsmo_loss <- loss(tnsr = cp_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP with imputation and smoothness, no orthogonality 
  cp_smooth_res <- cp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))
  
  cp_smooth_loss <- loss(tnsr = cp_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## CP2 with imputation and smoothness and orthogonality
  cp2_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(2), const=c("ortsmo", "uncons", "uncons"))
  
  cp2_ortsmo_loss <- loss(tnsr = cp2_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP2 with imputation and smoothness, no orthogonality 
  cp2_smooth_res <- cp_run(Mmiss@data, R_seq=c(2), const=c("smooth", "uncons", "uncons"))
  
  cp2_smooth_loss <- loss(tnsr = cp2_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## CP6 with imputation and smoothness, no orthogonality 
  cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=seq(2,6,by=1))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA-2
  mfpca2_res <- mfpca_run(Mmiss, M_seq=c(2))
  
  mfpca2_loss <- loss(tnsr = mfpca2_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA-3
  mfpca3_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca3_loss <- loss(tnsr = mfpca3_res$est, true_tnsr=Msmooth@data)
  
  ## MFPCA-6
  mfpca6_res <- mfpca_run(Mmiss, M_seq=c(6))
  
  mfpca6_loss <- loss(tnsr = mfpca6_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-flexible with imputation and smoothness, no orthogonality 
  pf_smooth_res <- parafac2_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("smooth", "uncons", "uncons"))
  
  pf_smooth_loss <- loss(tnsr = pf_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-2 with imputation and smoothness, no orthogonality 
  pf2_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(2), const=c("smooth", "uncons", "uncons"))
  
  pf2_smooth_loss <- loss(tnsr = pf2_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
  pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  pf3_smooth_loss <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-flexible with imputation and smoothness and orthogonality 
  pf_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  pf_ortsmo_loss <- loss(tnsr = pf_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-2 with imputation and smoothness and orthogonality
  pf2_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(2), const=c("ortsmo", "uncons", "uncons"))
  
  pf2_ortsmo_loss <- loss(tnsr = pf2_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## PARAFAC2-3 with imputation and smoothness and orthogonality
  pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  pf3_ortsmo_loss <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## SCA-ECP-2 with imputation and smoothness and orthogonality
  ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=seq(2,6,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  ecp_ortsmo_loss <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## SCA-ECP-2 with imputation and smoothness and orthogonality
  ecp2_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(2), const=c("ortsmo", "uncons", "uncons"))
  
  ecp2_ortsmo_loss <- loss(tnsr = ecp2_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  ## SCA-ECP with imputation and smoothness and orthogonality
  ecp3_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  ecp3_ortsmo_loss <- loss(tnsr = ecp3_ortsmo_res$est, true_tnsr=Msmooth@data)
  
  output <- list("oracle_rank"=oracle_opt$opt_para[1],"oracle_r2"=oracle_opt$opt_para[2], "oracle_loss"=oracle_loss,
                 "kcv_rank"=kcv_opt$opt_para[1],"kcv_r2"=kcv_opt$opt_para[2], "kcv_loss"=kcv_loss, 
                 "fpca_cf_rank"=fpca_cf_res$min_col, "fpca_cf_loss"=fpca_cf_loss, "fpca_ct_rank"=fpca_ct_res$min_col, "fpca_ct_loss"=fpca_ct_loss,
                 "cp_ortsmo_rank"=cp_ortsmo_res$best_R, "cp_ortsmo_loss"=cp_ortsmo_loss, 
                 "cp2_ortsmo_loss"=cp2_ortsmo_loss,"cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp6_ortsmo_loss"=cp6_ortsmo_loss,
                 "cp_smooth_rank"=cp_smooth_res$best_R, "cp_smooth_loss"=cp_smooth_loss, 
                 "cp2_smooth_loss"=cp2_smooth_loss,"cp3_smooth_loss"=cp3_smooth_loss, "cp6_smooth_loss"=cp6_smooth_loss,
                 "pf_smooth_loss"=pf_smooth_loss, "pf_smooth_rank"=pf_smooth_res$best_R,
                 "pf2_smooth_loss"=pf2_smooth_loss, "pf3_smooth_loss"=pf3_smooth_loss,
                 "pf_ortsmo_loss"=pf_ortsmo_loss, "pf_ortsmo_rank"=pf_ortsmo_res$best_R,
                 "pf2_ortsmo_loss"=pf2_ortsmo_loss, "pf3_ortsmo_loss"=pf3_ortsmo_loss,
                 "ecp_ortsmo_loss"=ecp_ortsmo_loss, "ecp_ortsmo_rank"=ecp_ortsmo_res$best_R,
                 "ecp2_ortsmo_loss"=ecp2_ortsmo_loss, "ecp3_ortsmo_loss"=ecp3_ortsmo_loss,
                 "mfpca_rank"=mfpca_res$best_M, "mfpca_true_rank"=mfpca_res$true_M, "mfpca_loss"=mfpca_loss, 
                 "mfpca2_loss"=mfpca2_loss, "mfpca2_true_rank"=mfpca2_res$true_M,
                 "mfpca3_loss"=mfpca3_loss, "mfpca3_true_rank"=mfpca3_res$true_M,
                 "mfpca6_loss"=mfpca6_loss, "mfpca6_true_rank"=mfpca6_res$true_M)
  
  output
}

save(full_p500, file="./full_new/full_p500.Rda")

stopCluster(cl)
