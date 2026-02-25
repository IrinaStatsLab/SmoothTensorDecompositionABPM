library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(doRNG)
library(rTensor)
library(MASS)
library(refund)
library(SmoothHOOI)

source("./fpca_run.R") # script for running univariate FPCA
source("./cp_run.R")   # script for running CP decomposition with imputation and smoothness
source("./mfpca_run.R") # script for running MFPCA

## Should run this on HPC 

nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)


load("./synthetic_raw.Rda") 

N <- 100

D2 <- SecDiffMat(366)

Ldense_1 <- spline(L_tilde[,1], n = 366)$y
Ldense_2 <- spline(L_tilde[,2], n = 366)$y
Ldense_3 <- spline(L_tilde[,3], n = 366)$y
Ldense <- cbind(Ldense_1, Ldense_2, Ldense_3)
Ldense_init <- cbind(Ldense_1, Ldense_2, Ldense_3)

Ldense_svd <- svd(Ldense_init)
Ldense <- Ldense_svd$u %*% t(Ldense_svd$v)


error_sd <- stats::sd(E, na.rm=T)

rank_grid <- as.matrix(expand.grid(r1<-seq(2,4,by=1), r2<-seq(2,4,by=1)))
lambda_seq <- seq(0,500,by=50)

set.seed(36612111)

ar1_large_full <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data4(L=Ldense, b=20, r2=3, p=30, pattern="random", percent=0.2, phi=0.5, ar1_noise_sd = 0.01)
  
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle_memeff(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                              init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv_memeff(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                        L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data)
  
  
  ## 10-fold h cross-validation
  kcv_hblock_opt <- kcv_hblock(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10, h=3,
                               L0=NULL, D=D2, tol=0.1, max_iter=500, init=0)
  
  kcv_hblock_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_hblock_opt$opt_para[1:2]), lambda=as.numeric(kcv_hblock_opt$opt_para[3]),
                           L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_hblock_tilde <- MakeIdent(L=kcv_hblock_res$L, G=kcv_hblock_res$G, R=kcv_hblock_res$R)
  
  kcv_hblock_loss <- loss(tnsr = kcv_hblock_res$est, smooth_tnsr=Msmooth@data)
  

  ## FPCA, allowing the number of principal components to vary
  fpca_res <- fpca_run(Mmiss)
  
  fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data)
  
  ## CP with imputation and smoothness and orthogonality
  cp_res <- cp_run(Mmiss@data, R_seq=seq(3,10,by=1), const=c("ortsmo", "uncons", "uncons"))
  
  cp_loss <- loss(tnsr = cp_res$est, smooth_tnsr=Msmooth@data)
  
  ## CP with imputation and smoothness, no orthogonality
  cp_no_res <- cp_run(Mmiss@data, R_seq=seq(3,10,by=1), const=c("smooth", "uncons", "uncons"))
  
  cp_no_loss <- loss(tnsr = cp_no_res$est, smooth_tnsr=Msmooth@data)
  
  ## MFPCA
  mfpca_res <- mfpca_run(Mmiss, M_seq=seq(3,10,by=1))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data)
  
  output <- list("oracle_loss"=oracle_loss$loss_M, "oracle_para"=oracle_opt$opt_para,
                 "kcv_loss"=kcv_loss$loss_M, "kcv_para"=kcv_opt$opt_para,
                 "kcv_hblock_loss"=kcv_hblock_loss$loss_M, "kcv_hblock_para"=kcv_hblock_opt$opt_para,
                 "fpca_loss"=fpca_loss$loss_M, "fpca_rank"=ncol(fpca_res$Lfpca),
                 "cp_loss"=cp_loss$loss_M,"cp_rank"=cp_res$best_R,
                 "cp_no_loss"=cp_no_loss$loss_M, "cp_no_rank"=cp_no_res$best_R,
                 "mfpca_loss"=mfpca_loss$loss_M, "mfpca_rank"=ncol(mfpca_res$Lmfpca)
  )
  output
}

save(ar1_large_full, file="./ar1_large_full.Rda")


stopCluster(cl)






