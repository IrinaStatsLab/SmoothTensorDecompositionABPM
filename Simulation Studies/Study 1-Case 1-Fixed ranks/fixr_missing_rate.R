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

### Should run this on HPC 

nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

# Results generated from running SmoothHOOI algorithm on real data, with hyperparameter optimization and identifiability correction
load("./synthetic_raw.Rda") 

N <- 100 # Number of replicates 
D2 <- SecDiffMat(24) # Second-order difference matrix

# search grid
rank_grid <- as.matrix(expand.grid(r1<-c(3), r2<-c(2)))
lambda_seq <- seq(0,50,by=1)

###################################################################

# no missingness
set.seed(98765)

fixr_miss0 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  ## Generate synthetic data
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=200, noise_level=1, pattern="random", percent=0)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)

  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)

  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)

  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)

  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)

  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)

  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)

  fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data,
                    L = fpca_res$Lfpca, true_L = L_tilde)

  ## CP3 with imputation and smoothness and orthogonality
  cp3_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_loss <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                   L = cp3_res$Lcp, true_L = L_tilde)
  
  cp3_loss1 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_loss2 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_loss3 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_top3 <- cp6_res$Lcp_ordered[,1:3]
  
  cp6_loss <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                   L = Lcp6_top3, true_L = L_tilde)
  
  cp6_loss1 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_loss2 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_loss3 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_no_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_no_res <- qr.Q(qr(cp3_no_res$Lcp))
  
  cp3_no_loss <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp3_no_res, true_L = L_tilde)
  
  cp3_no_loss1 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_no_loss2 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_no_loss3 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_no_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_no_top3 <- cp6_no_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_no_res <- qr.Q(qr(Lcp6_no_top3))
  
  cp6_no_loss <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp6_no_res, true_L = L_tilde)
  
  cp6_no_loss1 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_no_loss2 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_no_loss3 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data,
                     L = mfpca_res$Lmfpca, true_L = L_tilde)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss,
                 "cp3_loss"=cp3_loss, "cp3_loss1"=cp3_loss1, "cp3_loss2"=cp3_loss2, "cp3_loss3"=cp3_loss3, 
                 "cp6_loss"=cp6_loss, "cp6_loss1"=cp6_loss1,"cp6_loss2"=cp6_loss2, "cp6_loss3"=cp6_loss3,
                 "cp3_no_loss"=cp3_no_loss, "cp3_no_loss1"=cp3_no_loss1, "cp3_no_loss2"=cp3_no_loss2, "cp3_no_loss3"=cp3_no_loss3, 
                 "cp6_no_loss"=cp6_no_loss, "cp6_no_loss1"=cp6_no_loss1, "cp6_no_loss2"=cp6_no_loss2,"cp6_no_loss3"=cp6_no_loss3,
                 "mfpca_loss"=mfpca_loss)
  output
}

save(fixr_miss0, file="./fixr_new/fixr_miss0.Rda")

## missing rate 10% 
set.seed(97531)

fixr_miss10 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                       E = E, p=200, noise_level=1, pattern="random", percent=0.1)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data,
                    L = fpca_res$Lfpca, true_L = L_tilde)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_loss <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                   L = cp3_res$Lcp, true_L = L_tilde)
  
  cp3_loss1 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_loss2 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_loss3 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_top3 <- cp6_res$Lcp_ordered[,1:3]
  
  cp6_loss <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                   L = Lcp6_top3, true_L = L_tilde)
  
  cp6_loss1 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_loss2 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_loss3 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_no_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_no_res <- qr.Q(qr(cp3_no_res$Lcp))
  
  cp3_no_loss <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp3_no_res, true_L = L_tilde)
  
  cp3_no_loss1 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_no_loss2 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_no_loss3 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_no_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_no_top3 <- cp6_no_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_no_res <- qr.Q(qr(Lcp6_no_top3))
  
  cp6_no_loss <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp6_no_res, true_L = L_tilde)
  
  cp6_no_loss1 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_no_loss2 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_no_loss3 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data,
                     L = mfpca_res$Lmfpca, true_L = L_tilde)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss,
                 "cp3_loss"=cp3_loss, "cp3_loss1"=cp3_loss1, "cp3_loss2"=cp3_loss2, "cp3_loss3"=cp3_loss3, 
                 "cp6_loss"=cp6_loss, "cp6_loss1"=cp6_loss1,"cp6_loss2"=cp6_loss2, "cp6_loss3"=cp6_loss3,
                 "cp3_no_loss"=cp3_no_loss, "cp3_no_loss1"=cp3_no_loss1, "cp3_no_loss2"=cp3_no_loss2, "cp3_no_loss3"=cp3_no_loss3, 
                 "cp6_no_loss"=cp6_no_loss, "cp6_no_loss1"=cp6_no_loss1, "cp6_no_loss2"=cp6_no_loss2,"cp6_no_loss3"=cp6_no_loss3,
                 "mfpca_loss"=mfpca_loss)
  output
}

save(fixr_miss10, file="./fixr_new/fixr_miss10.Rda")

## missing rate 20% 
set.seed(86420)

fixr_miss20 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=200, noise_level=1, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data,
                    L = fpca_res$Lfpca, true_L = L_tilde)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_loss <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                   L = cp3_res$Lcp, true_L = L_tilde)
  
  cp3_loss1 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_loss2 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_loss3 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_top3 <- cp6_res$Lcp_ordered[,1:3]
  
  cp6_loss <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                   L = Lcp6_top3, true_L = L_tilde)
  
  cp6_loss1 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_loss2 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_loss3 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_no_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_no_res <- qr.Q(qr(cp3_no_res$Lcp))
  
  cp3_no_loss <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp3_no_res, true_L = L_tilde)
  
  cp3_no_loss1 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_no_loss2 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_no_loss3 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_no_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_no_top3 <- cp6_no_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_no_res <- qr.Q(qr(Lcp6_no_top3))
  
  cp6_no_loss <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp6_no_res, true_L = L_tilde)
  
  cp6_no_loss1 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_no_loss2 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_no_loss3 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data,
                     L = mfpca_res$Lmfpca, true_L = L_tilde)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss,
                 "cp3_loss"=cp3_loss, "cp3_loss1"=cp3_loss1, "cp3_loss2"=cp3_loss2, "cp3_loss3"=cp3_loss3, 
                 "cp6_loss"=cp6_loss, "cp6_loss1"=cp6_loss1,"cp6_loss2"=cp6_loss2, "cp6_loss3"=cp6_loss3,
                 "cp3_no_loss"=cp3_no_loss, "cp3_no_loss1"=cp3_no_loss1, "cp3_no_loss2"=cp3_no_loss2, "cp3_no_loss3"=cp3_no_loss3, 
                 "cp6_no_loss"=cp6_no_loss, "cp6_no_loss1"=cp6_no_loss1, "cp6_no_loss2"=cp6_no_loss2,"cp6_no_loss3"=cp6_no_loss3,
                 "mfpca_loss"=mfpca_loss,
                 "Lcp3_ortho"=cp3_res$Lcp, "Lcp6_ortho"=cp6_res$Lcp_ordered, "Lcp3"=cp3_no_res$Lcp, "Lcp6"=cp6_no_res$Lcp, 
                 "Lcp3_ordered"=cp3_no_res$Lcp_ordered, "Lcp6_ordered"=cp6_no_res$Lcp_ordered)
  output
}

save(fixr_miss20, file="./fixr_new/fixr_miss20.Rda")

## missing rate 50% 
set.seed(15151)

fixr_miss50 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
            E = E, p=200, noise_level=1, pattern="random", percent=0.5)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, smooth_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, smooth_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss <- loss(tnsr = fpca_res$est, smooth_tnsr=Msmooth@data,
                    L = fpca_res$Lfpca, true_L = L_tilde)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_loss <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                   L = cp3_res$Lcp, true_L = L_tilde)
  
  cp3_loss1 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_loss2 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_loss3 <- loss(tnsr = cp3_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp3_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_top3 <- cp6_res$Lcp_ordered[,1:3]
  
  cp6_loss <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                   L = Lcp6_top3, true_L = L_tilde)
  
  cp6_loss1 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_loss2 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_loss3 <- loss(tnsr = cp6_res$est, smooth_tnsr=Msmooth@data,
                    L =  as.matrix(cp6_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_no_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_no_res <- qr.Q(qr(cp3_no_res$Lcp))
  
  cp3_no_loss <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp3_no_res, true_L = L_tilde)
  
  cp3_no_loss1 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_no_loss2 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_no_loss3 <- loss(tnsr = cp3_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp3_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_no_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_no_top3 <- cp6_no_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_no_res <- qr.Q(qr(Lcp6_no_top3))
  
  cp6_no_loss <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                      L = Lcp_cp6_no_res, true_L = L_tilde)
  
  cp6_no_loss1 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_no_loss2 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_no_loss3 <- loss(tnsr = cp6_no_res$est, smooth_tnsr=Msmooth@data,
                       L =  as.matrix(cp6_no_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss <- loss(tnsr = mfpca_res$est, smooth_tnsr=Msmooth@data,
                     L = mfpca_res$Lmfpca, true_L = L_tilde)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, "fpca_loss"=fpca_loss,
                 "cp3_loss"=cp3_loss, "cp3_loss1"=cp3_loss1, "cp3_loss2"=cp3_loss2, "cp3_loss3"=cp3_loss3, 
                 "cp6_loss"=cp6_loss, "cp6_loss1"=cp6_loss1,"cp6_loss2"=cp6_loss2, "cp6_loss3"=cp6_loss3,
                 "cp3_no_loss"=cp3_no_loss, "cp3_no_loss1"=cp3_no_loss1, "cp3_no_loss2"=cp3_no_loss2, "cp3_no_loss3"=cp3_no_loss3, 
                 "cp6_no_loss"=cp6_no_loss, "cp6_no_loss1"=cp6_no_loss1, "cp6_no_loss2"=cp6_no_loss2,"cp6_no_loss3"=cp6_no_loss3,
                 "mfpca_loss"=mfpca_loss)
  output
}

save(fixr_miss50, file="./fixr_new/fixr_miss50.Rda")

stopCluster(cl)
