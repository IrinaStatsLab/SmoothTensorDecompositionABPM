library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(doRNG)
library(rTensor)
library(MASS)
library(refund)
library(SmoothHOOI)
library(MFPCA)
library(multiway)

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
D2 <- SecDiffMat(24)

rank_grid <- as.matrix(expand.grid(r1<-c(3), r2<-c(2)))
lambda_seq <- seq(0,50,by=1)

## noise level 0*empirical var
set.seed(11111)

fixr_noise0 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=200, noise_level=0, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss_list <- list() 
  for (k in 1:length(fpca_res$Lfpca_list)){
    Lfpca_k <- fpca_res$Lfpca_list[[k]]
    fpca_loss_list[[k]] <- loss(tnsr = fpca_res$est, true_tnsr=Msmooth@data,
                                L = Lfpca_k, true_L = L_tilde)
  }
  fpca_loss_unlist <- sapply(fpca_loss_list, unlist)            
  fpca_loss_avg <- rowMeans(fpca_loss_unlist, na.rm = TRUE)
  
  fpca_loss <- as.list(fpca_loss_avg)
  
  ## FPCA, fix number of principal components to 3
  fpca_ct_res <- fpca_run(Mmiss, npc=3, center=TRUE)
  
  fpca_ct_loss_list <- list() 
  for (k in 1:length(fpca_ct_res$Lfpca_list)){
    Lfpca_ct_k <- fpca_ct_res$Lfpca_list[[k]]
    fpca_ct_loss_list[[k]] <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data,
                                   L = Lfpca_ct_k, true_L = L_tilde)
  }
  fpca_ct_loss_unlist <- sapply(fpca_ct_loss_list, unlist)            
  fpca_ct_loss_avg <- rowMeans(fpca_ct_loss_unlist, na.rm = TRUE)
  
  fpca_ct_loss <- as.list(fpca_ct_loss_avg)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = cp3_ortsmo_res$Lcp, true_L = L_tilde)
  
  cp3_ortsmo_loss1 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_ortsmo_loss2 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_ortsmo_loss3 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_ortsmo_top3 <- cp6_ortsmo_res$Lcp_ordered[,1:3]
  
  cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = Lcp6_ortsmo_top3, true_L = L_tilde)
  
  cp6_ortsmo_loss1 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_ortsmo_loss2 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_ortsmo_loss3 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_smooth_res <- qr.Q(qr(cp3_smooth_res$Lcp))
  
  cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp3_smooth_res, true_L = L_tilde)
  
  cp3_smooth_loss1 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_smooth_loss2 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_smooth_loss3 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_smooth_top3 <- cp6_smooth_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_smooth_res <- qr.Q(qr(Lcp6_smooth_top3))
  
  cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp6_smooth_res, true_L = L_tilde)
  
  cp6_smooth_loss1 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_smooth_loss2 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_smooth_loss3 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss_list <- list() 
  for (k in 1:length(mfpca_res$Lmfpca_list)){
    Lmfpca_k <- mfpca_res$Lmfpca_list[[k]]
    mfpca_loss_list[[k]] <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data,
                                 L = Lmfpca_k, true_L = L_tilde)
  }
  mfpca_loss_unlist <- sapply(mfpca_loss_list, unlist)            
  mfpca_loss_avg <- rowMeans(mfpca_loss_unlist, na.rm = TRUE)
  
  mfpca_loss <- as.list(mfpca_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
  pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  pf3_smooth_loss_list <- list() 
  for (k in 1:length(pf3_smooth_res$Lparafac2_list)){
    Lpf3_smooth_k <- qr.Q(qr(pf3_smooth_res$Lparafac2_list[[k]]))
    pf3_smooth_loss_list[[k]] <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_smooth_k, true_L = L_tilde)
  }
  pf3_smooth_loss_unlist <- sapply(pf3_smooth_loss_list, unlist)            
  pf3_smooth_loss_avg <- rowMeans(pf3_smooth_loss_unlist, na.rm = TRUE)
  
  pf3_smooth_loss <- as.list(pf3_smooth_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness and orthogonality
  pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  pf3_ortsmo_loss_list <- list() 
  for (k in 1:length(pf3_ortsmo_res$Lparafac2_list)){
    Lpf3_ortsmo_k <- pf3_ortsmo_res$Lparafac2_list[[k]]
    pf3_ortsmo_loss_list[[k]] <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_ortsmo_k, true_L = L_tilde)
  }
  pf3_ortsmo_loss_unlist <- sapply(pf3_ortsmo_loss_list, unlist)            
  pf3_ortsmo_loss_avg <- rowMeans(pf3_ortsmo_loss_unlist, na.rm = TRUE)
  
  pf3_ortsmo_loss <- as.list(pf3_ortsmo_loss_avg)
  
  
  ## SCA-ECP with imputation and smoothness and orthogonality
  ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  ecp_ortsmo_loss_list <- list() 
  for (k in 1:length(ecp_ortsmo_res$Lecp_list)){
    Lecp_ortsmo_k <- ecp_ortsmo_res$Lecp_list[[k]]
    ecp_ortsmo_loss_list[[k]] <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lecp_ortsmo_k, true_L = L_tilde)
  }
  ecp_ortsmo_loss_unlist <- sapply(ecp_ortsmo_loss_list, unlist)            
  ecp_ortsmo_loss_avg <- rowMeans(ecp_ortsmo_loss_unlist, na.rm = TRUE)
  
  ecp_ortsmo_loss <- as.list(ecp_ortsmo_loss_avg)
  
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, 
                 "fpca_loss"=fpca_loss, "fpca_loss_unlist"=fpca_loss_unlist,"fpca_ct_loss"=fpca_ct_loss, "fpca_ct_loss_unlist"=fpca_ct_loss_unlist,
                 "cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp3_ortsmo_loss1"=cp3_ortsmo_loss1, "cp3_ortsmo_loss2"=cp3_ortsmo_loss2, "cp3_ortsmo_loss3"=cp3_ortsmo_loss3, 
                 "cp6_ortsmo_loss"=cp6_ortsmo_loss, "cp6_ortsmo_loss1"=cp6_ortsmo_loss1,"cp6_ortsmo_loss2"=cp6_ortsmo_loss2, "cp6_ortsmo_loss3"=cp6_ortsmo_loss3,
                 "cp3_smooth_loss"=cp3_smooth_loss, "cp3_smooth_loss1"=cp3_smooth_loss1, "cp3_smooth_loss2"=cp3_smooth_loss2, "cp3_smooth_loss3"=cp3_smooth_loss3, 
                 "cp6_smooth_loss"=cp6_smooth_loss, "cp6_smooth_loss1"=cp6_smooth_loss1, "cp6_smooth_loss2"=cp6_smooth_loss2,"cp6_smooth_loss3"=cp6_smooth_loss3,
                 "mfpca_loss"=mfpca_loss,
                 "pf3_smooth_loss"=pf3_smooth_loss,
                 "pf3_ortsmo_loss"=pf3_ortsmo_loss,
                 "ecp_ortsmo_loss"=ecp_ortsmo_loss)
  output
}

save(fixr_noise0, file="./fixr_new/fixr_noise0.Rda")

## noise level 0.1*empirical var
# set.seed(76543)
# 
# fixr_noise01 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
#   sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
#                         E = E, p=200, noise_level=0.1, pattern="random", percent=0.2)
#   Mmiss <- sim_data$sim_Mmiss
#   Msmooth <- sim_data$sim_Msmooth
#   
#   ## Oracle
#   oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
#                        init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
#   
#   oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
#                        L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
#   
#   oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
#   
#   oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data,
#                       L = oracle_tilde$L_tilde, true_L = L_tilde)
#   
#   ## 10-fold cross-validation
#   kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
#                  L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
#   
#   kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
#                     L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
#   
#   kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
#   
#   kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data,
#                    L = kcv_tilde$L_tilde, true_L = L_tilde)
#   
#   ## FPCA, fix number of principal components to 3
#   fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
#   
#   fpca_loss_list <- list() 
#   for (k in 1:length(fpca_res$Lfpca_list)){
#     Lfpca_k <- fpca_res$Lfpca_list[[k]]
#     fpca_loss_list[[k]] <- loss(tnsr = fpca_res$est, true_tnsr=Msmooth@data,
#                                 L = Lfpca_k, true_L = L_tilde)
#   }
#   fpca_loss_unlist <- sapply(fpca_loss_list, unlist)            
#   fpca_loss_avg <- rowMeans(fpca_loss_unlist, na.rm = TRUE)
#   
#   fpca_loss <- as.list(fpca_loss_avg)
#   
#   ## FPCA, fix number of principal components to 3
#   fpca_ct_res <- fpca_run(Mmiss, npc=3, center=TRUE)
#   
#   fpca_ct_loss_list <- list() 
#   for (k in 1:length(fpca_ct_res$Lfpca_list)){
#     Lfpca_ct_k <- fpca_ct_res$Lfpca_list[[k]]
#     fpca_ct_loss_list[[k]] <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data,
#                                    L = Lfpca_ct_k, true_L = L_tilde)
#   }
#   fpca_ct_loss_unlist <- sapply(fpca_ct_loss_list, unlist)            
#   fpca_ct_loss_avg <- rowMeans(fpca_ct_loss_unlist, na.rm = TRUE)
#   
#   fpca_ct_loss <- as.list(fpca_ct_loss_avg)
#   
#   ## CP3 with imputation and smoothness and orthogonality
#   cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
#   
#   cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
#                           L = cp3_ortsmo_res$Lcp, true_L = L_tilde)
#   
#   cp3_ortsmo_loss1 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
#   
#   cp3_ortsmo_loss2 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
#   
#   cp3_ortsmo_loss3 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
#   
#   ## CP6 with imputation and smoothness and orthogonality
#   cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
#   
#   Lcp6_ortsmo_top3 <- cp6_ortsmo_res$Lcp_ordered[,1:3]
#   
#   cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
#                           L = Lcp6_ortsmo_top3, true_L = L_tilde)
#   
#   cp6_ortsmo_loss1 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
#   
#   cp6_ortsmo_loss2 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
#   
#   cp6_ortsmo_loss3 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
#   
#   ## CP3 with imputation and smoothness, no orthogonality 
#   cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
#   
#   Lcp_cp3_smooth_res <- qr.Q(qr(cp3_smooth_res$Lcp))
#   
#   cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
#                           L = Lcp_cp3_smooth_res, true_L = L_tilde)
#   
#   cp3_smooth_loss1 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp3_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
#   
#   cp3_smooth_loss2 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp3_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
#   
#   cp3_smooth_loss3 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp3_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
#   
#   ## CP6 with imputation and smoothness, no orthogonality
#   cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
#   
#   Lcp6_smooth_top3 <- cp6_smooth_res$Lcp_ordered[,1:3]
#   
#   Lcp_cp6_smooth_res <- qr.Q(qr(Lcp6_smooth_top3))
#   
#   cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
#                           L = Lcp_cp6_smooth_res, true_L = L_tilde)
#   
#   cp6_smooth_loss1 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp6_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
#   
#   cp6_smooth_loss2 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp6_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
#   
#   cp6_smooth_loss3 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
#                            L =  as.matrix(cp6_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
#   
#   
#   ## MFPCA 
#   mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
#   
#   mfpca_loss_list <- list() 
#   for (k in 1:length(mfpca_res$Lmfpca_list)){
#     Lmfpca_k <- mfpca_res$Lmfpca_list[[k]]
#     mfpca_loss_list[[k]] <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data,
#                                  L = Lmfpca_k, true_L = L_tilde)
#   }
#   mfpca_loss_unlist <- sapply(mfpca_loss_list, unlist)            
#   mfpca_loss_avg <- rowMeans(mfpca_loss_unlist, na.rm = TRUE)
#   
#   mfpca_loss <- as.list(mfpca_loss_avg)
#   
#   
#   ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
#   pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
#   
#   pf3_smooth_loss_list <- list() 
#   for (k in 1:length(pf3_smooth_res$Lparafac2_list)){
#     Lpf3_smooth_k <- qr.Q(qr(pf3_smooth_res$Lparafac2_list[[k]]))
#     pf3_smooth_loss_list[[k]] <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data,
#                                       L = Lpf3_smooth_k, true_L = L_tilde)
#   }
#   pf3_smooth_loss_unlist <- sapply(pf3_smooth_loss_list, unlist)            
#   pf3_smooth_loss_avg <- rowMeans(pf3_smooth_loss_unlist, na.rm = TRUE)
#   
#   pf3_smooth_loss <- as.list(pf3_smooth_loss_avg)
#   
# 
#   ## PARAFAC2-3 with imputation and smoothness and orthogonality
#   pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
#   
#   pf3_ortsmo_loss_list <- list() 
#   for (k in 1:length(pf3_ortsmo_res$Lparafac2_list)){
#     Lpf3_ortsmo_k <- pf3_ortsmo_res$Lparafac2_list[[k]]
#     pf3_ortsmo_loss_list[[k]] <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data,
#                                       L = Lpf3_ortsmo_k, true_L = L_tilde)
#   }
#   pf3_ortsmo_loss_unlist <- sapply(pf3_ortsmo_loss_list, unlist)            
#   pf3_ortsmo_loss_avg <- rowMeans(pf3_ortsmo_loss_unlist, na.rm = TRUE)
#   
#   pf3_ortsmo_loss <- as.list(pf3_ortsmo_loss_avg)
#   
#   
#   ## SCA-ECP with imputation and smoothness and orthogonality
#   ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
#   
#   ecp_ortsmo_loss_list <- list() 
#   for (k in 1:length(ecp_ortsmo_res$Lecp_list)){
#     Lecp_ortsmo_k <- ecp_ortsmo_res$Lecp_list[[k]]
#     ecp_ortsmo_loss_list[[k]] <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data,
#                                       L = Lecp_ortsmo_k, true_L = L_tilde)
#   }
#   ecp_ortsmo_loss_unlist <- sapply(ecp_ortsmo_loss_list, unlist)            
#   ecp_ortsmo_loss_avg <- rowMeans(ecp_ortsmo_loss_unlist, na.rm = TRUE)
#   
#   ecp_ortsmo_loss <- as.list(ecp_ortsmo_loss_avg)
#   
#   
#   output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, 
#                  "fpca_loss"=fpca_loss, "fpca_loss_unlist"=fpca_loss_unlist,"fpca_ct_loss"=fpca_ct_loss, "fpca_ct_loss_unlist"=fpca_ct_loss_unlist,
#                  "cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp3_ortsmo_loss1"=cp3_ortsmo_loss1, "cp3_ortsmo_loss2"=cp3_ortsmo_loss2, "cp3_ortsmo_loss3"=cp3_ortsmo_loss3, 
#                  "cp6_ortsmo_loss"=cp6_ortsmo_loss, "cp6_ortsmo_loss1"=cp6_ortsmo_loss1,"cp6_ortsmo_loss2"=cp6_ortsmo_loss2, "cp6_ortsmo_loss3"=cp6_ortsmo_loss3,
#                  "cp3_smooth_loss"=cp3_smooth_loss, "cp3_smooth_loss1"=cp3_smooth_loss1, "cp3_smooth_loss2"=cp3_smooth_loss2, "cp3_smooth_loss3"=cp3_smooth_loss3, 
#                  "cp6_smooth_loss"=cp6_smooth_loss, "cp6_smooth_loss1"=cp6_smooth_loss1, "cp6_smooth_loss2"=cp6_smooth_loss2,"cp6_smooth_loss3"=cp6_smooth_loss3,
#                  "mfpca_loss"=mfpca_loss,
#                  "pf3_smooth_loss"=pf3_smooth_loss,
#                  "pf3_ortsmo_loss"=pf3_ortsmo_loss,
#                  "ecp_ortsmo_loss"=ecp_ortsmo_loss)
#   output
# }
# 
# save(fixr_noise01, file="./fixr_new/fixr_noise01.Rda")

## noise level 0.5*empirical var
set.seed(55555)

fixr_noise05 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=200, noise_level=0.5, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss_list <- list() 
  for (k in 1:length(fpca_res$Lfpca_list)){
    Lfpca_k <- fpca_res$Lfpca_list[[k]]
    fpca_loss_list[[k]] <- loss(tnsr = fpca_res$est, true_tnsr=Msmooth@data,
                                L = Lfpca_k, true_L = L_tilde)
  }
  fpca_loss_unlist <- sapply(fpca_loss_list, unlist)            
  fpca_loss_avg <- rowMeans(fpca_loss_unlist, na.rm = TRUE)
  
  fpca_loss <- as.list(fpca_loss_avg)
  
  ## FPCA, fix number of principal components to 3
  fpca_ct_res <- fpca_run(Mmiss, npc=3, center=TRUE)
  
  fpca_ct_loss_list <- list() 
  for (k in 1:length(fpca_ct_res$Lfpca_list)){
    Lfpca_ct_k <- fpca_ct_res$Lfpca_list[[k]]
    fpca_ct_loss_list[[k]] <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data,
                                   L = Lfpca_ct_k, true_L = L_tilde)
  }
  fpca_ct_loss_unlist <- sapply(fpca_ct_loss_list, unlist)            
  fpca_ct_loss_avg <- rowMeans(fpca_ct_loss_unlist, na.rm = TRUE)
  
  fpca_ct_loss <- as.list(fpca_ct_loss_avg)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = cp3_ortsmo_res$Lcp, true_L = L_tilde)
  
  cp3_ortsmo_loss1 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_ortsmo_loss2 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_ortsmo_loss3 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_ortsmo_top3 <- cp6_ortsmo_res$Lcp_ordered[,1:3]
  
  cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = Lcp6_ortsmo_top3, true_L = L_tilde)
  
  cp6_ortsmo_loss1 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_ortsmo_loss2 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_ortsmo_loss3 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_smooth_res <- qr.Q(qr(cp3_smooth_res$Lcp))
  
  cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp3_smooth_res, true_L = L_tilde)
  
  cp3_smooth_loss1 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_smooth_loss2 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_smooth_loss3 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_smooth_top3 <- cp6_smooth_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_smooth_res <- qr.Q(qr(Lcp6_smooth_top3))
  
  cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp6_smooth_res, true_L = L_tilde)
  
  cp6_smooth_loss1 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_smooth_loss2 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_smooth_loss3 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss_list <- list() 
  for (k in 1:length(mfpca_res$Lmfpca_list)){
    Lmfpca_k <- mfpca_res$Lmfpca_list[[k]]
    mfpca_loss_list[[k]] <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data,
                                 L = Lmfpca_k, true_L = L_tilde)
  }
  mfpca_loss_unlist <- sapply(mfpca_loss_list, unlist)            
  mfpca_loss_avg <- rowMeans(mfpca_loss_unlist, na.rm = TRUE)
  
  mfpca_loss <- as.list(mfpca_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
  pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  pf3_smooth_loss_list <- list() 
  for (k in 1:length(pf3_smooth_res$Lparafac2_list)){
    Lpf3_smooth_k <- qr.Q(qr(pf3_smooth_res$Lparafac2_list[[k]]))
    pf3_smooth_loss_list[[k]] <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_smooth_k, true_L = L_tilde)
  }
  pf3_smooth_loss_unlist <- sapply(pf3_smooth_loss_list, unlist)            
  pf3_smooth_loss_avg <- rowMeans(pf3_smooth_loss_unlist, na.rm = TRUE)
  
  pf3_smooth_loss <- as.list(pf3_smooth_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness and orthogonality
  pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  pf3_ortsmo_loss_list <- list() 
  for (k in 1:length(pf3_ortsmo_res$Lparafac2_list)){
    Lpf3_ortsmo_k <- pf3_ortsmo_res$Lparafac2_list[[k]]
    pf3_ortsmo_loss_list[[k]] <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_ortsmo_k, true_L = L_tilde)
  }
  pf3_ortsmo_loss_unlist <- sapply(pf3_ortsmo_loss_list, unlist)            
  pf3_ortsmo_loss_avg <- rowMeans(pf3_ortsmo_loss_unlist, na.rm = TRUE)
  
  pf3_ortsmo_loss <- as.list(pf3_ortsmo_loss_avg)
  
  
  ## SCA-ECP with imputation and smoothness and orthogonality
  ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  ecp_ortsmo_loss_list <- list() 
  for (k in 1:length(ecp_ortsmo_res$Lecp_list)){
    Lecp_ortsmo_k <- ecp_ortsmo_res$Lecp_list[[k]]
    ecp_ortsmo_loss_list[[k]] <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lecp_ortsmo_k, true_L = L_tilde)
  }
  ecp_ortsmo_loss_unlist <- sapply(ecp_ortsmo_loss_list, unlist)            
  ecp_ortsmo_loss_avg <- rowMeans(ecp_ortsmo_loss_unlist, na.rm = TRUE)
  
  ecp_ortsmo_loss <- as.list(ecp_ortsmo_loss_avg)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, 
                 "fpca_loss"=fpca_loss, "fpca_loss_unlist"=fpca_loss_unlist,"fpca_ct_loss"=fpca_ct_loss, "fpca_ct_loss_unlist"=fpca_ct_loss_unlist,
                 "cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp3_ortsmo_loss1"=cp3_ortsmo_loss1, "cp3_ortsmo_loss2"=cp3_ortsmo_loss2, "cp3_ortsmo_loss3"=cp3_ortsmo_loss3, 
                 "cp6_ortsmo_loss"=cp6_ortsmo_loss, "cp6_ortsmo_loss1"=cp6_ortsmo_loss1,"cp6_ortsmo_loss2"=cp6_ortsmo_loss2, "cp6_ortsmo_loss3"=cp6_ortsmo_loss3,
                 "cp3_smooth_loss"=cp3_smooth_loss, "cp3_smooth_loss1"=cp3_smooth_loss1, "cp3_smooth_loss2"=cp3_smooth_loss2, "cp3_smooth_loss3"=cp3_smooth_loss3, 
                 "cp6_smooth_loss"=cp6_smooth_loss, "cp6_smooth_loss1"=cp6_smooth_loss1, "cp6_smooth_loss2"=cp6_smooth_loss2,"cp6_smooth_loss3"=cp6_smooth_loss3,
                 "mfpca_loss"=mfpca_loss,
                 "pf3_smooth_loss"=pf3_smooth_loss,
                 "pf3_ortsmo_loss"=pf3_ortsmo_loss,
                 "ecp_ortsmo_loss"=ecp_ortsmo_loss)
  output
}

save(fixr_noise05, file="./fixr_new/fixr_noise05.Rda")

## noise level 1.5*empirical var
set.seed(1515)

fixr_noise15 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=200, noise_level=1.5, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss_list <- list() 
  for (k in 1:length(fpca_res$Lfpca_list)){
    Lfpca_k <- fpca_res$Lfpca_list[[k]]
    fpca_loss_list[[k]] <- loss(tnsr = fpca_res$est, true_tnsr=Msmooth@data,
                                L = Lfpca_k, true_L = L_tilde)
  }
  fpca_loss_unlist <- sapply(fpca_loss_list, unlist)            
  fpca_loss_avg <- rowMeans(fpca_loss_unlist, na.rm = TRUE)
  
  fpca_loss <- as.list(fpca_loss_avg)
  
  ## FPCA, fix number of principal components to 3
  fpca_ct_res <- fpca_run(Mmiss, npc=3, center=TRUE)
  
  fpca_ct_loss_list <- list() 
  for (k in 1:length(fpca_ct_res$Lfpca_list)){
    Lfpca_ct_k <- fpca_ct_res$Lfpca_list[[k]]
    fpca_ct_loss_list[[k]] <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data,
                                   L = Lfpca_ct_k, true_L = L_tilde)
  }
  fpca_ct_loss_unlist <- sapply(fpca_ct_loss_list, unlist)            
  fpca_ct_loss_avg <- rowMeans(fpca_ct_loss_unlist, na.rm = TRUE)
  
  fpca_ct_loss <- as.list(fpca_ct_loss_avg)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = cp3_ortsmo_res$Lcp, true_L = L_tilde)
  
  cp3_ortsmo_loss1 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_ortsmo_loss2 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_ortsmo_loss3 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_ortsmo_top3 <- cp6_ortsmo_res$Lcp_ordered[,1:3]
  
  cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = Lcp6_ortsmo_top3, true_L = L_tilde)
  
  cp6_ortsmo_loss1 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_ortsmo_loss2 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_ortsmo_loss3 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_smooth_res <- qr.Q(qr(cp3_smooth_res$Lcp))
  
  cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp3_smooth_res, true_L = L_tilde)
  
  cp3_smooth_loss1 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_smooth_loss2 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_smooth_loss3 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_smooth_top3 <- cp6_smooth_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_smooth_res <- qr.Q(qr(Lcp6_smooth_top3))
  
  cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp6_smooth_res, true_L = L_tilde)
  
  cp6_smooth_loss1 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_smooth_loss2 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_smooth_loss3 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss_list <- list() 
  for (k in 1:length(mfpca_res$Lmfpca_list)){
    Lmfpca_k <- mfpca_res$Lmfpca_list[[k]]
    mfpca_loss_list[[k]] <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data,
                                 L = Lmfpca_k, true_L = L_tilde)
  }
  mfpca_loss_unlist <- sapply(mfpca_loss_list, unlist)            
  mfpca_loss_avg <- rowMeans(mfpca_loss_unlist, na.rm = TRUE)
  
  mfpca_loss <- as.list(mfpca_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
  pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  pf3_smooth_loss_list <- list() 
  for (k in 1:length(pf3_smooth_res$Lparafac2_list)){
    Lpf3_smooth_k <- qr.Q(qr(pf3_smooth_res$Lparafac2_list[[k]]))
    pf3_smooth_loss_list[[k]] <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_smooth_k, true_L = L_tilde)
  }
  pf3_smooth_loss_unlist <- sapply(pf3_smooth_loss_list, unlist)            
  pf3_smooth_loss_avg <- rowMeans(pf3_smooth_loss_unlist, na.rm = TRUE)
  
  pf3_smooth_loss <- as.list(pf3_smooth_loss_avg)
  
  ## PARAFAC2-3 with imputation and smoothness and orthogonality
  pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  pf3_ortsmo_loss_list <- list() 
  for (k in 1:length(pf3_ortsmo_res$Lparafac2_list)){
    Lpf3_ortsmo_k <- pf3_ortsmo_res$Lparafac2_list[[k]]
    pf3_ortsmo_loss_list[[k]] <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_ortsmo_k, true_L = L_tilde)
  }
  pf3_ortsmo_loss_unlist <- sapply(pf3_ortsmo_loss_list, unlist)            
  pf3_ortsmo_loss_avg <- rowMeans(pf3_ortsmo_loss_unlist, na.rm = TRUE)
  
  pf3_ortsmo_loss <- as.list(pf3_ortsmo_loss_avg)
  
  
  ## SCA-ECP with imputation and smoothness and orthogonality
  ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  ecp_ortsmo_loss_list <- list() 
  for (k in 1:length(ecp_ortsmo_res$Lecp_list)){
    Lecp_ortsmo_k <- ecp_ortsmo_res$Lecp_list[[k]]
    ecp_ortsmo_loss_list[[k]] <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lecp_ortsmo_k, true_L = L_tilde)
  }
  ecp_ortsmo_loss_unlist <- sapply(ecp_ortsmo_loss_list, unlist)            
  ecp_ortsmo_loss_avg <- rowMeans(ecp_ortsmo_loss_unlist, na.rm = TRUE)
  
  ecp_ortsmo_loss <- as.list(ecp_ortsmo_loss_avg)
  
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, 
                 "fpca_loss"=fpca_loss, "fpca_loss_unlist"=fpca_loss_unlist,"fpca_ct_loss"=fpca_ct_loss, "fpca_ct_loss_unlist"=fpca_ct_loss_unlist,
                 "cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp3_ortsmo_loss1"=cp3_ortsmo_loss1, "cp3_ortsmo_loss2"=cp3_ortsmo_loss2, "cp3_ortsmo_loss3"=cp3_ortsmo_loss3, 
                 "cp6_ortsmo_loss"=cp6_ortsmo_loss, "cp6_ortsmo_loss1"=cp6_ortsmo_loss1,"cp6_ortsmo_loss2"=cp6_ortsmo_loss2, "cp6_ortsmo_loss3"=cp6_ortsmo_loss3,
                 "cp3_smooth_loss"=cp3_smooth_loss, "cp3_smooth_loss1"=cp3_smooth_loss1, "cp3_smooth_loss2"=cp3_smooth_loss2, "cp3_smooth_loss3"=cp3_smooth_loss3, 
                 "cp6_smooth_loss"=cp6_smooth_loss, "cp6_smooth_loss1"=cp6_smooth_loss1, "cp6_smooth_loss2"=cp6_smooth_loss2,"cp6_smooth_loss3"=cp6_smooth_loss3,
                 "mfpca_loss"=mfpca_loss,
                 "pf3_smooth_loss"=pf3_smooth_loss,
                 "pf3_ortsmo_loss"=pf3_ortsmo_loss,
                 "ecp_ortsmo_loss"=ecp_ortsmo_loss)
  output
}

save(fixr_noise15, file="./fixr_new/fixr_noise15.Rda")

## noise level 2*empirical var
set.seed(87654)

fixr_noise2 <- foreach(i = 1:N, .packages = c("rTensor","MASS","refund","multiway","funData","MFPCA","SmoothHOOI")) %dorng%{
  sim_data <- sim_data2(L = L_tilde, R= R_tilde, mean_G = mean_G, cov_G = cov_G, 
                        E = E, p=200, noise_level=2, pattern="random", percent=0.2)
  Mmiss <- sim_data$sim_Mmiss
  Msmooth <- sim_data$sim_Msmooth
  
  ## Oracle
  oracle_opt <- oracle(tnsr=Mmiss@data, smooth_tnsr=Msmooth@data, rank_grid=rank_grid, lambda_seq=lambda_seq,
                       init=0, D = D2, max_iter = 500, tol = 0.1, L0 = NULL)
  
  oracle_res <- mglram(Mmiss@data, ranks=as.numeric(oracle_opt$opt_para[1:2]), lambda=as.numeric(oracle_opt$opt_para[3]),
                       L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  oracle_tilde <- MakeIdent(L=oracle_res$L, G=oracle_res$G, R=oracle_res$R)
  
  oracle_loss <- loss(tnsr = oracle_res$est, true_tnsr=Msmooth@data,
                      L = oracle_tilde$L_tilde, true_L = L_tilde)
  
  ## 10-fold cross-validation
  kcv_opt <- kcv(tnsr=Mmiss@data, rank_grid=rank_grid, lambda_seq=lambda_seq, k=10,
                 L0 = NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_res <- mglram(Mmiss@data, ranks=as.numeric(kcv_opt$opt_para[1:2]), lambda=as.numeric(kcv_opt$opt_para[3]),
                    L0=NULL, D = D2, tol = 0.1, max_iter = 500, init = 0)
  
  kcv_tilde <- MakeIdent(L=kcv_res$L, G=kcv_res$G, R=kcv_res$R)
  
  kcv_loss <- loss(tnsr = kcv_res$est, true_tnsr=Msmooth@data,
                   L = kcv_tilde$L_tilde, true_L = L_tilde)
  
  ## FPCA, fix number of principal components to 3
  fpca_res <- fpca_run(Mmiss, npc=3, center=FALSE)
  
  fpca_loss_list <- list() 
  for (k in 1:length(fpca_res$Lfpca_list)){
    Lfpca_k <- fpca_res$Lfpca_list[[k]]
    fpca_loss_list[[k]] <- loss(tnsr = fpca_res$est, true_tnsr=Msmooth@data,
                                L = Lfpca_k, true_L = L_tilde)
  }
  fpca_loss_unlist <- sapply(fpca_loss_list, unlist)            
  fpca_loss_avg <- rowMeans(fpca_loss_unlist, na.rm = TRUE)
  
  fpca_loss <- as.list(fpca_loss_avg)
  
  ## FPCA, fix number of principal components to 3
  fpca_ct_res <- fpca_run(Mmiss, npc=3, center=TRUE)
  
  fpca_ct_loss_list <- list() 
  for (k in 1:length(fpca_ct_res$Lfpca_list)){
    Lfpca_ct_k <- fpca_ct_res$Lfpca_list[[k]]
    fpca_ct_loss_list[[k]] <- loss(tnsr = fpca_ct_res$est, true_tnsr=Msmooth@data,
                                   L = Lfpca_ct_k, true_L = L_tilde)
  }
  fpca_ct_loss_unlist <- sapply(fpca_ct_loss_list, unlist)            
  fpca_ct_loss_avg <- rowMeans(fpca_ct_loss_unlist, na.rm = TRUE)
  
  fpca_ct_loss <- as.list(fpca_ct_loss_avg)
  
  ## CP3 with imputation and smoothness and orthogonality
  cp3_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  cp3_ortsmo_loss <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = cp3_ortsmo_res$Lcp, true_L = L_tilde)
  
  cp3_ortsmo_loss1 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_ortsmo_loss2 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_ortsmo_loss3 <- loss(tnsr = cp3_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness and orthogonality
  cp6_ortsmo_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("ortsmo", "uncons", "uncons"))
  
  Lcp6_ortsmo_top3 <- cp6_ortsmo_res$Lcp_ordered[,1:3]
  
  cp6_ortsmo_loss <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                          L = Lcp6_ortsmo_top3, true_L = L_tilde)
  
  cp6_ortsmo_loss1 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_ortsmo_loss2 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_ortsmo_loss3 <- loss(tnsr = cp6_ortsmo_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_ortsmo_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP3 with imputation and smoothness, no orthogonality 
  cp3_smooth_res <- cp_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  Lcp_cp3_smooth_res <- qr.Q(qr(cp3_smooth_res$Lcp))
  
  cp3_smooth_loss <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp3_smooth_res, true_L = L_tilde)
  
  cp3_smooth_loss1 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp3_smooth_loss2 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp3_smooth_loss3 <- loss(tnsr = cp3_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp3_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  ## CP6 with imputation and smoothness, no orthogonality
  cp6_smooth_res <- cp_run(Mmiss@data, R_seq=c(6), const=c("smooth", "uncons", "uncons"))
  
  Lcp6_smooth_top3 <- cp6_smooth_res$Lcp_ordered[,1:3]
  
  Lcp_cp6_smooth_res <- qr.Q(qr(Lcp6_smooth_top3))
  
  cp6_smooth_loss <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                          L = Lcp_cp6_smooth_res, true_L = L_tilde)
  
  cp6_smooth_loss1 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,1]), true_L =  as.matrix(L_tilde[,1]))
  
  cp6_smooth_loss2 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,2]), true_L =  as.matrix(L_tilde[,2]))
  
  cp6_smooth_loss3 <- loss(tnsr = cp6_smooth_res$est, true_tnsr=Msmooth@data,
                           L =  as.matrix(cp6_smooth_res$Lcp_ordered[,3]), true_L =  as.matrix(L_tilde[,3]))
  
  
  ## MFPCA 
  mfpca_res <- mfpca_run(Mmiss, M_seq=c(3))
  
  mfpca_loss_list <- list() 
  for (k in 1:length(mfpca_res$Lmfpca_list)){
    Lmfpca_k <- mfpca_res$Lmfpca_list[[k]]
    mfpca_loss_list[[k]] <- loss(tnsr = mfpca_res$est, true_tnsr=Msmooth@data,
                                 L = Lmfpca_k, true_L = L_tilde)
  }
  mfpca_loss_unlist <- sapply(mfpca_loss_list, unlist)            
  mfpca_loss_avg <- rowMeans(mfpca_loss_unlist, na.rm = TRUE)
  
  mfpca_loss <- as.list(mfpca_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness, no orthogonality 
  pf3_smooth_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("smooth", "uncons", "uncons"))
  
  pf3_smooth_loss_list <- list() 
  for (k in 1:length(pf3_smooth_res$Lparafac2_list)){
    Lpf3_smooth_k <- qr.Q(qr(pf3_smooth_res$Lparafac2_list[[k]]))
    pf3_smooth_loss_list[[k]] <- loss(tnsr = pf3_smooth_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_smooth_k, true_L = L_tilde)
  }
  pf3_smooth_loss_unlist <- sapply(pf3_smooth_loss_list, unlist)            
  pf3_smooth_loss_avg <- rowMeans(pf3_smooth_loss_unlist, na.rm = TRUE)
  
  pf3_smooth_loss <- as.list(pf3_smooth_loss_avg)
  
  
  ## PARAFAC2-3 with imputation and smoothness and orthogonality
  pf3_ortsmo_res <- parafac2_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  pf3_ortsmo_loss_list <- list() 
  for (k in 1:length(pf3_ortsmo_res$Lparafac2_list)){
    Lpf3_ortsmo_k <- pf3_ortsmo_res$Lparafac2_list[[k]]
    pf3_ortsmo_loss_list[[k]] <- loss(tnsr = pf3_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lpf3_ortsmo_k, true_L = L_tilde)
  }
  pf3_ortsmo_loss_unlist <- sapply(pf3_ortsmo_loss_list, unlist)            
  pf3_ortsmo_loss_avg <- rowMeans(pf3_ortsmo_loss_unlist, na.rm = TRUE)
  
  pf3_ortsmo_loss <- as.list(pf3_ortsmo_loss_avg)
  
  ## SCA-ECP with imputation and smoothness and orthogonality
  ecp_ortsmo_res <- sca_ecp_run(Mmiss@data, R_seq=c(3), const=c("ortsmo", "uncons", "uncons"))
  
  ecp_ortsmo_loss_list <- list() 
  for (k in 1:length(ecp_ortsmo_res$Lecp_list)){
    Lecp_ortsmo_k <- ecp_ortsmo_res$Lecp_list[[k]]
    ecp_ortsmo_loss_list[[k]] <- loss(tnsr = ecp_ortsmo_res$est, true_tnsr=Msmooth@data,
                                      L = Lecp_ortsmo_k, true_L = L_tilde)
  }
  ecp_ortsmo_loss_unlist <- sapply(ecp_ortsmo_loss_list, unlist)            
  ecp_ortsmo_loss_avg <- rowMeans(ecp_ortsmo_loss_unlist, na.rm = TRUE)
  
  ecp_ortsmo_loss <- as.list(ecp_ortsmo_loss_avg)
  
  output <- list("oracle_loss"=oracle_loss, "kcv_loss"=kcv_loss, 
                 "fpca_loss"=fpca_loss, "fpca_loss_unlist"=fpca_loss_unlist,"fpca_ct_loss"=fpca_ct_loss, "fpca_ct_loss_unlist"=fpca_ct_loss_unlist,
                 "cp3_ortsmo_loss"=cp3_ortsmo_loss, "cp3_ortsmo_loss1"=cp3_ortsmo_loss1, "cp3_ortsmo_loss2"=cp3_ortsmo_loss2, "cp3_ortsmo_loss3"=cp3_ortsmo_loss3, 
                 "cp6_ortsmo_loss"=cp6_ortsmo_loss, "cp6_ortsmo_loss1"=cp6_ortsmo_loss1,"cp6_ortsmo_loss2"=cp6_ortsmo_loss2, "cp6_ortsmo_loss3"=cp6_ortsmo_loss3,
                 "cp3_smooth_loss"=cp3_smooth_loss, "cp3_smooth_loss1"=cp3_smooth_loss1, "cp3_smooth_loss2"=cp3_smooth_loss2, "cp3_smooth_loss3"=cp3_smooth_loss3, 
                 "cp6_smooth_loss"=cp6_smooth_loss, "cp6_smooth_loss1"=cp6_smooth_loss1, "cp6_smooth_loss2"=cp6_smooth_loss2,"cp6_smooth_loss3"=cp6_smooth_loss3,
                 "mfpca_loss"=mfpca_loss,
                 "pf3_smooth_loss"=pf3_smooth_loss,
                 "pf3_ortsmo_loss"=pf3_ortsmo_loss,
                 "ecp_ortsmo_loss"=ecp_ortsmo_loss)
  output
}

save(fixr_noise2, file="./fixr_new/fixr_noise2.Rda")

stopCluster(cl)



