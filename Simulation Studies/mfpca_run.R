mfpca_run <- function(tnsr, M_seq){
  nmiss_idx <- which(!is.na(tnsr@data))
  
  a <- tnsr@modes[1]
  b <- tnsr@modes[2]
  
  argvals_grid <- seq(1, a, by = 1)
  
  fd_list <- list()
  # iterate over the mode 2 of the tensor and make them into matrices
  for (i in 1:b){
    assign(paste0("matrix_", i), t(tnsr@data[,i,]))
    assign(paste0("fd_", i), funData(argvals = list(argvals_grid), X = get(paste0("matrix_", i))))
    fd_list[[i]] <- get(paste0("fd_", i))
  }
  
  mfd_all <- multiFunData(fd_list)
  
  uniExpansions <- replicate(b, list(type = "uFPCA"), simplify = FALSE)
  
  if (length(M_seq)==1){
    MFPCAres <- MFPCA(mfd_all, M = M_seq[1], uniExpansions = uniExpansions, fit=TRUE)
    eigenfunctions <- MFPCAres$functions
    reconstruction <- MFPCAres$fit
    
    est <- array(NA, dim=tnsr@modes)
    for (i in 1:b){
      est[,i,] <- t(reconstruction[[i]]@X)
    }
    
    L_list <- list()
    for (i in 1:b){
      Lmfpca_i <- t(eigenfunctions[[i]]@X)
      L_list[[i]] <- qr.Q(qr(Lmfpca_i))
    }
    
    Lmfpca <- Reduce("+", L_list) / length(L_list)
    Lmfpca <- qr.Q(qr(Lmfpca))
    
    out=list(est=est, Lmfpca=Lmfpca)
    
  } else if (length(M_seq)>1){
    exp_var = c()
    est = list()
    Lmfpca_list = list()
    for (i in 1:length(M_seq)){
      M_i <- M_seq[i]
      MFPCAres <- MFPCA(mfd_all, M = M_i, uniExpansions = uniExpansions, fit=TRUE)
      eigenfunctions <- MFPCAres$functions
      reconstruction <- MFPCAres$fit
      
      est_i <- array(NA, dim=tnsr@modes)
      for (j in 1:b){
        est_i[,j,] <- t(reconstruction[[j]]@X)
      }
      
      exp_var[i] <- sum((est_i[nmiss_idx])^2)/sum((tnsr@data[nmiss_idx])^2)
      est[[i]] <- est_i
      
      L_list <- list()
      for (j in 1:b){
        Lmfpca_i <- t(eigenfunctions[[j]]@X)
        L_list[[j]] <- qr.Q(qr(Lmfpca_i))
      }
      
      Lmfpca_i <- Reduce("+", L_list) / length(L_list)
      Lmfpca_i <- qr.Q(qr(Lmfpca_i))
      Lmfpca_list[[i]] <- Lmfpca_i
    }
    best_M_idx <- which.max(exp_var)
    best_M <- M_seq[best_M_idx]
    best_est <- est[[best_M_idx]]
    best_Lmfpca <- Lmfpca_list[[best_M_idx]]
    
    out=list(est=best_est, Lmfpca=best_Lmfpca)
  }
}


