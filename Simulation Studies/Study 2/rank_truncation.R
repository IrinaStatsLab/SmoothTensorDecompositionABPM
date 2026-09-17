rank_truncation <- function(tnsr_array, mode, true_rank_L=NULL, true_rank_R=NULL, L_tilde, G_tilde, R_tilde){
  n <- dim(tnsr_array@data)[3]
  if (mode=="L"){
    L_truncated <- L_tilde[,1:true_rank_L]
    
    est_truncated <- array(NA, dim(tnsr_array@data))
    for (i in 1:n){
      est_truncated[ , , i] <- matrix(L_tilde[, 1:true_rank_L], ncol=true_rank_L) %*% matrix(G_tilde[1:true_rank_L, , i], nrow=true_rank_L) %*% t(R_tilde)
    }
    
    out=list(est_truncated=est_truncated, L_truncated=L_truncated)
    return(out)
    
  } else if (mode=="R"){
    R_truncated <- R_tilde[,1:true_rank_R]
    est_truncated <- array(NA, dim(tnsr_array@data))
    for (i in 1:n){
      est_truncated[ , , i] <- L_tilde %*% matrix(G_tilde[,1:true_rank_R, i], ncol=true_rank_R) %*% matrix(t(R_tilde[, 1:true_rank_R]), nrow=true_rank_R)
    }
    
    out=list(est_truncated=est_truncated, R_truncated=R_truncated)
    return(out)
    
  } else if (mode=="both"){
    L_truncated = L_tilde[,1:true_rank_L]
    R_truncated = R_tilde[,1:true_rank_R]
    
    est_truncated <- array(NA, dim(tnsr_array@data))
    for (i in 1:n){
      est_truncated[ , , i] <- L_truncated %*% G_tilde[1:true_rank_L,1:true_rank_R, i] %*% t(R_truncated)
    }
    
    out=list(est_truncated=est_truncated, L_truncated=L_truncated, R_truncated=R_truncated)
    return(out)
    
  }
  
}

