calc_EV_ssesst <- function(M_denom, L_tilde, G_tilde, R_tilde) {
  
  r1 <- ncol(L_tilde)
  r2 <- ncol(R_tilde)
  n  <- dim(G_tilde)[3]
  
  ## SST: uncentered, over observed entries only
  denom <- sum(M_denom[!is.na(M_denom)]^2)
  
  cum_EV_L <- numeric(r1)
  
  for (k in seq_len(r1)) {
    
    sse <- 0
    
    for (i in seq_len(n)) {
      
      M_i <- M_denom[, , i]
      
      L_k <- L_tilde[, seq_len(k), drop = FALSE]
      
      G_k <- matrix(
        G_tilde[seq_len(k), , i],
        nrow = k,
        ncol = r2
      )
      
      fit_i <- L_k %*% G_k %*% t(R_tilde)
      
      obs_i <- which(!is.na(M_i))
      sse <- sse + sum((M_i[obs_i] - fit_i[obs_i])^2)
      
    }
    
    cum_EV_L[k] <- 1 - sse / denom
  }
  
  sep_EV_L <- diff(c(0, cum_EV_L))
  
  
  cum_EV_R <- numeric(r2)
  
  for (k in seq_len(r2)) {
    
    sse <- 0
    
    for (i in seq_len(n)) {
      
      M_i <- M_denom[, , i]
      
      R_k <- R_tilde[, seq_len(k), drop = FALSE]
      
      G_k <- matrix(
        G_tilde[, seq_len(k), i],
        nrow = r1,
        ncol = k
      )
      
      fit_i <- L_tilde %*% G_k %*% t(R_k)
      
      obs_i <- which(!is.na(M_i))
      sse <- sse + sum((M_i[obs_i] - fit_i[obs_i])^2)
      
    }
    
    cum_EV_R[k] <- 1 - sse / denom
  }
  
  sep_EV_R <- diff(c(0, cum_EV_R))
  
  return(
    list(
      cum_EV_L = cum_EV_L,
      sep_EV_L = sep_EV_L,
      cum_EV_R = cum_EV_R,
      sep_EV_R = sep_EV_R
    )
  )
}