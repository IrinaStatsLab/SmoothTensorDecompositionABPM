library(rTensor)

mglram_freeC <- function(tnsr, ranks, lambda, L0 = NULL, C0 = NULL, D,
                         tol = 1e-5, max_iter = 500, init = 0, identify = TRUE) {
  
  ## tnsr: a x b x n (time x measurement x subject)
  ## ranks: c(r1, r2, r3)
  
  if (length(dim(tnsr)) != 3) stop("tnsr must be a three-way array.")
  if (length(ranks) != 3) stop("ranks must be c(r1, r2, r3).")
  if (nrow(D) != ncol(D)) stop("D must be square.")
  
  a <- dim(tnsr)[1]; b <- dim(tnsr)[2]; n <- dim(tnsr)[3]
  r1 <- ranks[1]; r2 <- ranks[2]; r3 <- ranks[3]
  
  if (a != nrow(D))
    stop(sprintf("Time must be mode 1: D is %d x %d but dim(tnsr) = %s.",
                 nrow(D), ncol(D), paste(dim(tnsr), collapse = " x ")))
  if (r1 > a || r2 > b || r3 > n) stop("Specified ranks exceed tensor dimensions.")
  
  ## rank(K) <= min(n, r1*r2)
  r3_max <- min(n, r1 * r2)
  if (r3 > r3_max) {
    warning(sprintf("r3 = %d exceeds min(n, r1*r2) = %d; capping at %d.",
                    r3, r3_max, r3_max))
    r3 <- r3_max
  }
  
  Omega <- !is.na(tnsr)
  M_fill <- tnsr
  M_fill[!Omega] <- init
  has_miss <- !all(Omega)
  
  ## smoothing operators
  D_tilde <- rbind(diag(a), sqrt(lambda) * D)
  A <- diag(a) + lambda * crossprod(D)
  eg <- eigen(A, symmetric = TRUE)
  A_sqrt <- eg$vectors %*% diag(sqrt(eg$values)) %*% t(eg$vectors)
  A_inv_sqrt <- eg$vectors %*% diag(1 / sqrt(eg$values)) %*% t(eg$vectors)
  
  slice <- function(A3, i) A3[, , i, drop = FALSE][, , 1]
  
  ## initialization
  if (is.null(L0)) L0 <- diag(a)[, seq_len(r1), drop = FALSE]
  if (!all(dim(L0) == c(a, r1))) stop("L0 has incorrect dimensions.")
  
  U <- svd(A_sqrt %*% L0, nu = r1, nv = 0)$u
  
  if (is.null(C0)) {
    M_tilde <- array(A_inv_sqrt %*% matrix(M_fill, nrow = a), dim = c(a, b, n))
    unf_subj <- t(matrix(M_tilde, nrow = a * b, ncol = n))
    C0 <- svd(unf_subj, nu = r3, nv = 0)$u
  }
  if (!all(dim(C0) == c(n, r3))) stop("C0 has incorrect dimensions.")
  C <- C0
  
  ## complete-data alternating updates
  inner <- function(Mf, U, C, max_it, tl) {
    o_old <- Inf
    
    for (it in seq_len(max_it)) {
      
      ## N_l = sum_i C_il M_i
      Z <- array(0, dim = c(a, b, r3))
      for (l in seq_len(r3)) {
        for (i in seq_len(n))
          Z[, , l] <- Z[, , l] + C[i, l] * slice(Mf, i)
      }
      
      ## R-step
      SU <- A_inv_sqrt %*% U
      MR <- matrix(0, b, b)
      for (l in seq_len(r3)) {
        H <- crossprod(Z[, , l], SU)
        MR <- MR + tcrossprod(H)
      }
      MR <- (MR + t(MR)) / 2
      R <- eigen(MR, symmetric = TRUE)$vectors[, seq_len(r2), drop = FALSE]
      
      ## U-step
      MU <- matrix(0, a, a)
      for (l in seq_len(r3)) {
        H <- A_inv_sqrt %*% Z[, , l] %*% R
        MU <- MU + tcrossprod(H)
      }
      MU <- (MU + t(MU)) / 2
      U <- eigen(MU, symmetric = TRUE)$vectors[, seq_len(r1), drop = FALSE]
      
      ## C-step
      AU <- A_inv_sqrt %*% U
      K <- matrix(0, n, r1 * r2)
      for (i in seq_len(n))
        K[i, ] <- as.vector(crossprod(AU, slice(Mf, i)) %*% R)
      
      C <- svd(K, nu = r3, nv = 0)$u
      
      ## recover L
      M_L <- A_inv_sqrt %*% tcrossprod(U) %*% A_inv_sqrt
      M_L <- (M_L + t(M_L)) / 2
      L <- eigen(M_L, symmetric = TRUE)$vectors[, seq_len(r1), drop = FALSE]
      
      ## recover core
      Bi <- solve(crossprod(D_tilde %*% L))
      G <- array(0, dim = c(r3, r1, r2))
      
      for (l in seq_len(r3)) {
        Mbar <- matrix(0, a, b)
        for (i in seq_len(n))
          Mbar <- Mbar + C[i, l] * slice(Mf, i)
        G[l, , ] <- Bi %*% crossprod(L, Mbar) %*% R
      }
      
      ## fitted values and penalty
      fit <- pen <- array(0, dim = c(a, b, n))
      DL <- D %*% L
      
      for (i in seq_len(n)) {
        Gam <- matrix(0, r1, r2)
        for (l in seq_len(r3))
          Gam <- Gam + C[i, l] * G[l, , ]
        
        fit[, , i] <- L %*% Gam %*% t(R)
        pen[, , i] <- DL %*% Gam %*% t(R)
      }
      
      o <- sum((Mf - fit)^2) + lambda * sum(pen^2)
      if (is.finite(o_old) && abs(o_old - o) <= tl * max(1, abs(o_old))) break
      o_old <- o
    }
    
    list(L = L, R = R, C = C, G = G, U = U, fit = fit, pen = pen,
         sv_K = svd(K, nu = 0, nv = 0)$d, inner_iter = it)
  }
  
  ## MLSVD identification
  identify_MLSVD <- function(C, L, R, G) {
    unf <- function(G, k) {
      if (k == 1) matrix(G, nrow = r3)
      else if (k == 2) matrix(aperm(G, c(2, 1, 3)), nrow = r1)
      else matrix(aperm(G, c(3, 1, 2)), nrow = r2)
    }
    
    s1 <- svd(unf(G, 1)); s2 <- svd(unf(G, 2)); s3 <- svd(unf(G, 3))
    U1 <- s1$u[, seq_len(r3), drop = FALSE]
    U2 <- s2$u[, seq_len(r1), drop = FALSE]
    U3 <- s3$u[, seq_len(r2), drop = FALSE]
    
    sgn <- function(M) {
      idx <- apply(abs(M), 2, which.max)
      s <- sign(M[cbind(idx, seq_len(ncol(M)))])
      s[s == 0] <- 1
      s
    }
    
    U1 <- sweep(U1, 2, sgn(C %*% U1), "*")
    U2 <- sweep(U2, 2, sgn(L %*% U2), "*")
    U3 <- sweep(U3, 2, sgn(R %*% U3), "*")
    
    Gt <- ttl(as.tensor(G), list(t(U1), t(U2), t(U3)), ms = c(1, 2, 3))@data
    
    relgap <- function(d) {
      if (length(d) < 2) return(NA_real_)
      d <- sort(d, decreasing = TRUE)
      if (max(d) == 0) return(0)
      min(-diff(d)) / max(d)
    }
    
    list(C = C %*% U1, L = L %*% U2, R = R %*% U3, G = Gt,
         sv = list(subject = s1$d, time = s2$d, meas = s3$d),
         relgap = c(subject = relgap(s1$d),
                    time = relgap(s2$d),
                    meas = relgap(s3$d)))
  }
  
  finalize <- function(res, obj_vec, converged) {
    sv_core <- relgap <- NULL
    
    if (identify) {
      idn <- identify_MLSVD(res$C, res$L, res$R, res$G)
      res$C <- idn$C; res$L <- idn$L; res$R <- idn$R; res$G <- idn$G
      sv_core <- idn$sv
      relgap <- idn$relgap
    }
    
    list(L = res$L, R = res$R, C = res$C, G = res$G, U = res$U,
         fitted = res$fit, obj_vec = obj_vec, converged = converged,
         r3 = r3, sv_K = res$sv_K, sv_core = sv_core,
         relgap = relgap, inner_iter = res$inner_iter)
  }
  
  ## complete data
  if (!has_miss) {
    res <- inner(M_fill, U, C, max_iter, tol)
    obj <- sum((M_fill - res$fit)^2) + lambda * sum(res$pen^2)
    return(finalize(res, obj, TRUE))
  }
  
  ## missing-data outer loop
  obj_vec <- rep(NA_real_, max_iter)
  converged <- FALSE
  o_old <- Inf
  
  for (t in seq_len(max_iter)) {
    res <- inner(M_fill, U, C, max_iter, tol)
    U <- res$U; C <- res$C
    
    o <- sum((tnsr[Omega] - res$fit[Omega])^2) + lambda * sum(res$pen^2)
    obj_vec[t] <- o
    
    if (is.finite(o_old) && o > o_old + 1e-8 * max(1, abs(o_old)))
      warning(sprintf("Objective increased at outer iteration %d.", t))
    
    M_fill[!Omega] <- res$fit[!Omega]
    
    if (is.finite(o_old) && abs(o_old - o) <= tol * max(1, abs(o_old))) {
      converged <- TRUE
      break
    }
    
    o_old <- o
  }
  
  finalize(res, obj_vec[seq_len(t)], converged)
}

load("./missing_normalized_3.Rda")

D2 <- SecDiffMat(24)
res <- mglram_freeC(tnsr = missing_normalized_3@data, ranks = c(3, 2, 6), lambda = 4,
                    D = D2, tol = 1e-5, max_iter = 500, init = 0)

res$obj_vec

dataL = data.frame(loading = c(res$L[, 1], res$L[, 2], res$L[, 3]),
                   component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24)),
                   hour = rep(12:35, 3))

pL = dataL %>%
  ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  facet_grid(~ component) + 
  ylab("Coefficient") +
  theme(text = element_text(size = 12)) 
print(pL)

res_SmoothHOOI <- mglram(tnsr = missing_normalized_3@data, ranks = c(3, 2), init=0, D = D2,
                         lambda = 4, max_iter = 500, tol = 1e-5, L0 = NULL)

loss(tnsr=res_SmoothHOOI$est, true_tnsr=res$fitted, 
     L = res_SmoothHOOI$L, true_L = res$L)

res <- mglram_freeC(tnsr = missing_normalized_3@data, ranks = c(3, 2, 3), lambda = 4,
                    D = D2, tol = 1e-5, max_iter = 500, init = 0)

dataL = data.frame(loading = c(res$L[, 1], res$L[, 2], res$L[, 3]),
                   component = c(rep("1st", 24), rep("2nd", 24), rep("3rd", 24)),
                   hour = rep(12:35, 3))

pL = dataL %>%
  ggplot(aes(x = hour, y = loading)) + geom_line(alpha = 0.5) + geom_point() + geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  facet_grid(~ component) + 
  ylab("Coefficient") +
  theme(text = element_text(size = 12)) 
print(pL)

#ggsave("./compress_3rdmode_L.pdf", pL, width = 9, height = 3, dpi = 600)

loss(tnsr=res_SmoothHOOI$est, true_tnsr=res$fitted, 
     L = res_SmoothHOOI$L, true_L = res$L)

nmiss_idx <- which(!is.na(missing_normalized_3@data))
1 - sum((missing_normalized_3@data[nmiss_idx]-res$fitted[nmiss_idx])^2)/sum(missing_normalized_3@data[nmiss_idx]^2)



