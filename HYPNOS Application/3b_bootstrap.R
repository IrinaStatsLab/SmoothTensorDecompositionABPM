library(SmoothHOOI)
library(ggplot2)
library(dplyr)

load("./missing_normalized_3.Rda")
load("./AlgorithmResIdent.Rda") # L_tilde, G_tilde, R_tilde

X <- missing_normalized_3@data

n <- dim(X)[3]
ranks <- c(3, 2)
lambda <- 4
B <- 500              
D2 <- SecDiffMat(24)


# Helper function for aligning booststrap columns to original components 
align_factor <- function(ref, est) {
  r <- ncol(ref)
  
  # All permutations, we don't actuallt need to order because we have addressed identifiability
  # perms <- as.matrix(expand.grid(rep(list(seq_len(r)), r)))
  # perms <- perms[apply(perms, 1, function(p) length(unique(p)) == r),
  #                , drop = FALSE] # remove rows with duplicated numbers
  # 
  # score <- apply(perms, 1, function(p) {
  #   sum(abs(diag(crossprod(ref, est[, p, drop = FALSE]))))
  # }) # inner product of matched columns, abs for magnitude only
  # 
  # p <- perms[which.max(score), ]
  p <- seq_len(r)
  s <- sign(diag(crossprod(ref, est[, p, drop = FALSE]))) # same sign, same direction
  s[s == 0] <- 1
  
  P <- matrix(0, r, r)
  P[cbind(p, seq_len(r))] <- s # align both order and sign 
  
  aligned <- est %*% P
  
  list(aligned = aligned, P = P)
}

# Explained variability on observed entries 
calc_EV <- function(X, L, G, R) {
  
  r1 <- ncol(L)
  r2 <- ncol(R)
  n  <- dim(G)[3]
  
  denom <- sum(X[!is.na(X)]^2)
  
  cum_EV_L <- numeric(r1)
  
  for (k in seq_len(r1)) {
    sse <- 0
    for (i in seq_len(n)) {
      M_i <- X[, , i]
      L_k <- L[, seq_len(k), drop = FALSE]
      G_k <- matrix(G[seq_len(k), , i], nrow = k, ncol = r2)
      fit_i <- L_k %*% G_k %*% t(R)
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
      M_i <- X[, , i]
      R_k <- R[, seq_len(k), drop = FALSE]
      G_k <- matrix(G[, seq_len(k), i], nrow = r1, ncol = k) 
      fit_i <- L %*% G_k %*% t(R_k)
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

ev_original <- calc_EV(X = missing_normalized_3@data, L = L_tilde, G = G_tilde, R = R_tilde)

# Bootstrap 
set.seed(42)

L_boot <- array(NA, c(B, 24, ranks[1]))
R_boot <- array(NA, c(B, 3, ranks[2]))
G_boot <- vector("list", B)

ev_boot <- rep(NA, B)
ev_Lcum_boot <- matrix(NA, B, ranks[1])
ev_Lsep_boot <- matrix(NA, B, ranks[1])
ev_Rcum_boot <- matrix(NA, B, ranks[2])
ev_Rsep_boot <- matrix(NA, B, ranks[2])
boot_indices <- matrix(NA, B, n)
success <- logical(B)

for (b in seq_len(B)) {
  idx <- sample(n, n, replace = TRUE)
  boot_indices[b,] <- idx
  Xb <- X[,,idx, drop = FALSE]
  
  res <- tryCatch({
    fit <- mglram(tnsr = Xb, ranks = ranks, init = 0, D = D2,
                  lambda = lambda, max_iter = 1500,
                  tol = 1e-5, L0 = NULL)
    
    if (!isTRUE(fit$conv)) stop("mglram did not converge")
    
    tilde <- MakeIdent(L = fit$L, G = fit$G, R = fit$R)
    
    # Align L and R to the original fit
    aL <- align_factor(L_tilde, tilde$L_tilde)
    aR <- align_factor(R_tilde, tilde$R_tilde)
    
    # Transform G consistently
    G_aligned <- array(NA, dim(tilde$G_tilde))
    for (i in seq_len(n)) {
      G_aligned[,,i] <- t(aL$P) %*% tilde$G_tilde[,,i] %*% aR$P
    }
    
    ev <- calc_EV(Xb, L = aL$aligned, G = G_aligned, R = aR$aligned)
    
    list(aL = aL$aligned, aR = aR$aligned, aG = G_aligned, ev = ev)
  }, error = function(e) e)
  
  if (inherits(res, "error")) {
    message("Bootstrap ", b, ": ", conditionMessage(res))
    next
  }
  
  L_boot[b,,] <- res$aL
  R_boot[b,,] <- res$aR
  G_boot[[b]] <- res$aG
  ev_boot[b] <- tail(res$ev$cum_EV_L,1)
  ev_Lcum_boot[b,] <- res$ev$cum_EV_L
  ev_Lsep_boot[b,] <- res$ev$sep_EV_L
  ev_Rcum_boot[b,] <- res$ev$cum_EV_R
  ev_Rsep_boot[b,] <- res$ev$sep_EV_R
  success[b] <- TRUE
  
  if (b %% 20 == 0) message("Completed ", b, "/", B)
}

cat("Successful bootstrap fits:", sum(success), "/", B, "\n") 


# save(
#   L_boot, R_boot, G_boot,
#   ev_boot,
#   ev_Lcum_boot, ev_Lsep_boot,
#   ev_Rcum_boot, ev_Rsep_boot,
#   success, boot_indices,
#   file="bootstrap_SmoothHOOI.Rda")

# load("./bootstrap_SmoothHOOI.Rda")

sum(success) # 500

# 95% bootstrap confidence intervals for L and EV
L_ci <- apply(L_boot, c(2, 3), quantile,
              probs = c(0.025, 0.975), na.rm = TRUE)

# Plot coefficients for L components with 95% bootstrap bands
dataL <- data.frame(
  loading = c(L_tilde[,1], L_tilde[,2], L_tilde[,3]),
  lower =  as.vector(L_ci[1,,]),
  upper = as.vector(L_ci[2,,]),
  component = factor(rep(c("1st", "2nd", "3rd"), each = 24),
                     levels = c("1st", "2nd", "3rd")),
  hour = rep(12:35, 3)
)

pL <- dataL %>%
  ggplot(aes(x = hour, y = loading)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#6BAED6", alpha = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, color = "grey45", linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36),
    labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)
  ) +
  xlab("Hour of the day") +
  facet_grid(~ component) +
  ylab("Coefficient") +
  theme(text = element_text(size = 15)) +
  coord_cartesian(ylim = c(-0.4, 0.4))

pL

# ggsave("./L_comp_95CI.pdf", dpi=600, width=8, height=2.5)

ev_ci <- quantile(ev_boot[success], probs = c(0.025, 0.975), na.rm = TRUE)

ev_Lsep_ci <- apply(ev_Lsep_boot[success, ], 2, quantile, probs = c(0.025, 0.975))
ev_Rsep_ci <- apply(ev_Rsep_boot[success, ], 2, quantile, probs = c(0.025, 0.975))

round(ev_ci, 4)
round(ev_Lsep_ci, 4)
round(ev_Rsep_ci, 4)


