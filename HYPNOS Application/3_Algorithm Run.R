load("./Application/missing_normalized_3.Rda")

library(SmoothHOOI)

# Second Order Difference matrix 
D2 <- SecDiffMat(24)

# Run algorithm with r1=3, r2=2, lambda=4
res <- mglram(tnsr = missing_normalized_3@data, ranks = c(3, 2), init=0, D = D2,
       lambda = 4, max_iter = 500, tol = 1e-5, L0 = NULL)

res$conv

tilde <- MakeIdent(L=res$L, G=res$G, R=res$R)

# Identifiability
L_tilde <- tilde$L_tilde 
R_tilde <- tilde$R_tilde
G_tilde <- tilde$G_tilde

G_mat <- cbind(G_tilde[1,1,], G_tilde[1,2,], 
               G_tilde[2,1,], G_tilde[2,2,],
               G_tilde[3,1,], G_tilde[3,2,])

# mean and covariance matrices of G scores
mean_G <- colMeans(G_mat)
cov_G <- cov(G_mat)

# residuals
E <- missing_normalized_3@data - res$est

#save(L_tilde, R_tilde, G_tilde, E, file="AlgorithmResIdent.Rda")

#save(L_tilde, R_tilde, mean_G, cov_G, E, file="synthetic_base.Rda")
