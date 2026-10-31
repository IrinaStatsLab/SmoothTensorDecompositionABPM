library(multiway)
library(ggplot2)

load("./missing_normalized_3.Rda")
X <- missing_normalized_3@data
X[is.nan(X)] <- NA

nrm <- function(M){                     
  s <- sqrt(colSums(M^2)); s[s < 1e-12] <- 1
  sweep(M, 2, s, "/")
}

match_perm <- function(A, Ref){
  M <- abs(crossprod(nrm(Ref), nrm(A)))
  perm  <- integer(ncol(Ref))
  avail <- seq_len(ncol(A))
  for (j in seq_len(ncol(Ref))) {
    k <- avail[which.max(M[j, avail])]
    perm[j] <- k
    avail <- setdiff(avail, k)
  }
  perm
}

align_to_ref <- function(A, Ref){
  perm   <- match_perm(A, Ref)
  A_perm <- nrm(A)[, perm, drop = FALSE]
  sgn <- sign(diag(crossprod(A_perm, nrm(Ref))))
  sgn[sgn == 0] <- 1
  sweep(A_perm, 2, sgn, "*")
}


# R = 3, ortsmo
R <- 3
nrun <- 10
CONST <- c("ortsmo", "uncons", "uncons")
set.seed(2)
res <- lapply(1:nrun, function(i) parafac(X, nfac = R, const = CONST, output = "best"))

Ref  <- res[[1]]$A
Aaligned  <- lapply(res, function(f) align_to_ref(f$A, Ref))

R <- ncol(Ref)
nrun <- length(Aaligned)

df <- do.call(rbind, lapply(seq_len(nrun), function(i) data.frame(
  hour      = rep(12:35, R),
  loading   = as.vector(Aaligned[[i]]),
  component = rep(paste0("Component ", 1:R), each = 24),
  run       = factor(i)
)))

p <- ggplot(df, aes(hour, loading, group = run, colour = run)) +
  geom_line(alpha = .85) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_grid(~ component) +
  scale_x_continuous(breaks = seq(12, 36, 6), labels = c(12, 18, 0, 6, 12)) +
  labs(x = "Hour of the day", y = "Coefficient") +
  theme_bw() + theme(text = element_text(size = 12))

print(p)

#ggsave("CP3_ortsmo_aligned_10runs.pdf", p, width = 9, height = 3)

# R = 3, smooth
R <- 3
nrun <- 10
CONST <- c("smooth", "uncons", "uncons")
set.seed(1)
res <- lapply(1:nrun, function(i) parafac(X, nfac = R, const = CONST, output = "best"))

Ref  <- res[[1]]$A
Aaligned  <- lapply(res, function(f) align_to_ref(f$A, Ref))

R <- ncol(Ref)
nrun <- length(Aaligned)

df <- do.call(rbind, lapply(seq_len(nrun), function(i) data.frame(
  hour      = rep(12:35, R),
  loading   = as.vector(Aaligned[[i]]),
  component = rep(paste0("Component ", 1:R), each = 24),
  run       = factor(i)
)))

p <- ggplot(df, aes(hour, loading, group = run, colour = run)) +
  geom_line(alpha = .85) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_grid(~ component) +
  scale_x_continuous(breaks = seq(12, 36, 6), labels = c(12, 18, 0, 6, 12)) +
  labs(x = "Hour of the day", y = "Coefficient") +
  theme_bw() + theme(text = element_text(size = 12))

print(p)

##ggsave("CP3_smooth_aligned_10runs.pdf", p, width =9 , height = 3)

# R = 6, ortsmo

R <- 6
nrun <- 10
CONST <- c("ortsmo", "uncons", "uncons")
set.seed(3)
res <- lapply(1:nrun, function(i) parafac(X, nfac = R, const = CONST, output = "best"))

Ref  <- res[[1]]$A
Aaligned  <- lapply(res, function(f) align_to_ref(f$A, Ref))

R <- ncol(Ref)
nrun <- length(Aaligned)

df <- do.call(rbind, lapply(seq_len(nrun), function(i) data.frame(
  hour      = rep(12:35, R),
  loading   = as.vector(Aaligned[[i]]),
  component = rep(paste0("Component ", 1:R), each = 24),
  run       = factor(i)
)))

p <- ggplot(df, aes(hour, loading, group = run, colour = run)) +
  geom_line(alpha = .85) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_grid(~ component) +
  scale_x_continuous(breaks = seq(12, 36, 6), labels = c(12, 18, 0, 6, 12)) +
  labs(x = "Hour of the day", y = "Coefficient") +
  theme_bw() + theme(text = element_text(size = 12))

print(p)

#ggsave("CP6_ortsmo_aligned_10runs.pdf", p, width = 9, height = 3)

# R = 6, smooth
R <- 6
nrun <- 10
CONST <- c("smooth", "uncons", "uncons")
set.seed(4)
res <- lapply(1:nrun, function(i) parafac(X, nfac = R, const = CONST, output = "best"))

Ref  <- res[[1]]$A
Aaligned  <- lapply(res, function(f) align_to_ref(f$A, Ref))

R <- ncol(Ref)
nrun <- length(Aaligned)

df <- do.call(rbind, lapply(seq_len(nrun), function(i) data.frame(
  hour      = rep(12:35, R),
  loading   = as.vector(Aaligned[[i]]),
  component = rep(paste0("Component ", 1:R), each = 24),
  run       = factor(i)
)))

p <- ggplot(df, aes(hour, loading, group = run, colour = run)) +
  geom_line(alpha = .85) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_grid(~ component) +
  scale_x_continuous(breaks = seq(12, 36, 6), labels = c(12, 18, 0, 6, 12)) +
  labs(x = "Hour of the day", y = "Coefficient") +
  theme_bw() + theme(text = element_text(size = 12))

print(p)

#ggsave("CP6_smooth_aligned_10runs.pdf", p, width = 9, height = 3)



