load("./filtered_ABPM_new.Rda") # filtered_abpm_new
load("./combinedCGMmetricsCovariates.Rda") # combined
load("./AlgorithmResIdent.Rda")

library(rTensor)
library(dplyr)
library(ggplot2)
library(car)
library(tidyverse)
library(broom)
library(readxl)
library(gtsummary)
library(knitr)
library(kableExtra)

# Data preparation for linear regression
id_207 <- unique(filtered_abpm_new$ID) # 207 ids
id_cov <- unique(combined$id) # 204 ids
match <- match(id_cov, id_207) 
matched_id <- id_207[match[!is.na(match)]] # which ids are in both data

## indexes of matached ids in the tensor
matched_idx <- which(unique(filtered_abpm_new$ID) %in% matched_id)

## covariates data for matched ids
covariates <- combined[combined$id %in% matched_id, ] 
xpred <- covariates[,-c(2:4)]

## Make race=white and gender=male as reference levels
xpred$race_White_inv <- 1 - xpred$race_White 
xpred$gender_inv <- 1 - xpred$gender

# Regression Model 1
scoredata1  <- cbind(G_tilde[1, 1, matched_idx], xpred)
names(scoredata1)[1] = "TensorScore1"
out1 = lm(TensorScore1 ~ age + BMI + gender_inv + race_White_inv + odi4, data = scoredata1)
summary(out1) 

# Regression Model 2
scoredata2 = cbind(G_tilde[1, 2, matched_idx], xpred)
names(scoredata2)[1] = "TensorScore2"
out2 = lm(TensorScore2 ~ age + BMI + gender_inv + race_White_inv + odi4, data = scoredata2)
summary(out2) 


Ri_191 =  ttm(as.tensor(G_tilde[, , matched_idx]), R_tilde, 2)

RiL1_191 = Ri_191@data[1, , ]

RiL2_191 = Ri_191@data[2, , ]

RiL3_191 = Ri_191@data[3, , ]

DBPmean = rep(70, 24) + (L_tilde[, 1] * mean(RiL1_191[1,]) + L_tilde[, 2] * mean(RiL2_191[1,]) + L_tilde[, 3] * mean(RiL3_191[1,])) * 12 # 12 is the sd of DBP, doing this because we normalized the data initially
SBPmean = rep(129, 24) + (L_tilde[, 1] * mean(RiL1_191[2,]) + L_tilde[, 2] * mean(RiL2_191[2,]) + L_tilde[, 3] * mean(RiL3_191[2,])) * 18 # 18 is the sd of SBP
HRmean = rep(76, 24) + (L_tilde[, 1] * mean(RiL1_191[3,]) + L_tilde[, 2] * mean(RiL2_191[3,]) + L_tilde[, 3] * mean(RiL3_191[3,])) * 13 # 13 is the sd of HR 

# Scaling
sdD = c(12, 18, 13) 

# Visualize regression model results (Age)

# Main effects for each model 
C1 = (L_tilde[, 1] %*% t(R_tilde[, 1])) %*% diag(sdD)
C2 = (L_tilde[, 1] %*% t(R_tilde[, 2])) %*% diag(sdD)

effectEst_all <- lapply(1:5, function(x) coef(out1)[x+1]*C1 + coef(out2)[x+1]*C2)

fullSigma1_all <- lapply(1:5, function(x) (summary(out1)$coefficients[, "Std. Error"][x+1]^2) * (C1[ , 1] %*% t(C1[, 1])) + (summary(out2)$coefficients[, "Std. Error"][x+1]^2) * (C2[ , 1] %*% t(C2[, 1])))

fullSigma2_all <- lapply(1:5, function(x) (summary(out1)$coefficients[, "Std. Error"][x+1]^2) * (C1[ , 2] %*% t(C1[, 2])) + (summary(out2)$coefficients[, "Std. Error"][x+1]^2) * (C2[ , 2] %*% t(C2[, 2])))    

fullSigma3_all <- lapply(1:5, function(x) (summary(out1)$coefficients[, "Std. Error"][x+1]^2) * (C1[ , 3] %*% t(C1[, 3])) + (summary(out2)$coefficients[, "Std. Error"][x+1]^2) * (C2[ , 3] %*% t(C2[, 3])))    

time_points <- (1:nrow(effectEst_all[[1]])) + 12

df_long <- data.frame(
  time = rep(time_points, 3),
  outcome = rep(c("DBP", "SBP", "HR"), each = length(time_points)),
  effect = c(effectEst_all[[1]][,1], effectEst_all[[1]][,2], effectEst_all[[1]][,3]),
  lower = c(effectEst_all[[1]][,1] - 1.96 * sqrt(diag(fullSigma1_all[[1]])),
            effectEst_all[[1]][,2] - 1.96 * sqrt(diag(fullSigma2_all[[1]])),
            effectEst_all[[1]][,3] - 1.96 * sqrt(diag(fullSigma3_all[[1]]))),
  upper = c(effectEst_all[[1]][,1] + 1.96 * sqrt(diag(fullSigma1_all[[1]])),
            effectEst_all[[1]][,2] + 1.96 * sqrt(diag(fullSigma2_all[[1]])),
            effectEst_all[[1]][,3] + 1.96 * sqrt(diag(fullSigma3_all[[1]])))
)

df_long$outcome = factor(df_long$outcome, levels = c("DBP","SBP","HR"))

## effect size plot
ggplot(df_long, aes(x = time, y = effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(color = "blue") +
  facet_wrap(~ outcome, nrow = 1) +
  ylim(-0.6, 0.3) +
  labs(x = "Hour of the day", y = "Effect Size of Age") +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
    labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)
  ) +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15))

# ggsave("effect_of_age.pdf", dpi=600, width = 8, height=3)

# Visualize regression model results (BMI)

df_long <- data.frame(
  hour = rep(seq(12,35, by=1), 3),
  outcome = rep(c("DBP", "SBP", "HR"), each = length(seq(12,35, by=1))),
  effect = c(effectEst_all[[2]][,1], effectEst_all[[2]][,2], effectEst_all[[2]][,3]),
  lower = c(effectEst_all[[2]][,1] - 1.96 * sqrt(diag(fullSigma1_all[[2]])),
            effectEst_all[[2]][,2] - 1.96 * sqrt(diag(fullSigma2_all[[2]])),
            effectEst_all[[2]][,3] - 1.96 * sqrt(diag(fullSigma3_all[[2]]))),
  upper = c(effectEst_all[[2]][,1] + 1.96 * sqrt(diag(fullSigma1_all[[2]])),
            effectEst_all[[2]][,2] + 1.96 * sqrt(diag(fullSigma2_all[[2]])),
            effectEst_all[[2]][,3] + 1.96 * sqrt(diag(fullSigma3_all[[2]])))
)

df_long$outcome = factor(df_long$outcome, levels = c("DBP", "SBP", "HR"))

## effect size plot
ggplot(df_long, aes(x = hour, y = effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(color = "blue") +
  facet_wrap(~ outcome, nrow = 1) +
  ylim(-0.5, 1) +
  labs(x = "Hour of the day", y = "Effect Size of BMI") +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
    labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)
  ) +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15))

#ggsave("effect_of_bmi.pdf", dpi=600, width = 8, height=3)


# Visualize regression model results (Sex)

df_long <- data.frame(
  hour = rep(seq(12,35, by=1), 3),
  outcome = rep(c("DBP", "SBP", "HR"), each = length(seq(12,35, by=1))),
  effect = c(effectEst_all[[3]][,1], effectEst_all[[3]][,2], effectEst_all[[3]][,3]),
  lower = c(effectEst_all[[3]][,1] - 1.96 * sqrt(diag(fullSigma1_all[[3]])),
            effectEst_all[[3]][,2] - 1.96 * sqrt(diag(fullSigma2_all[[3]])),
            effectEst_all[[3]][,3] - 1.96 * sqrt(diag(fullSigma3_all[[3]]))),
  upper = c(effectEst_all[[3]][,1] + 1.96 * sqrt(diag(fullSigma1_all[[3]])),
            effectEst_all[[3]][,2] + 1.96 * sqrt(diag(fullSigma2_all[[3]])),
            effectEst_all[[3]][,3] + 1.96 * sqrt(diag(fullSigma3_all[[3]])))
)

df_long$outcome = factor(df_long$outcome, levels = c("DBP", "SBP", "HR"))

## effect size plot
ggplot(df_long, aes(x = hour, y = effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(color = "blue") +
  facet_wrap(~ outcome, nrow = 1) +
  ylim(-7.5, 7.5) +
  labs(x = "Hour of the day", y = "Effect Size of Sex") +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
    labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)
  ) +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15))

# ggsave("effect_of_sex.pdf", dpi=600, width = 8, height=3)


# Visualize regression model results (Race)

df_long <- data.frame(
  hour = rep(seq(12,35, by=1), 3),
  outcome = rep(c("DBP", "SBP", "HR"), each = length(seq(12,35, by=1))),
  effect = c(effectEst_all[[4]][,1], effectEst_all[[4]][,2], effectEst_all[[4]][,3]),
  lower = c(effectEst_all[[4]][,1] - 1.96 * sqrt(diag(fullSigma1_all[[4]])),
            effectEst_all[[4]][,2] - 1.96 * sqrt(diag(fullSigma2_all[[4]])),
            effectEst_all[[4]][,3] - 1.96 * sqrt(diag(fullSigma3_all[[4]]))),
  upper = c(effectEst_all[[4]][,1] + 1.96 * sqrt(diag(fullSigma1_all[[4]])),
            effectEst_all[[4]][,2] + 1.96 * sqrt(diag(fullSigma2_all[[4]])),
            effectEst_all[[4]][,3] + 1.96 * sqrt(diag(fullSigma3_all[[4]])))
)

df_long$outcome = factor(df_long$outcome, levels = c("DBP", "SBP", "HR"))

## effect size plot
ggplot(df_long, aes(x = hour, y = effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(color = "blue") +
  facet_wrap(~ outcome, nrow = 1) +
  ylim(-2.5, 10) +
  labs(x = "Hour of the day", y = "Effect Size of Race") +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
    labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)
  ) +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15))

# ggsave("effect_of_race.pdf", dpi=600, width = 8, height=3)


# Visualize regression model results (ODI4)

df_long <- data.frame(
  hour = rep(seq(12,35, by=1), 3),
  outcome = rep(c("DBP", "SBP", "HR"), each = length(seq(12,35, by=1))),
  effect = c(effectEst_all[[5]][,1], effectEst_all[[5]][,2], effectEst_all[[5]][,3]),
  lower = c(effectEst_all[[5]][,1] - 1.96 * sqrt(diag(fullSigma1_all[[5]])),
            effectEst_all[[5]][,2] - 1.96 * sqrt(diag(fullSigma2_all[[5]])),
            effectEst_all[[5]][,3] - 1.96 * sqrt(diag(fullSigma3_all[[5]]))),
  upper = c(effectEst_all[[5]][,1] + 1.96 * sqrt(diag(fullSigma1_all[[5]])),
            effectEst_all[[5]][,2] + 1.96 * sqrt(diag(fullSigma2_all[[5]])),
            effectEst_all[[5]][,3] + 1.96 * sqrt(diag(fullSigma3_all[[5]])))
)

df_long$outcome = factor(df_long$outcome, levels = c("DBP", "SBP", "HR"))

# effect size plot
ggplot(df_long, aes(x = hour, y = effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_line(color = "blue") +
  facet_wrap(~ outcome, nrow = 1) +
  ylim(-0.5, 7.5) +
  labs(x = "Hour of the day", y = "Effect Size of log(ODI4)") +
  scale_x_continuous(
    breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
    labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)
  ) +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15))

# ggsave("effect_of_odi4.pdf", dpi=600, width = 8, height=3)



