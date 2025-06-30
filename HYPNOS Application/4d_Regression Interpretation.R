load("./Application/filtered_ABPM_new.Rda") # filtered_abpm_new
load("./Application/combinedCGMmetricsCovariates.Rda") # combined
load("./Application/AlgorithmResIdent.Rda")

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

# ggsave("./Application/effect_of_age.pdf", dpi=600, width = 8, height=3)

DBPmean_age_low =  DBPmean - effectEst_all[[1]][,1]*10
DBPmean_age_up =  DBPmean + effectEst_all[[1]][,1]*10

SBPmean_age_low =  SBPmean - effectEst_all[[1]][,2]*10
SBPmean_age_up =  SBPmean + effectEst_all[[1]][,2]*10

HRmean_age_low = HRmean - effectEst_all[[1]][,3]*10
HRmean_age_up = HRmean + effectEst_all[[1]][,3]*10

age_data = data.frame(
  mean = c(rep(DBPmean,2), rep(SBPmean,2), rep(HRmean, 2)),
  value = c(DBPmean_age_low, DBPmean_age_up, 
            SBPmean_age_low, SBPmean_age_up, 
            HRmean_age_low, HRmean_age_up),
  age_level = rep(rep(c("-10 years", "+10 years"), each = 24), 6),
  hour = rep(12:35, 12),
  group = c(rep("DBP",24*2),rep("SBP",24*2),rep("HR",24*2)))
age_data$group <- factor(age_data$group, levels = c("DBP", "SBP", "HR"))

pdataL <- age_data %>%
  ggplot(aes(x = hour, y = value, group = age_level, col= age_level)) +
  geom_line(alpha = 1) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(color = "Age Change") +
  geom_line(aes(x = hour, y = mean), color="black") +
  facet_wrap(vars(group), scales="free_y") +
  theme(text = element_text(size = 15)) 
pdataL 

# ggsave("./Application/age_change_plot.pdf", dpi=600, width = 8, height=3)

## estimation: median age + others at mean levels
g1_median = -5.258 - 0.057*median(xpred$age) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)
g2_median = 10.129 -0.110*median(xpred$age) - 0.124*mean(xpred$BMI) + 1.384*mean(xpred$gender_inv) - 0.505*mean(xpred$race_White_inv) + 0.060*mean(xpred$odi4)

g1_q1 = -5.258 - 0.057*quantile(xpred$age, 0.25) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)
g2_q1 = 10.129 -0.110*quantile(xpred$age, 0.25) - 0.124*mean(xpred$BMI) + 1.384*mean(xpred$gender_inv) - 0.505*mean(xpred$race_White_inv) + 0.060*mean(xpred$odi4)

g1_q3 = -5.258 - 0.057*quantile(xpred$age, 0.75) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)
g2_q3 = -5.258 - 0.057*quantile(xpred$age, 0.75) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)

age_median <- g1_median * C1 + g2_median * C2
age_q1 <- g1_q1 * C1 + g2_q1 * C2
age_q3 <- g1_q3 * C1 + g2_q3 * C2

DBPmean_noL1 = rep(70, 24) + (L_tilde[, 2] * mean(RiL2_191[1,]) + L_tilde[, 3] * mean(RiL3_191[1,])) * 12
SBPmean_noL1 = rep(129, 24) + (L_tilde[, 2] * mean(RiL2_191[2,]) + L_tilde[, 3] * mean(RiL3_191[2,])) * 18 # 18 is the sd of SBP
HRmean_noL1 = rep(76, 24) + (L_tilde[, 2] * mean(RiL2_191[3,]) + L_tilde[, 3] * mean(RiL3_191[3,])) * 13 # 13 is the sd of HR 

age_data = data.frame(
  value = c(DBPmean_noL1+age_median[,1], DBPmean_noL1+age_q1[,1], DBPmean_noL1+age_q3[,1], 
            SBPmean_noL1+age_median[,2], SBPmean_noL1+age_q1[,2], SBPmean_noL1+age_q3[,2], 
            HRmean_noL1+age_median[,3], HRmean_noL1+age_q1[,3], HRmean_noL1+age_q3[,3]),
  age_level = rep(rep(c("median","1st quartile", "3rd quartile"), each = 24), 3),
  hour = rep(12:35, 9),
  group = c(rep("DBP",24*3),rep("SBP",24*3),rep("HR",24*3)))
age_data$group <- factor(age_data$group, levels = c("DBP", "SBP", "HR"))
age_data$odi4_level <- factor(age_data$age_level, levels = c("1st quartile", "median", "3rd quartile"))

p <- age_data %>%
  ggplot(aes(x = hour, y = value, group = age_level, col = odi4_level)) +
  geom_line(alpha = 1) + 
  scale_color_manual(values = c("1st quartile" = "skyblue", "median" = "#4292C6", "3rd quartile" = "#2066A8"),
                     labels = c("1st quartile" = "Q1: 53", "median" = "median: 61", "3rd quartile" = "Q3: 67")) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(color = "Quantile") +
  facet_wrap(vars(group), scales = "free_y") +
  theme(text = element_text(size = 15)) 
p

# ggsave("./Application/age_change_plot.pdf", dpi=600, width = 9, height=3)

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

#ggsave("./Application/effect_of_bmi.pdf", dpi=600, width = 8, height=3)

DBPmean_bmi_low =  DBPmean - effectEst_all[[2]][,1]*10
DBPmean_bmi_up =  DBPmean + effectEst_all[[2]][,1]*10

SBPmean_bmi_low =  SBPmean - effectEst_all[[2]][,2]*10
SBPmean_bmi_up =  SBPmean + effectEst_all[[2]][,2]*10

HRmean_bmi_low = HRmean - effectEst_all[[2]][,3]*10
HRmean_bmi_up = HRmean + effectEst_all[[2]][,3]*10

bmi_data = data.frame(
  mean = c(rep(DBPmean,2), rep(SBPmean,2), rep(HRmean, 2)),
  value = c(DBPmean_bmi_low, DBPmean_bmi_up, 
            SBPmean_bmi_low, SBPmean_bmi_up, 
            HRmean_bmi_low, HRmean_bmi_up),
  bmi_level = rep(rep(c("-10", "+10"), each = 24), 6),
  hour = rep(12:35, 12),
  group = c(rep("DBP",24*2),rep("SBP",24*2),rep("HR",24*2)))
bmi_data$group <- factor(bmi_data$group, levels = c("DBP", "SBP", "HR"))

pdataL <- bmi_data %>%
  ggplot(aes(x = hour, y = value, group = bmi_level, col= bmi_level)) +
  geom_line(alpha = 1) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(color = "BMI Change") +
  geom_line(aes(x = hour, y = mean), color="black") +
  facet_wrap(vars(group), scales="free_y") +
  theme(text = element_text(size = 15)) 
pdataL

# ggsave("./Application/bmi_change_plot.pdf", dpi=600, width = 8, height=3)

## estimation: median BMI + others at mean levels
g1_median = -5.258 - 0.057*mean(xpred$age) + 0.106*median(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)
g2_median = 10.129 -0.110*mean(xpred$age) - 0.124*median(xpred$BMI) + 1.384*mean(xpred$gender_inv) - 0.505*mean(xpred$race_White_inv) + 0.060*mean(xpred$odi4)

g1_q1 = -5.258 - 0.057*mean(xpred$age) + 0.106*quantile(xpred$BMI, 0.25) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)
g2_q1 = 10.129 -0.110*mean(xpred$age) - 0.124*quantile(xpred$BMI, 0.25) + 1.384*mean(xpred$gender_inv) - 0.505*mean(xpred$race_White_inv) + 0.060*mean(xpred$odi4)

g1_q3 = -5.258 - 0.057*mean(xpred$age) + 0.106*quantile(xpred$BMI, 0.75) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)
g2_q3 = -5.258 - 0.057*mean(xpred$age) + 0.106*quantile(xpred$BMI, 0.75) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*mean(xpred$odi4)

bmi_median <- g1_median * C1 + g2_median * C2
bmi_q1 <- g1_q1 * C1 + g2_q1 * C2
bmi_q3 <- g1_q3 * C1 + g2_q3 * C2

DBPmean_noL1 = rep(70, 24) + (L_tilde[, 2] * mean(RiL2_191[1,]) + L_tilde[, 3] * mean(RiL3_191[1,])) * 12
SBPmean_noL1 = rep(129, 24) + (L_tilde[, 2] * mean(RiL2_191[2,]) + L_tilde[, 3] * mean(RiL3_191[2,])) * 18 # 18 is the sd of SBP
HRmean_noL1 = rep(76, 24) + (L_tilde[, 2] * mean(RiL2_191[3,]) + L_tilde[, 3] * mean(RiL3_191[3,])) * 13 # 13 is the sd of HR 

bmi_data = data.frame(
  value = c(DBPmean_noL1+bmi_median[,1], DBPmean_noL1+bmi_q1[,1], DBPmean_noL1+bmi_q3[,1], 
            SBPmean_noL1+bmi_median[,2], SBPmean_noL1+bmi_q1[,2], SBPmean_noL1+bmi_q3[,2], 
            HRmean_noL1+bmi_median[,3], HRmean_noL1+bmi_q1[,3], HRmean_noL1+bmi_q3[,3]),
  bmi_level = rep(rep(c("median","1st quartile", "3rd quartile"), each = 24), 3),
  hour = rep(12:35, 9),
  group = c(rep("DBP",24*3),rep("SBP",24*3),rep("HR",24*3)))
bmi_data$group <- factor(bmi_data$group, levels = c("DBP", "SBP", "HR"))
bmi_data$odi4_level <- factor(bmi_data$bmi_level, levels = c("1st quartile", "median", "3rd quartile"))

p <- bmi_data %>%
  ggplot(aes(x = hour, y = value, group = bmi_level, col = odi4_level)) +
  geom_line(alpha = 1) + 
  scale_color_manual(values = c("1st quartile" = "skyblue", "median" = "#4292C6", "3rd quartile" = "#2066A8"),
                     labels = c("1st quartile" = "Q1: 30.16", "median" = "median: 32.89", "3rd quartile" = "Q3: 36.95")) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(color = "Quantile") +
  facet_wrap(vars(group), scales = "free_y") +
  theme(text = element_text(size = 15)) 
p

# ggsave("./Application/bmi_change_plot.pdf", dpi=600, width = 9, height=3)

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

# ggsave("./Application/effect_of_sex.pdf", dpi=600, width = 8, height=3)

## Estimation for male vs female
g1_male = -5.258 - 0.057*mean(xpred$age[xpred$gender_inv==0]) + 0.106*mean(xpred$BMI[xpred$gender_inv==0]) + 1.948*mean(xpred$race_White_inv[xpred$gender_inv==0]) + 1.615*mean(xpred$odi4[xpred$gender_inv==0])
g2_male = 10.129 -0.110*mean(xpred$age[xpred$gender_inv==0]) - 0.124*mean(xpred$BMI[xpred$gender_inv==0]) - 0.505*mean(xpred$race_White_inv[xpred$gender_inv==0]) + 0.060*mean(xpred$odi4[xpred$gender_inv==0])
g1_female = -5.258 - 0.057*mean(xpred$age[xpred$gender_inv==1]) + 0.106*mean(xpred$BMI[xpred$gender_inv==1]) - 0.275 + 1.948*mean(xpred$race_White_inv[xpred$gender_inv==1]) + 1.615*mean(xpred$odi4[xpred$gender_inv==1]) 
g2_female = 10.129 -0.110*mean(xpred$age[xpred$gender_inv==1]) - 0.124*mean(xpred$BMI[xpred$gender_inv==1]) + 1.384 - 0.505*mean(xpred$race_White_inv[xpred$gender_inv==0]) + 0.060*mean(xpred$odi4[xpred$gender_inv==1]) 

male <- g1_male * C1 + g2_male * C2
female <- g1_female * C1 + g2_female * C2

sex_data = data.frame(
  mean = c(rep(DBPmean,2), rep(SBPmean,2), rep(HRmean, 2)),
  value = c(DBPmean_noL1 + male[,1], DBPmean_noL1 + female[,1], 
            SBPmean_noL1 + male[,2], SBPmean_noL1 + female[,2],
            HRmean_noL1 + male[,3], HRmean_noL1 + female[,3]),
  sex_level = rep(rep(c("Male", "Female"), each = 24), 6),
  hour = rep(12:35, 12),
  group = c(rep("DBP",24*2),rep("SBP",24*2),rep("HR",24*2)))
sex_data$group <- factor(sex_data$group, levels = c("DBP", "SBP", "HR"))

pdataL <- sex_data %>%
  ggplot(aes(x = hour, y = value, group = sex_level, col = sex_level)) +
  geom_line(alpha = 1) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(col="Sex") +
  #geom_line(aes(x = hour, y = mean), color="black") +
  facet_wrap(vars(group), scales="free_y") +
  theme(text = element_text(size = 15)) 
pdataL

# ggsave("./Application/sex_change_plot.pdf", dpi=600, width = 9, height=3)

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

# ggsave("./Application/effect_of_race.pdf", dpi=600, width = 8, height=3)

g1_white = -5.258 - 0.057*mean(xpred$age[xpred$race_White_inv==0]) + 0.106*mean(xpred$BMI[xpred$race_White_inv==0]) + 1.615*mean(xpred$odi4[xpred$race_White_inv==0]) - 0.275*mean(xpred$gender_inv[xpred$race_White_inv==0])
g2_white = 10.129 -0.110*mean(xpred$age[xpred$race_White_inv==0]) - 0.124*mean(xpred$BMI[xpred$race_White_inv==0]) + 1.384*mean(xpred$gender_inv[xpred$race_White_inv==0]) + 0.060*mean(xpred$odi4[xpred$race_White_inv==0])
g1_others = -5.258 - 0.057*mean(xpred$age[xpred$race_White_inv==1]) + 0.106*mean(xpred$BMI[xpred$race_White_inv==1]) + 1.948 + 1.615*mean(xpred$odi4[xpred$race_White_inv==1]) - 0.275*mean(xpred$gender_inv[xpred$race_White_inv==1])
g2_others = 10.129 -0.110*mean(xpred$age[xpred$race_White_inv==1]) - 0.124*mean(xpred$BMI[xpred$race_White_inv==1]) - 0.505 + 1.384*mean(xpred$gender_inv[xpred$race_White_inv==1]) + 0.060*mean(xpred$odi4[xpred$race_White_inv==1])

white <- g1_white * C1 + g2_white * C2
others <- g1_others * C1 + g2_others * C2

race_data = data.frame(
  #mean = c(rep(DBPmean,2), rep(SBPmean,2), rep(HRmean, 2)),
  value = c(DBPmean + white[,1], DBPmean + others[,1], 
            SBPmean + white[,2], SBPmean + others[,2],
            HRmean + white[,3], HRmean + others[,3]),
  race_level = rep(rep(c("White", "Others"), each = 24), 6),
  hour = rep(12:35, 12),
  group = c(rep("DBP",24*2),rep("SBP",24*2),rep("HR",24*2)))
race_data$group <- factor(race_data$group, levels = c("DBP", "SBP", "HR"))

## Estimation for White vs others
pdataL <- race_data %>%
  ggplot(aes(x = hour, y = value, group = race_level, col = race_level)) +
  geom_line(alpha = 1) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(col="Race") +
  #geom_line(aes(x = hour, y = mean), color="black") +
  facet_wrap(vars(group), scales="free_y") +
  theme(text = element_text(size = 15)) 
pdataL

#ggsave("./Application/race_change_plot.pdf", dpi=600, width = 9, height=3)

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

# ggsave("./Application/effect_of_odi4.pdf", dpi=600, width = 8, height=3)

DBPmean_odi4_low =  DBPmean - effectEst_all[[5]][,1]*1
DBPmean_odi4_up =  DBPmean + effectEst_all[[5]][,1]*1

SBPmean_odi4_low =  SBPmean - effectEst_all[[5]][,2]*1
SBPmean_odi4_up =  SBPmean + effectEst_all[[5]][,2]*1

HRmean_odi4_low = HRmean - effectEst_all[[5]][,3]*1
HRmean_odi4_up = HRmean + effectEst_all[[5]][,3]*1

odi4_data = data.frame(
  mean = c(rep(DBPmean,2), rep(SBPmean,2), rep(HRmean, 2)),
  value = c(DBPmean_odi4_low, DBPmean_odi4_up, 
            SBPmean_odi4_low, SBPmean_odi4_up, 
            HRmean_odi4_low, HRmean_odi4_up),
  odi4_level = rep(rep(c("-1", "+1"), each = 24), 6),
  hour = rep(12:35, 12),
  group = c(rep("DBP",24*2),rep("SBP",24*2),rep("HR",24*2)))
odi4_data$group <- factor(odi4_data$group, levels = c("DBP", "SBP", "HR"))

pdataL <- odi4_data %>%
  ggplot(aes(x = hour, y = value, group = odi4_level, col= odi4_level)) +
  geom_line(alpha = 1) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(color = "log(ODI4)\nChange") +
  geom_line(aes(x = hour, y = mean), color="black") +
  facet_wrap(vars(group), scales="free_y") +
  theme(text = element_text(size = 15)) 
pdataL

# ggsave("./Application/odi4_change_plot.pdf", dpi=600, width = 8, height=3)

## Estimation for median ODI4 + others at mean levels
g1_median = -5.258 - 0.057*mean(xpred$age) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*median(xpred$odi4)
g2_median = 10.129 -0.110*mean(xpred$age) - 0.124*mean(xpred$BMI) + 1.384*mean(xpred$gender_inv) - 0.505*mean(xpred$race_White_inv) + 0.060*median(xpred$odi4)

g1_q1 = -5.258 - 0.057*mean(xpred$age) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*quantile(xpred$odi4, 0.25)
g2_q1 = 10.129 -0.110*mean(xpred$age) - 0.124*mean(xpred$BMI) + 1.384*mean(xpred$gender_inv) - 0.505*mean(xpred$race_White_inv) + 0.060*quantile(xpred$odi4, 0.25)

g1_q3 = -5.258 - 0.057*mean(xpred$age) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*quantile(xpred$odi4, 0.75)
g2_q3 = -5.258 - 0.057*mean(xpred$age) + 0.106*mean(xpred$BMI) - 0.275*mean(xpred$gender_inv) + 1.948*mean(xpred$race_White_inv) + 1.615*quantile(xpred$odi4, 0.75)

odi4_median <- g1_median * C1 + g2_median * C2
odi4_q1 <- g1_q1 * C1 + g2_q1 * C2
odi4_q3 <- g1_q3 * C1 + g2_q3 * C2

DBPmean_noL1 = rep(70, 24) + (L_tilde[, 2] * mean(RiL2_191[1,]) + L_tilde[, 3] * mean(RiL3_191[1,])) * 12
SBPmean_noL1 = rep(129, 24) + (L_tilde[, 2] * mean(RiL2_191[2,]) + L_tilde[, 3] * mean(RiL3_191[2,])) * 18 # 18 is the sd of SBP
HRmean_noL1 = rep(76, 24) + (L_tilde[, 2] * mean(RiL2_191[3,]) + L_tilde[, 3] * mean(RiL3_191[3,])) * 13 # 13 is the sd of HR 

odi4_data = data.frame(
  value = c(DBPmean_noL1+odi4_median[,1], DBPmean_noL1+odi4_q1[,1], DBPmean_noL1+odi4_q3[,1], 
            SBPmean_noL1+odi4_median[,2], SBPmean_noL1+odi4_q1[,2], SBPmean_noL1+odi4_q3[,2], 
            HRmean_noL1+odi4_median[,3], HRmean_noL1+odi4_q1[,3], HRmean_noL1+odi4_q3[,3]),
  odi4_level = rep(rep(c("median","1st quartile", "3rd quartile"), each = 24), 3),
  hour = rep(12:35, 9),
  group = c(rep("DBP",24*3),rep("SBP",24*3),rep("HR",24*3)))
odi4_data$group <- factor(odi4_data$group, levels = c("DBP", "SBP", "HR"))
odi4_data$odi4_level <- factor(odi4_data$odi4_level, levels = c("1st quartile", "median", "3rd quartile"))

p <- odi4_data %>%
  ggplot(aes(x = hour, y = value, group = odi4_level, col = odi4_level)) +
  geom_line(alpha = 1) + 
  scale_color_manual(values = c("1st quartile" = "skyblue", "median" = "#4292C6", "3rd quartile" = "#2066A8"),
                     labels = c("1st quartile" = "Q1: 10.44", "median" = "median: 14.82", "3rd quartile" = "Q3: 22.93")) +
  scale_x_continuous(breaks = c(12, 15, 18, 21, 24, 27, 30, 33, 36), 
                     labels = c(12, 15, 18, 21, 0, 3, 6, 9, 12)) +
  xlab("Hour of the day") + 
  labs(color = "Quantile") +
  facet_wrap(vars(group), scales = "free_y") +
  theme(text = element_text(size = 15)) 
p

#ggsave("./Application/odi4_change_plot.pdf", dpi=600, width = 9, height=3)


