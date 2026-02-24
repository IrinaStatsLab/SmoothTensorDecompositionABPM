load("./AlgorithmResIdent.Rda") # L_tilde, R_tilde, G_tilde
load("./filtered_ABPM_new.Rda") # filtered_abpm_new
load("./combinedCGMmetricsCovariates.Rda") # combined

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

## check the orders of the ids of the tensor and the covariates correspond to each other
# unique(filtered_abpm_new$ID)[matched_idx]==covariates$id

## check the distribution of the covariates
# pdf("covariates_plots.pdf", width = 8, height = 6)
# 
# # Iterate through each covariates
# for (i in 2:34) {
#   hist(xpred[[i]], main=colnames(xpred)[i])
# }
# 
# # Close the PDF device
# dev.off()

## Make race=white and gender=male as reference levels
xpred$race_White_inv <- 1 - xpred$race_White 
xpred$gender_inv <- 1 - xpred$gender # in gender, male = 1. So we create a new variable gender_inv where female = 1

## Check outliers
summary(xpred$odi4)
boxplot(xpred$odi4)

summary(xpred$age)
boxplot(xpred$age)

## Summary Tables

summary_table <- xpred %>%
  mutate(
    gender = factor(gender, labels = c("Female","Male")),
    race_White = factor(race_White, labels = c("Others", "White")),
    exp_odi4 = exp(odi4)
  ) %>%
  select(age, gender, race_White, BMI, exp_odi4) %>% 
  tbl_summary(
    by = gender,         
    statistic = list(all_continuous() ~ "{mean} ({sd}); {median}[{p25}, {p75}]",
                     all_categorical() ~ "{n} ({p}%)"),  
    digits = all_continuous() ~ 1
  ) 


# Regression Model 1
scoredata1  <- cbind(G_tilde[1, 1, matched_idx], xpred)
names(scoredata1)[1] = "TensorScore1"
out1 = lm(TensorScore1 ~ age + BMI + gender_inv + race_White_inv + odi4, data = scoredata1)
summary(out1) 

confint(out1)

# Regression Model 2
scoredata2 = cbind(G_tilde[1, 2, matched_idx], xpred)
names(scoredata2)[1] = "TensorScore2"
out2 = lm(TensorScore2 ~ age + BMI + gender_inv + race_White_inv + odi4, data = scoredata2)
summary(out2) 

confint(out2)

# Model Diagnostics
# plot(out1)
# plot(out2)
# vif(out1)
# vif(out2)

# Regression Model 3
scoredata3  <- cbind(G_tilde[2, 1, matched_idx], xpred)
names(scoredata3)[1] = "TensorScore3"
out3 = lm(TensorScore3 ~ age + BMI + gender_inv + race_White_inv + odi4, data = scoredata3)
summary(out3) 

confint(out3)

# Regression Model 4
scoredata4 = cbind(G_tilde[2, 2, matched_idx], xpred)
names(scoredata4)[1] = "TensorScore4"
out4 = lm(TensorScore4 ~ age + BMI + gender_inv + race_White_inv + odi4, data = scoredata4)
summary(out4) 

confint(out4)

# Bonferroni Adjustment for model 1 and 2
scoresdata = cbind(G_tilde[1, 1, matched_idx], G_tilde[1, 2, matched_idx], xpred)
names(scoresdata)[1] = "TensorScore1"
names(scoresdata)[2] = "TensorScore2"
out = lm(cbind(TensorScore1, TensorScore2) ~ age + BMI + gender_inv + race_White_inv + odi4,
         data = scoresdata)
Manova(out, test.statistic = "Pillai")
for (i in c(2:5, 0)){
  print(tidy(out) %>% 
          filter(row_number() %% 6 == i) %>% 
          mutate(p.value_adjusted=p.adjust(p.value, 
                                           method="bonferroni")))
}

# Regression Models for Day/Night means

## Take means of measurements for day and night
means_bynight = filtered_abpm_new %>%
  mutate(type = ifelse(TIME_OF_DAY == "Night", "Night", "Day")) %>%
  group_by(ID, type) %>%
  summarize(SBPmean = mean(SBP, na.rm = TRUE), DBPmean = mean(DBP, na.rm = TRUE),
            MAPmean = mean(MAP, na.rm = TRUE), HRmean = mean(HR, na.rm = TRUE), .groups = "drop")

## Visualize means_bynight

df_long <- means_bynight %>%
  select(ID, type, SBPmean, DBPmean, HRmean) %>%
  pivot_longer(
    cols = c(DBPmean, SBPmean, HRmean),
    names_to = "Measure",
    values_to = "Value"
  ) %>% 
  mutate(Measure = factor(Measure, levels = c("DBPmean", "SBPmean", "HRmean")))
  
ggplot(df_long, aes(x = type, y = Value)) +
  geom_boxplot(aes(fill = type),
               width = 0.35,
               outlier.shape = NA,
               alpha = 0.4) +
  geom_line(aes(group = ID), alpha = 0.2) +
  geom_point(size = 0.5) +
  facet_wrap(~ Measure, ncol = 3, scales = "free_y", 
             labeller = as_labeller(c(
               "DBPmean" = "DBP (mmHg)",
               "SBPmean" = "SBP (mmHg)",
               "HRmean" = "HR (bpm)"
             ))) +
  labs(x = "Time of day", y = "Value") +
  theme_bw() +
  theme(legend.position = "none", 
        strip.text = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

# ggsave("Day_Night_Boxplots.pdf", dpi=600, width=9, height=5)

## Correlation between Day and Night
wide_df <- df_long %>%
  pivot_wider(names_from = type, values_from = Value)

correlations <- wide_df %>%
  group_by(Measure) %>%
  summarize(
    correlation = cor(Day, Night, use = "complete.obs"),
    p_value = cor.test(Day, Night)$p.value
  )
correlations

##DBPmean       0.651 
##SBPmean       0.723 
##HRmean        0.731

## correlation between 6 measures
means_bynight_wide = pivot_wider(means_bynight, id_cols = "ID", names_from = "type", values_from = c("SBPmean", "DBPmean", "MAPmean", "HRmean"))
means_bynight_wide$ID = as.numeric(as.character(means_bynight_wide$ID))

cor_matrix <- cor(means_bynight_wide[, c("DBPmean_Day","SBPmean_Day", "HRmean_Day", 
                                         "DBPmean_Night","SBPmean_Night", "HRmean_Night")], 
                  use = "complete.obs")
cor_matrix

library(ggcorrplot)

new_names <- c("DBP_Day","SBP_Day","HR_Day","DBP_Night","SBP_Night","HR_Night")
colnames(cor_matrix) <- new_names
rownames(cor_matrix) <- new_names

ggcorrplot(cor_matrix, 
           method = "square",       
           type = "lower",          
           lab = TRUE,              
           lab_size = 5,           
           colors = c("blue", "white", "firebrick"), 
           # title = "Correlation Matrix of ABPM Parameters",
           ggtheme = theme_minimal()) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), 
    #axis.text.y = element_text(size = 12),
    axis.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
  )

#ggsave("Day_Night_Measure_Corr.pdf", dpi=600, width=7, height=7)

alldata = left_join(xpred, means_bynight_wide, by= c("id" = "ID"))

outDBPDay = lm(DBPmean_Day ~ age + BMI + gender_inv + race_White_inv + odi4, data = alldata)
summary(outDBPDay)

outDBPNight = lm(DBPmean_Night ~ age + BMI + gender_inv + race_White_inv + odi4, data = alldata)
summary(outDBPNight)

outSBPDay = lm(SBPmean_Day ~ age + BMI + gender_inv + race_White_inv + odi4, data = alldata)
summary(outSBPDay)

outSBPNight = lm(SBPmean_Night ~ age + BMI + gender_inv + race_White_inv + odi4, data = alldata)
summary(outSBPNight)

outHRDay = lm(HRmean_Day ~ age + BMI + gender_inv + race_White_inv + odi4, data = alldata)
summary(outHRDay)

outHRNight = lm(HRmean_Night ~ age + BMI + gender_inv + race_White_inv + odi4, data = alldata)
summary(outHRNight)

## Bonferroni Adjustment
out = lm(cbind(DBPmean_Day, DBPmean_Night, SBPmean_Day, SBPmean_Night,  HRmean_Day, HRmean_Night) ~ age + BMI + gender_inv + race_White_inv + odi4,
         data = alldata)
Manova(out, test.statistic = "Pillai")
for (i in c(2:5, 0)){
  print(tidy(out) %>% 
          filter(row_number() %% 6 == i) %>% 
          mutate(p.value_adjusted=p.adjust(p.value, 
                                           method="bonferroni")))
}

