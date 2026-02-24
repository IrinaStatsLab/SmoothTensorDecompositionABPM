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
library(patchwork)
library(ggh4x)
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

# Scaling
sdD = c(12, 18, 13) 

# Main effects for each model 
C1 = (L_tilde[, 1] %*% t(R_tilde[, 1])) %*% diag(sdD)
C2 = (L_tilde[, 1] %*% t(R_tilde[, 2])) %*% diag(sdD)

# mean level from L2 and L3
DBPmean_noL1 = rep(70, 24) + (L_tilde[, 2] * mean(RiL2_191[1,]) + L_tilde[, 3] * mean(RiL3_191[1,])) * 12
SBPmean_noL1 = rep(129, 24) + (L_tilde[, 2] * mean(RiL2_191[2,]) + L_tilde[, 3] * mean(RiL3_191[2,])) * 18 
HRmean_noL1 = rep(76, 24) + (L_tilde[, 2] * mean(RiL2_191[3,]) + L_tilde[, 3] * mean(RiL3_191[3,])) * 13 


models <- list(out1, out2)
vars <- names(coef(out1))

# Extract model coefficients
B <- sapply(models, function(m) coef(m)[vars])

# Calculate ODI quantiles
odi_quantiles <- quantile(xpred$odi4, probs = c(0.25, 0.5, 0.75))
names(odi_quantiles) <- c("1st quartile", "median", "3rd quartile")

generate_predictions <- function(gender_val, race_val, odi_level, odi_val) {
  
  x_vec <- setNames(rep(NA, length(vars)), vars)
  x_vec[c("(Intercept)", "age", "BMI", "gender_inv", "race_White_inv", "odi4")] <-
    c(1, mean(xpred$age), mean(xpred$BMI), gender_val, race_val, odi_val)

  g <- drop(x_vec %*% B)
  q <- g[1] * C1 + g[2] * C2
  
  data.frame(
    hour = 12:35,
    DBP = DBPmean_noL1 + q[, 1],
    SBP = SBPmean_noL1 + q[, 2],
    HR = HRmean_noL1 + q[, 3],
    gender = ifelse(gender_val == 0, "Male", "Female"),
    race = ifelse(race_val == 0, "White", "Other"),
    odi4_level = odi_level
  )
}

# Use expand_grid to create all combinations of inputs we need
param_grid <- expand.grid(gender_val = c(0, 1), race_val = c(0, 1), odi4_level = names(odi_quantiles))

results_list <- list()

for (i in 1:nrow(param_grid)) {
  current_params <- param_grid[i, ]
  odi_value <- odi_quantiles[current_params$odi4_level]
  
  temp_df <- generate_predictions(
    gender_val = current_params$gender_val,
    race_val = current_params$race_val,
    odi_level = current_params$odi4_level,
    odi_val = odi_value
  )
  
  results_list[[i]] <- temp_df
}

results_list

all_data_inverted <- bind_rows(results_list) %>% 
  pivot_longer(
    cols = c(DBP, SBP, HR),
    names_to = "group",
    values_to = "value"
  ) %>%
  mutate(
    # Create a new combined demographic variable
    demographic = factor(paste(gender, race, sep = ", "), 
                         levels = c("Male, White", "Male, Other", "Female, White", "Female, Other")),
    group = factor(group, levels = c("DBP", "SBP", "HR")),
    odi4_level = factor(odi4_level, levels = c("1st quartile", "median", "3rd quartile"))
  )

all_data_inverted

# Define a consistent theme for all three main plots
plot_theme_inverted <- list(
  geom_line(alpha = 1, linewidth = 0.8),
  facet_wrap(~ demographic, nrow = 1), 
  scale_color_manual(
    name = "ODI4 Quantile",
    values = c("1st quartile" = "skyblue", "median" = "#4292C6", "3rd quartile" = "#2066A8"),
    labels = c(
      "1st quartile" = paste("Q1:", round(exp(odi_quantiles[1]),2)),
      "median" = paste("Median:", round(exp(odi_quantiles[2]),2)),
      "3rd quartile" = paste("Q3:", round(exp(odi_quantiles[3]),2))
    )
  ),
  scale_x_continuous(
    breaks = c(12,  18,  24,  30,  36),
    labels = c(12,  18,  0, 6, 12)
  ),
  theme_minimal(base_size = 14),
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
)


p_dbp <- all_data_inverted %>%
  filter(group == "DBP") %>%
  ggplot(aes(x = hour, y = value, group = odi4_level, color = odi4_level)) +
  plot_theme_inverted +
  labs(title = "Diastolic Blood Pressure (DBP)", x = NULL, y = "value") +
  theme(axis.text.y = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(hjust=0.5, size=16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),       
        legend.title = element_text(size = 16))


p_sbp <- all_data_inverted %>%
  filter(group == "SBP") %>%
  ggplot(aes(x = hour, y = value, group = odi4_level, color = odi4_level)) +
  plot_theme_inverted +
  labs(title = "Systolic Blood Pressure (SBP)", x = NULL, y = "value") +
  theme(axis.text.y = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(hjust=0.5, size=16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),       
        legend.title = element_text(size = 16))


p_hr <- all_data_inverted %>%
  filter(group == "HR") %>%
  ggplot(aes(x = hour, y = value, group = odi4_level, color = odi4_level)) +
  plot_theme_inverted +
  labs(title = "Heart Rate (HR)", x = "Hour of the Day", y = "value") +
  theme(axis.text.y = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(hjust=0.5, size=16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),       
        legend.title = element_text(size = 16))


# Stack the three plots vertically
p_dbp / p_sbp / p_hr +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# p_dbp / p_hr +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom")

#ggsave("odi4_sex_race.pdf", dpi=600, width = 10, height=8)

