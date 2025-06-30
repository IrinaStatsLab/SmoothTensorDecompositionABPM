# load simulation results for Case 2-Flexible ranks
load("Simulation Studies/full/sim_res_miss0_full.Rda")
load("Simulation Studies/full/sim_res_miss10_full.Rda")
load("Simulation Studies/full/sim_res_miss20_full.Rda")
load("Simulation Studies/full/sim_res_miss50_full.Rda")
load("Simulation Studies/full/sim_res_missstruc_full.Rda")
load("Simulation Studies/full/sim_res_p50_full.Rda")
load("Simulation Studies/full/sim_res_p500_full.Rda")
load("Simulation Studies/full/sim_res_noise01_full.Rda")
load("Simulation Studies/full/sim_res_noise2_full.Rda")

library(ggplot2)
library(ggh4x)

Setting <- c("m0", "m10", "m20", "m50", "struc", "p50", "p500", "n01", "n2")
full_results_list <- list(m0 = results_miss0_full, m10 = results_miss10_full,
                     m20 = results_miss20_full, m50 = results_miss50_full,
                     struc = results_missstruc_full,
                     p50 = results_p50_full, p500 = results_p500_full,
                     n01 = results_noise01_full, n2 = results_noise2_full)

oracle_M_full <- list()
oracle_L_full <- list()
kcv_M_full <- list()
kcv_L_full <- list()
fpca_M_full <- list()
fpca_L_full <- list()
oracle_opt_para_full <- list()
kcv_opt_para_full <- list()
Lfpca_full <- list()

# Reorganize the simulation results into lists 
for (set in Setting) {
    oracle_M_full[[set]] <- oracle_L_full[[set]] <- 
    kcv_M_full[[set]] <- kcv_L_full[[set]] <- 
    fpca_M_full[[set]] <- fpca_L_full[[set]] <- 
    Lfpca_full[[set]] <- rep(NA, 100)
  
    oracle_opt_para_full[[set]] <- kcv_opt_para_full[[set]] <- matrix(NA, 100, 3)
  
    for (i in 1:100) {
      dat <- full_results_list[[set]][[i]]
      oracle_M_full[[set]][i] <- dat$oracle_loss$loss_M
      oracle_L_full[[set]][i] <- dat$oracle_loss$loss_L
      kcv_M_full[[set]][i] <- dat$kcv_loss$loss_M
      kcv_L_full[[set]][i] <- dat$kcv_loss$loss_L
      fpca_M_full[[set]][i] <- dat$fpca_loss$loss_M
      fpca_L_full[[set]][i] <- dat$fpca_loss$loss_L
      oracle_opt_para_full[[set]][i, ] <- dat$oracle_opt$opt_para
      kcv_opt_para_full[[set]][i, ] <- dat$kcv_opt$opt_para
      Lfpca_full[[set]][i] <- ncol(dat$fpca_loss$Lfpca)
    }
}

Setting_codes <- c("m0", "m10", "m20", "m50", "struc", "p50", "m20", "p500", "n01", "m20", "n2")
Setting_labels <- c("0%", "10%", "20%", "50%", "structured", 
                    "n=50", "n=200", "n=500", 
                    "Noise: 0.1", "Noise: 1", "Noise: 2")
panel_labels <- c(rep("Missingness", 5), rep("Sample size", 3), rep("Noise level", 3))

methods_full <- list(Oracle_Case2 = oracle_M_full, k_fold_cv_Case2 = kcv_M_full, FPCA_Case2 = fpca_M_full)

data_list_full <- list()

for (method in names(methods_full)) {
  for (i in seq_along(Setting_codes)) {
    data_list_full[[length(data_list_full) + 1]] <- data.frame(
      group = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = methods_full[[method]][[Setting_codes[i]]]
    )
  }
}

data_full <- do.call(rbind, data_list_full)

# load simulation results for Case 1-Fixed ranks
load("Simulation Studies/cf3/cf3_res_miss0.Rda")
load("Simulation Studies/cf3/cf3_res_miss10.Rda")
load("Simulation Studies/cf3/cf3_res_miss20.Rda")
load("Simulation Studies/cf3/cf3_res_miss50.Rda")
load("Simulation Studies/cf3/cf3_res_struc.Rda")
load("Simulation Studies/cf3/cf3_res_p50.Rda")
load("Simulation Studies/cf3/cf3_res_p500.Rda")
load("Simulation Studies/cf3/cf3_res_noise01.Rda")
load("Simulation Studies/cf3/cf3_res_noise2.Rda")

fixr_results_list <- list(m0 = fixr_miss0, m10 = fixr_miss10,
                          m20 = fixr_miss20, m50 = fixr_miss50,
                          struc = fixr_missstruc,
                          p50 = fixr_p50, p500 = fixr_p500,
                          n01 = fixr_noise01, n2 = fixr_noise2)

oracle_M_fixr <- list()
oracle_L_fixr <- list()
kcv_M_fixr <- list()
kcv_L_fixr <- list()
fpca_M_fixr <- list()
fpca_L_fixr <- list()

# Reorganize the simulation results into lists
for (set in Setting) {
  oracle_M_fixr[[set]] <- oracle_L_fixr[[set]] <- 
    kcv_M_fixr[[set]] <- kcv_L_fixr[[set]] <- 
    fpca_M_fixr[[set]] <- fpca_L_fixr[[set]]  <- rep(NA, 100)

  
  for (i in 1:100) {
    dat <- fixr_results_list[[set]][[i]]
    oracle_M_fixr[[set]][i] <- dat$oracle_loss$loss_M
    oracle_L_fixr[[set]][i] <- dat$oracle_loss$loss_L
    kcv_M_fixr[[set]][i] <- dat$kcv_loss$loss_M
    kcv_L_fixr[[set]][i] <- dat$kcv_loss$loss_L
    fpca_M_fixr[[set]][i] <- dat$fpca_loss$loss_M
    fpca_L_fixr[[set]][i] <- dat$fpca_loss$loss_L
  }
}

methods_fixr <- list(Oracle_Case1 = oracle_M_fixr, k_fold_cv_Case1 = kcv_M_fixr, FPCA_Case1 = fpca_M_fixr)

data_list_fixr <- list()

for (method in names(methods_fixr)) {
  for (i in seq_along(Setting_codes)) {
    data_list_fixr[[length(data_list_fixr) + 1]] <- data.frame(
      group = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = methods_fixr[[method]][[Setting_codes[i]]]
    )
  }
}

data_fixr <- do.call(rbind, data_list_fixr)

data <- data.frame(rbind(data_fixr, data_full))

data$group <- factor(data$group, levels=c("Oracle_Case1", "k_fold_cv_Case1", "FPCA_Case1",
                                          "Oracle_Case2", "k_fold_cv_Case2", "FPCA_Case2"))
data$Setting <- factor(data$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0.1","Noise: 1","Noise: 2"))
data$panel <- factor(data$panel, levels=c("Missingness", "Sample size", "Noise level"))

data$Case <- ifelse(data$group %in% c("Oracle_Case1", "k_fold_cv_Case1", "FPCA_Case1"), "Case 1", "Case 2")
data$Case <- factor(data$Case, levels = c("Case 1", "Case 2"))

# Make data for ggplot
panel_limits <- data %>%
  group_by(panel) %>%
  summarise(ymin = min(y_value, na.rm = TRUE),
            ymax = max(y_value, na.rm = TRUE))

data_plot <- data %>%
  left_join(panel_limits, by = "panel") 

data_plot$method_group[data_plot$group %in% c("Oracle_Case1","Oracle_Case2")] <- "Oracle"
data_plot$method_group[data_plot$group %in% c("k_fold_cv_Case1","k_fold_cv_Case2")] <- "SmoothHOOI"
data_plot$method_group[data_plot$group %in% c("FPCA_Case1","FPCA_Case2")] <- "FPCA"
data_plot$method_group <- factor(data_plot$method_group, levels=c("Oracle","SmoothHOOI","FPCA"))

# Plot of losses of M
ggplot(data_plot, aes(x = Setting, y = y_value, fill = method_group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  facet_nested(rows = NULL, cols = vars(panel, Case), scales = "free", independent = "y", space = "free_x") +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  labs(y = "Loss of M", fill = "Group", title = "") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(hjust = 0.5, size = 18),
    axis.text.y = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

#ggsave("./Simulation Studies/boxM.pdf", dpi=600, width=13, height = 6)


lossL <- list(Oracle = oracle_L_fixr, SmoothHOOI = kcv_L_fixr, FPCA = fpca_L_fixr)

data_list_lossL <- list()

for (method in names(lossL)) {
  for (i in seq_along(Setting_codes)) {
    data_list_lossL[[length(data_list_lossL) + 1]] <- data.frame(
      group = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = lossL[[method]][[Setting_codes[i]]]
    )
  }
}

data_lossL <- do.call(rbind, data_list_lossL)

data_lossL$group <- factor(data_lossL$group, levels=c("Oracle", "SmoothHOOI", "FPCA"))
data_lossL$Setting <- factor(data_lossL$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0.1","Noise: 1","Noise: 2"))
data_lossL$panel <- factor(data_lossL$panel, levels=c("Missingness", "Sample size", "Noise level"))

# Plot of losses of L
ggplot(data_lossL, aes(x = Setting, y = y_value, fill=group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  facet_nested(~panel, scales = "free", space ="free_x", independent = "y") +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  #labs(y = "Loss of L", fill = "Group", title = "Boxplots for Loss of L") +
  labs(y = "Loss of L", fill = "Group", title = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size=15),
        #units="px"plot.title = element_text(hjust=0.5, size=15),
        axis.title = element_text(hjust=0.5, size=18),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16),       
        legend.title = element_text(size = 16))

#ggsave("./Simulation Studies/boxL.pdf", dpi=600, width=11, height = 5)

# Plot selected ranks
my_selected_r <- data.frame(
  Value = rep(c(2, 3, 4, 5, 6), 11),
  Count = c(0, 98, 1, 0, 1, 
            0, 100, 0 , 0, 0,
            0, 92, 4, 4, 0,
            0, 50, 23, 16, 11,
            0, 61, 22, 6, 11,
            0, 65, 7, 13, 15, 
            0, 92, 4, 4, 0,
            0, 99, 1, 0, 0,
            0, 100, 0, 0, 0, 
            0, 92, 4, 4, 0,
            3, 32, 26, 21, 18),
  Case = rep(c("0%", "10%", "20%", "50%", "Structured", "n=50", "n=200","n=500", "Noise: 0.1", "Noise: 1","Noise: 2"), each = 5),
  Panel = c(rep("Missingness",5*5), rep("Sample size",3*5), rep("Noise level",3*5)),
  Method = rep("SmoothHOOI", 55)
)

my_selected_r$Case <- factor(my_selected_r$Case, levels=c("0%", "10%", "20%", "50%", "Structured", "n=50", "n=200","n=500", "Noise: 0.1", "Noise: 1", "Noise: 2"))
my_selected_r$Value <- factor(my_selected_r$Value, levels=c(6,5,4,3,2))
my_selected_r$Panel <- factor(my_selected_r$Panel, levels=c("Missingness", "Sample size", "Noise level"))


fpca_selected_r <- data.frame(
  Value = rep(c(2, 3, 4, 5, 6), 11),
  Count = c(0, 99, 1, 0, 0, 
            0, 90, 10 , 0, 0,
            0, 72, 28, 0, 0,
            0, 56, 44, 0, 0,
            0, 72, 28, 0, 0,
            0, 66, 34, 0, 0, 
            0, 72, 28, 0, 0,
            0, 92, 8, 0, 0,
            0, 100, 0, 0, 0, 
            0, 72, 28, 0, 0,
            0, 53, 46, 1, 0),
  Case = rep(c("0%", "10%", "20%", "50%", "Structured", "n=50", "n=200","n=500", "Noise: 0.1", "Noise: 1","Noise: 2"), each = 5),
  Panel = c(rep("Missingness",5*5), rep("Sample size",3*5), rep("Noise level",3*5)),
  Method = rep("FPCA", 55)
)

fpca_selected_r$Case <- factor(fpca_selected_r$Case, levels=c("0%", "10%", "20%", "50%", "Structured", "n=50", "n=200","n=500", "Noise: 0.1", "Noise: 1", "Noise: 2"))
fpca_selected_r$Value <- factor(fpca_selected_r$Value, levels=c(6,5,4,3,2))
fpca_selected_r$Panel <- factor(fpca_selected_r$Panel, levels=c("Missingness", "Sample size", "Noise level"))

selected_r <- data.frame(rbind(fpca_selected_r, my_selected_r))

selected_r$Method <- factor(selected_r$Method, levels=c("SmoothHOOI", "FPCA"))

custom_colors <- c("2" = "#CCCCCC", "3" = "#B3CDE3", "4" = "#CCEBC5", "5" = "#FED9A6", "6" = "#FDDAEC")

custom_strips <- strip_nested(
  background_x = elem_list_rect(
    fill = c("grey85", "transparent"), 
    color = c("transparent", "transparent")         
  ),
  text_x = elem_list_text(
    face = c("plain", "plain"), 
    size = 13
  ),
  by_layer_x = TRUE,
  by_layer_y = TRUE 
)

ggplot(selected_r, aes(x = Case, y = Count, fill = Value)) +
  scale_fill_manual(values = custom_colors, name = "r1 values") +
  geom_bar(stat = "identity") + 
  facet_nested(rows = NULL, cols = vars(Panel, factor(Method)), scales = "free_x", strip = custom_strips, space="free_x") +
  labs(x = "Setting", y = "Counts", fill = "r1 values", 
       title = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=30, hjust = 1, vjust = 1, size=15),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15)) +
  guides(fill = guide_legend(reverse = TRUE)) 

# ggsave("./Simulation Studies/rankselect.pdf", dpi=600, width = 12, height=6)

