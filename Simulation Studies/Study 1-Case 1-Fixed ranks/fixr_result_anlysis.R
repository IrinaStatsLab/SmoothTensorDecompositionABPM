library(ggplot2)
library(ggh4x)
library(dplyr)

Setting <- c("m0", "m10", "m20", "m50", "struc", "p50", "p500", "n0", "n01", "n05", "n15", "n2")
Setting_codes <- c("m0", "m10", "m20", "m50", "struc", "p50", "m20", "p500", "n0", "n01","n05", "m20", "n15", "n2")
Setting_labels <- c("0%", "10%", "20%", "50%", "structured", 
                    "n=50", "n=200", "n=500", 
                    "Noise: 0", "Noise: 0.1", "Noise: 0.5",
                    "Noise: 1", "Noise: 1.5", "Noise: 2")

panel_labels <- c(rep("Missingness", 5), rep("Sample size", 3), rep("Noise level", 6))

rda_files <- list.files("./fixr_new", pattern = "\\.Rda$", full.names = TRUE)
for (f in rda_files) {
  load(f)
}

fixr_results_list <- list(m0 = fixr_miss0, m10 = fixr_miss10,
                          m20 = fixr_miss20, m50 = fixr_miss50,
                          struc = fixr_missstruc,
                          p50 = fixr_p50, p500 = fixr_p500,
                          n0 = fixr_noise0, n01 = fixr_noise01,
                          n05 = fixr_noise05, n15 = fixr_noise15,
                          n2 = fixr_noise2)

oracle_M_fixr <- list()
oracle_L_fixr <- list()
kcv_M_fixr <- list()
kcv_L_fixr <- list()
fpca_M_fixr <- list()
fpca_L_fixr <- list()
cp3_M_fixr <- list()
cp3_L_fixr <- list()
cp6_M_fixr <- list()
cp6_L_fixr <- list()
cp3_no_M_fixr <- list()
cp3_no_L_fixr <- list()
cp6_no_M_fixr <- list()
cp6_no_L_fixr <- list()
mfpca_M_fixr <- list()
mfpca_L_fixr <- list()

# Reorganize the simulation results into lists
for (set in Setting) {
  oracle_M_fixr[[set]] <- oracle_L_fixr[[set]] <- 
    kcv_M_fixr[[set]] <- kcv_L_fixr[[set]] <- 
    fpca_M_fixr[[set]] <- fpca_L_fixr[[set]]  <- 
    cp3_M_fixr[[set]] <- cp3_L_fixr[[set]]  <- 
    cp6_M_fixr[[set]] <- cp6_L_fixr[[set]]  <-
    cp3_no_M_fixr[[set]] <- cp3_no_L_fixr[[set]]  <- 
    cp6_no_M_fixr[[set]] <- cp6_no_L_fixr[[set]]  <-
    mfpca_M_fixr[[set]] <- mfpca_L_fixr[[set]]  <- 
    rep(NA, 100)
  
  
  for (i in 1:100) {
    dat <- fixr_results_list[[set]][[i]]
    oracle_M_fixr[[set]][i] <- dat$oracle_loss$loss_M
    oracle_L_fixr[[set]][i] <- dat$oracle_loss$loss_L
    kcv_M_fixr[[set]][i] <- dat$kcv_loss$loss_M
    kcv_L_fixr[[set]][i] <- dat$kcv_loss$loss_L
    fpca_M_fixr[[set]][i] <- dat$fpca_loss$loss_M
    fpca_L_fixr[[set]][i] <- dat$fpca_loss$loss_L
    cp3_M_fixr[[set]][i] <- dat$cp3_loss$loss_M
    cp3_L_fixr[[set]][i] <- dat$cp3_loss$loss_L
    cp6_M_fixr[[set]][i] <- dat$cp6_loss$loss_M
    cp6_L_fixr[[set]][i] <- dat$cp6_loss$loss_L
    cp3_no_M_fixr[[set]][i] <- dat$cp3_no_loss$loss_M
    cp3_no_L_fixr[[set]][i] <- dat$cp3_no_loss$loss_L
    cp6_no_M_fixr[[set]][i] <- dat$cp6_no_loss$loss_M
    cp6_no_L_fixr[[set]][i] <- dat$cp6_no_loss$loss_L
    mfpca_M_fixr[[set]][i] <- dat$mfpca_loss$loss_M
    mfpca_L_fixr[[set]][i] <- dat$mfpca_loss$loss_L
    
  }
}

methods_fixr <- list(Oracle = oracle_M_fixr, SmoothHOOI = kcv_M_fixr, FPCA = fpca_M_fixr, CP3_Ortsmo = cp3_M_fixr, CP6_Ortsmo = cp6_M_fixr, CP3_Smooth = cp3_no_M_fixr, CP6_Smooth = cp6_no_M_fixr, MFPCA = mfpca_M_fixr)

data_list_fixr <- list()

for (method in names(methods_fixr)) {
  for (i in seq_along(Setting_codes)) {
    data_list_fixr[[length(data_list_fixr) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = methods_fixr[[method]][[Setting_codes[i]]]
    )
  }
}

data_fixr <- do.call(rbind, data_list_fixr)

data <- data.frame(data_fixr)

data$method <- factor(data$method, levels=c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP3_Ortsmo", "CP6_Ortsmo", "CP3_Smooth", "CP6_Smooth"))
data$Setting <- factor(data$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0", "Noise: 0.1", "Noise: 0.5","Noise: 1", "Noise: 1.5", "Noise: 2"))
data$panel <- factor(data$panel, levels=c("Missingness", "Sample size", "Noise level"))


panel_limits <- data %>%
  group_by(panel) %>%
  summarise(ymin = min(y_value, na.rm = TRUE),
            ymax = max(y_value, na.rm = TRUE))


data <- data %>%
  left_join(panel_limits, by = "panel") 

ggplot(data[which(data$Setting!="Noise: 0.1"),], aes(x = Setting, y = y_value, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  facet_nested(~panel, scales = "free", space ="free_x", independent = "y") +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  labs(y = "Loss of M", fill = "Method", title = "") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
    plot.title = element_text(hjust = 0.5, size = 13),
    axis.title = element_text(hjust = 0.5, size = 13),
    axis.text.y = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )
#ggsave("./fixr_lossM_boxplot_complete.pdf", width=12, height=10)

lossL <- list(Oracle = oracle_L_fixr, SmoothHOOI = kcv_L_fixr, FPCA = fpca_L_fixr, MFPCA = mfpca_L_fixr, CP3_Ortsmo = cp3_L_fixr, CP6_Ortsmo = cp6_L_fixr, CP3_Smooth = cp3_no_L_fixr, CP6_Smooth = cp6_no_L_fixr)

data_list_lossL <- list()

for (method in names(lossL)) {
  for (i in seq_along(Setting_codes)) {
    data_list_lossL[[length(data_list_lossL) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = lossL[[method]][[Setting_codes[i]]]
    )
  }
}

data_lossL <- do.call(rbind, data_list_lossL)

data_lossL$method <- factor(data_lossL$method, levels=c("Oracle", "SmoothHOOI","FPCA","MFPCA", "CP3_Ortsmo", "CP6_Ortsmo", "CP3_Smooth", "CP6_Smooth"))
data_lossL$Setting <- factor(data_lossL$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0", "Noise: 0.1", "Noise: 0.5","Noise: 1", "Noise: 1.5", "Noise: 2"))
data_lossL$panel <- factor(data_lossL$panel, levels=c("Missingness", "Sample size", "Noise level"))

# Plot of losses of L
ggplot(data_lossL[which(data_lossL$Setting!="Noise: 0.1"),], aes(x = Setting, y = y_value, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  labs(y = "Loss of L", fill = "Method", title = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size=15),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15),
        strip.placement = "outside") +
  facet_nested(~panel, scales = "free", space ="free_x", independent = "y") 

#ggsave("./fixr_lossL_boxplot_complete.pdf", width=15, height=6)

missing_rate_vec_fixr <- c()
for (i in 1:100){
  missing_rate_vec_fixr[i] <- fixr_missstruc[[i]]$missing_rate
}
summary(missing_rate_vec_fixr)


data_subset <- data.frame(data_fixr) 
data_subset <- subset(data_subset, method %in% c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP3_Smooth", "CP6_Smooth"))

data_subset$method <- factor(data_subset$method, levels=c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP3_Smooth", "CP6_Smooth"))
data_subset$Setting <- factor(data_subset$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0", "Noise: 0.1", "Noise: 0.5","Noise: 1", "Noise: 1.5", "Noise: 2"))
data_subset$panel <- factor(data_subset$panel, levels=c("Missingness", "Sample size", "Noise level"))

panel_limits <- data_subset %>%
  group_by(panel) %>%
  summarise(ymin = min(y_value, na.rm = TRUE),
            ymax = max(y_value, na.rm = TRUE))


data_subset <- data_subset %>%
  left_join(panel_limits, by = "panel") 

ggplot(data_subset[which(data_subset$Setting!="Noise: 0.1"),], aes(x = Setting, y = y_value, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  facet_nested(~panel, scales = "free", space ="free_x", independent = "y") +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  labs(y = "Loss of M", fill = "Method", title = "") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(hjust = 0.5, size = 18),
    axis.text.y = element_text(hjust = 0.5, size = 15),
    strip.text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
#ggsave("./fixr_lossM_boxplot_noOrtho.pdf", width=12, height=8)

lossL <- list(Oracle = oracle_L_fixr, SmoothHOOI = kcv_L_fixr, FPCA = fpca_L_fixr, MFPCA = mfpca_L_fixr, CP3_Smooth = cp3_no_L_fixr, CP6_Smooth = cp6_no_L_fixr)

data_list_lossL <- list()

for (method in names(lossL)) {
  for (i in seq_along(Setting_codes)) {
    data_list_lossL[[length(data_list_lossL) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = lossL[[method]][[Setting_codes[i]]]
    )
  }
}

data_lossL <- do.call(rbind, data_list_lossL)

data_lossL$method <- factor(data_lossL$method, levels=c("Oracle", "SmoothHOOI","FPCA","MFPCA", "CP3_Smooth", "CP6_Smooth"))
data_lossL$Setting <- factor(data_lossL$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0", "Noise: 0.1", "Noise: 0.5","Noise: 1", "Noise: 1.5", "Noise: 2"))
data_lossL$panel <- factor(data_lossL$panel, levels=c("Missingness", "Sample size", "Noise level"))


# Plot of losses of L
ggplot(data_lossL[which(data_lossL$Setting!="Noise: 0.1"),], aes(x = Setting, y = y_value, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  labs(y = "Loss of L", fill = "Method", title = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size=15),
        axis.title = element_text(hjust=0.5, size=18),
        axis.text.y = element_text(hjust=0.5, size=15),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16),       
        legend.title = element_text(size = 16),
        strip.placement = "outside") +
  facet_nested(~panel, scales = "free", space ="free_x", independent = "y") 

#ggsave("./fixr_lossL_boxplot_noOrtho.pdf", width=12, height=6)

cp3_ortho_L1_fixr <- list()
cp3_ortho_L2_fixr <- list()
cp3_ortho_L3_fixr <- list()
cp6_ortho_L1_fixr <- list()
cp6_ortho_L2_fixr <- list()
cp6_ortho_L3_fixr <- list()
cp3_ortho_L_fixr <- list()
cp6_ortho_L_fixr <- list()

cp3_L1_fixr <- list()
cp3_L2_fixr <- list()
cp3_L3_fixr <- list()
cp6_L1_fixr <- list()
cp6_L2_fixr <- list()
cp6_L3_fixr <- list()
cp3_L_fixr <- list()
cp6_L_fixr <- list()


# Reorganize the simulation results into lists
for (set in Setting) {
  cp3_ortho_L1_fixr[[set]] <- cp3_ortho_L2_fixr[[set]] <- cp3_ortho_L3_fixr[[set]] <- 
    cp6_ortho_L1_fixr[[set]] <- cp6_ortho_L2_fixr[[set]] <- cp6_ortho_L3_fixr[[set]]  <- 
    cp3_L1_fixr[[set]] <- cp3_L2_fixr[[set]] <- cp3_L3_fixr[[set]] <- 
    cp6_L1_fixr[[set]] <- cp6_L2_fixr[[set]] <- cp6_L3_fixr[[set]]  <-
    cp3_ortho_L_fixr[[set]] <- cp6_ortho_L_fixr[[set]]  <- 
    cp3_L_fixr[[set]] <- cp6_L_fixr[[set]]  <-
    rep(NA, 100)
  
  
  for (i in 1:100) {
    dat <- fixr_results_list[[set]][[i]]
    cp3_ortho_L1_fixr[[set]][i] <- dat$cp3_loss1$loss_L
    cp3_ortho_L2_fixr[[set]][i] <- dat$cp3_loss2$loss_L
    cp3_ortho_L3_fixr[[set]][i] <- dat$cp3_loss3$loss_L
    cp6_ortho_L1_fixr[[set]][i] <- dat$cp6_loss1$loss_L
    cp6_ortho_L2_fixr[[set]][i] <- dat$cp6_loss2$loss_L
    cp6_ortho_L3_fixr[[set]][i] <- dat$cp6_loss3$loss_L
    cp3_ortho_L_fixr[[set]][i] <- dat$cp3_loss$loss_L
    cp6_ortho_L_fixr[[set]][i] <- dat$cp6_loss$loss_L
    
    cp3_L1_fixr[[set]][i] <- dat$cp3_no_loss1$loss_L
    cp3_L2_fixr[[set]][i] <- dat$cp3_no_loss2$loss_L
    cp3_L3_fixr[[set]][i] <- dat$cp3_no_loss3$loss_L
    cp6_L1_fixr[[set]][i] <- dat$cp6_no_loss1$loss_L
    cp6_L2_fixr[[set]][i] <- dat$cp6_no_loss2$loss_L
    cp6_L3_fixr[[set]][i] <- dat$cp6_no_loss3$loss_L
    cp3_L_fixr[[set]][i] <- dat$cp3_no_loss$loss_L
    cp6_L_fixr[[set]][i] <- dat$cp6_no_loss$loss_L
    
  }
}

lossL <- list(CP3_Ortho_L1 = cp3_ortho_L1_fixr, Cp3_Ortho_L2 = cp3_ortho_L2_fixr, CP3_Ortho_L3 = cp3_ortho_L3_fixr, 
              CP3_Ortho = cp3_ortho_L_fixr,
              CP3_L1 = cp3_L1_fixr, CP3_L2 = cp3_L2_fixr, CP3_L3 = cp3_L3_fixr, CP3 = cp3_L_fixr)

data_list_lossL <- list()

for (method in names(lossL)) {
  for (i in seq_along(Setting_codes)) {
    data_list_lossL[[length(data_list_lossL) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = lossL[[method]][[Setting_codes[i]]]
    )
  }
}

data_lossL <- do.call(rbind, data_list_lossL)

data_lossL$method <- factor(data_lossL$method, levels=c("CP3_Ortho_L1", "Cp3_Ortho_L2", "CP3_Ortho_L3",
                                                        "CP3_L1", "CP3_L2", "CP3_L3","CP3_Ortho", "CP3"))
data_lossL$Setting <- factor(data_lossL$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0", "Noise: 0.1", "Noise: 0.5","Noise: 1", "Noise: 1.5", "Noise: 2"))
data_lossL$panel <- factor(data_lossL$panel, levels=c("Missingness", "Sample size", "Noise level"))

lossL <- list(
  # CP3
  Ortsmo_L1 = cp3_ortho_L1_fixr, 
  Ortsmo_L2 = cp3_ortho_L2_fixr, 
  Ortsmo_L3 = cp3_ortho_L3_fixr, 
  Ortsmo_L    = cp3_ortho_L_fixr,
  Smooth_L1       = cp3_L1_fixr, 
  Smooth_L2       = cp3_L2_fixr, 
  Smooth_L3       = cp3_L3_fixr, 
  Smooth_L          = cp3_L_fixr
)

data_list_lossL <- list()
for (method in names(lossL)) {
  for (i in seq_along(Setting_codes)) {
    data_list_lossL[[length(data_list_lossL) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = lossL[[method]][[Setting_codes[i]]]
    )
  }
}
data_lossL <- do.call(rbind, data_list_lossL)

data_lossL$Setting <- factor(data_lossL$Setting, levels=c("0%", "10%", "20%", "50%", "structured","n=50", "n=200","n=500","Noise: 0", "Noise: 0.1", "Noise: 0.5","Noise: 1", "Noise: 1.5", "Noise: 2"))
data_lossL$panel <- factor(data_lossL$panel, levels=c("Missingness", "Sample size", "Noise level"))

data_lossL$methodpanel  <- ifelse(grepl("Ort", data_lossL$method), "Ortsmo", "Smooth")
data_lossL$Lpanel  <- ifelse(grepl("L1", data_lossL$method), "L1",
                             ifelse(grepl("L2", data_lossL$method), "L2",
                                    ifelse(grepl("L3", data_lossL$method), "L3", "L")))
data_lossL$Lpanel <- factor(data_lossL$Lpanel, levels = c("L1", "L2", "L3", "L"))

data_lossL$method <- factor(data_lossL$method, levels = c(
  "Ortsmo_L1", "Smooth_L1", 
  "Ortsmo_L2", "Smooth_L2", 
  "Ortsmo_L3", "Smooth_L3",
  "Ortsmo_L", "Smooth_L"
))

# Define specific colors for your variable names
my_colors <- c(
  "Ortsmo_L1" = "#C6DBEF", "Ortsmo_L2" = "#6BAED6", "Ortsmo_L3" = "#2171B5", "Ortsmo_L" = "#08306B",
  "Smooth_L1"       = "#FCBBA1", "Smooth_L2"       = "#FB6A4A", "Smooth_L3"       = "#CB181D", "Smooth_L"       = "#67000D"
  
)

library(scales)

# Plot for CP3
ggplot(data_lossL[which(data_lossL$methodpanel=="Smooth" & data_lossL$Setting!="Noise: 0.1"),], aes(x = Setting, y = y_value, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, outlier.size = 0.5) +
  
  # Use the manual gradient colors defined above
  scale_fill_manual(values = my_colors) +
  
  facet_nested(~ panel + Lpanel, scales = "free", space = "free_x") +
  
  labs(y = "Loss of L", fill = "Method") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size=12),
        #units="px"plot.title = element_text(hjust=0.5, size=15),
        axis.title = element_text(hjust=0.5, size=15),
        axis.text.y = element_text(hjust=0.5, size=12),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),       
        legend.title = element_text(size = 15),
        strip.placement = "outside",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color = NA)) 
#ggsave("./CP3_singleL_vs_entireL.pdf", width=15, height=6)

library(ggplot2)
library(ggh4x)
library(scales)

# Plot for CP3
ggplot(data_lossL[which(data_lossL$methodpanel=="Ortsmo"& data_lossL$Setting!="Noise: 0.1"),], aes(x = Setting, y = y_value, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, outlier.size = 0.5) +
  
  # Use the manual gradient colors defined above
  scale_fill_manual(values = my_colors) +
  
  # facet_nested(~ panel + Subpanel, scales = "free", space = "free_x") +
  
  facet_nested(~ panel + Lpanel, scales = "free", space = "free_x") +
  
  labs(y = "Loss of L", fill = "Method") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size=10),
        #units="px"plot.title = element_text(hjust=0.5, size=15),
        axis.title = element_text(hjust=0.5, size=12),
        axis.text.y = element_text(hjust=0.5, size=12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),       
        legend.title = element_text(size = 12),
        strip.placement = "outside",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color = NA)) 
#ggsave("./CP3Ortho_singleL_vs_entireL.pdf", width=15, height=5)
