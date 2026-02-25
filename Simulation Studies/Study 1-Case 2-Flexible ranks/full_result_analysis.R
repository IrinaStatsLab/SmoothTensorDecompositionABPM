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

rda_files <- list.files("./full_new_v2", pattern = "\\.Rda$", full.names = TRUE)
for (f in rda_files) {
  load(f)
}

full_results_list <- list(m0 = full_miss0, m10 = full_miss10,
                          m20 = full_miss20, m50 = full_miss50,
                          struc = full_missstruc,
                          p50 = full_p50, p500 = full_p500,
                          n0 = full_noise0, n01 = full_noise01,
                          n05 = full_noise05, n15 = full_noise15,
                          n2 = full_noise2)

oracle_M_full <- list()
kcv_M_full <- list()
fpca_M_full <- list()
cp_M_full <- list()
cp_no_M_full <- list()
mfpca_M_full <- list()

oracle_rank_full <- list()
kcv_rank_full <- list()
fpca_rank_full <- list()
cp_rank_full <- list()
cp_no_rank_full <- list()
mfpca_rank_full <- list()

# Reorganize the simulation results into lists
for (set in Setting) {
  oracle_M_full[[set]] <- oracle_rank_full[[set]] <- 
    kcv_M_full[[set]] <- kcv_rank_full[[set]] <- 
    fpca_M_full[[set]] <- fpca_rank_full[[set]]  <- 
    cp_M_full[[set]] <- cp_rank_full[[set]]  <- 
    cp_no_M_full[[set]] <- cp_no_rank_full[[set]]  <-
    mfpca_M_full[[set]] <- mfpca_rank_full[[set]]  <- 
    rep(NA, 100)
  
  
  for (i in 1:100) {
    dat <- full_results_list[[set]][[i]]
    oracle_M_full[[set]][i] <- dat$oracle_loss$loss_M
    oracle_rank_full[[set]][i] <- dat$oracle_rank
    kcv_M_full[[set]][i] <- dat$kcv_loss$loss_M
    kcv_rank_full[[set]][i] <- dat$kcv_rank
    fpca_M_full[[set]][i] <- dat$fpca_loss$loss_M
    fpca_rank_full[[set]][i] <- dat$fpca_rank
    cp_M_full[[set]][i] <- dat$cp_loss$loss_M
    cp_rank_full[[set]][i] <- dat$cp_rank
    cp_no_M_full[[set]][i] <- dat$cp_no_loss$loss_M
    cp_no_rank_full[[set]][i] <- dat$cp_no_rank
    mfpca_M_full[[set]][i] <- dat$mfpca_loss$loss_M
    mfpca_rank_full[[set]][i] <- dat$mfpca_rank
    
  }
}

methods_full <- list(Oracle = oracle_M_full, SmoothHOOI = kcv_M_full, FPCA = fpca_M_full, CP_Ortsmo = cp_M_full, CP_Smooth = cp_no_M_full, MFPCA = mfpca_M_full)

data_list_full <- list()

for (method in names(methods_full)) {
  for (i in seq_along(Setting_codes)) {
    data_list_full[[length(data_list_full) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = methods_full[[method]][[Setting_codes[i]]]
    )
  }
}

data_full <- do.call(rbind, data_list_full)

data <- data.frame(data_full)

data$method <- factor(data$method, levels=c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP_Ortsmo", "CP_Smooth"))
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
#ggsave("./full_lossM_boxplot_complete.pdf", width=12, height=10)

data_subset <- data.frame(data_full) 
data_subset <- subset(data_subset, method %in% c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP_Smooth"))

data_subset$method <- factor(data_subset$method, levels=c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP_Smooth"))
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
    axis.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(hjust = 0.5, size = 15),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )
#ggsave("./full_lossM_boxplot_noOrtho.pdf", width=12, height=6)

# table() rank for all the methods in c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP_Ortho", "CP")
rank_tables <- list()
for (method in c("oracle_rank_full", "kcv_rank_full", "fpca_rank_full", "mfpca_rank_full", "cp_rank_full", "cp_no_rank_full")) {
  rank_tables[[method]] <- list()
  for (set in Setting) {
    rank_tables[[method]][[set]] <- table(get(method)[[set]]) 
  }
} 

# organize rank_tables into a data frame for plotting
rank_data_list <- list()
for (method in names(rank_tables)) {
  for (i in seq_along(Setting_codes)) {
    ranks <- names(rank_tables[[method]][[Setting_codes[i]]])
    counts <- as.numeric(rank_tables[[method]][[Setting_codes[i]]])
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      Rank = ranks,
      Count = counts
    )
  }
}

rank_data_list_full <- do.call(rbind, rank_data_list)
rank_data_list_full

rank_data_list_full <- rank_data_list_full %>% 
  filter(method %in% c("kcv_rank_full", "fpca_rank_full", "cp_rank_full")) %>% 
  filter(Setting %in% setdiff(unique(rank_data_list_full$Setting), c("Noise: 0.1"))) %>% 
  mutate(panel=factor(panel, levels=c("Missingness", "Sample size", "Noise level")))

rank_data_list_full$method <- dplyr::recode(rank_data_list_full$method,
                                            kcv_rank_full = "SmoothHOOI",
                                            fpca_rank_full = "FPCA",
                                            cp_rank_full = "CP_Ortsmo")

rank_data_list_full$method <- factor(rank_data_list_full$method, levels=c("SmoothHOOI", "FPCA", "CP_Ortsmo"))

rank_data_list_full$Rank <- factor(rank_data_list_full$Rank,
                                   levels = c("2","3","4","5","6"))

pref <- c("n = 50", "n = 200", "n = 500")
others <- setdiff(unique(rank_data_list_full$Setting), pref)

rank_data_list_full$Setting <- factor(rank_data_list_full$Setting,
                                      levels = c(pref, others))

rank_data_list_full

sum(rank_data_list_full$Count[which(rank_data_list_full$method=="SmoothHOOI" & rank_data_list_full$Rank=="3")])

sum(rank_data_list_full$Count[which(rank_data_list_full$method=="FPCA" & rank_data_list_full$Rank=="3")])

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

ggplot(rank_data_list_full, aes(x = Setting, y = Count, fill = Rank)) +
  scale_fill_manual(values = custom_colors, name = "r1 values") +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + 
  facet_nested(rows = NULL, cols = vars(panel, factor(method)), scales = "free_x", strip = custom_strips, space="free_x") +
  labs(x = "Setting", y = "Counts", fill = "r1 values", 
       title = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=30, hjust = 1, vjust = 1, size=12),
        axis.title = element_text(hjust=0.5, size=12),
        axis.text.y = element_text(hjust=0.5, size=12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12),       
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(reverse = FALSE)) 

#ggsave("./rankselect.pdf", dpi=600, width = 12, height=6)

missing_rate_vec_full <- c()
for (i in 1:100){
  missing_rate_vec_full[i] <- full_missstruc[[i]]$missing_rate
}
summary(missing_rate_vec_full)

summary(c(missing_rate_vec_fixr,missing_rate_vec_full))

