library(ggplot2)
library(ggh4x)
library(dplyr)

Setting <- c( "m20")
Setting_codes <- c("m20")
Setting_labels <- c("20%")

panel_labels <- c(rep("Missingness", 1))

## flexible ranks

load("./ar1_large_full.Rda")

full_results_list <- list(m20 = ar1_large_full)

oracle_M_full <- list()
kcv_M_full <- list()
kcv_hblock_M_full <- list()
fpca_M_full <- list()
cp_M_full <- list()
cp_no_M_full <- list()
mfpca_M_full <- list()

oracle_rank_full <- list()
kcv_rank_full <- list()
kcv_hblock_rank_full <- list()
fpca_rank_full <- list()
cp_rank_full <- list()
cp_no_rank_full <- list()
mfpca_rank_full <- list()

# Reorganize the simulation results into lists
for (set in Setting) {
  oracle_M_full[[set]] <- oracle_rank_full[[set]] <- 
    kcv_M_full[[set]] <- kcv_rank_full[[set]] <- 
    kcv_hblock_M_full[[set]] <- kcv_hblock_rank_full[[set]] <-
    fpca_M_full[[set]] <- fpca_rank_full[[set]]  <- 
    cp_M_full[[set]] <- cp_rank_full[[set]]  <- 
    cp_no_M_full[[set]] <- cp_no_rank_full[[set]]  <-
    mfpca_M_full[[set]] <- mfpca_rank_full[[set]]  <- 
    rep(NA, 100)
  
  
  for (i in 1:100) {
    dat <- full_results_list[[set]][[i]]
    oracle_M_full[[set]][i] <- dat$oracle_loss
    oracle_rank_full[[set]][i] <- dat$oracle_para[1]
    kcv_M_full[[set]][i] <- dat$kcv_loss
    kcv_rank_full[[set]][i] <- dat$kcv_para[1]
    kcv_hblock_M_full[[set]][i] <- dat$kcv_hblock_loss
    kcv_hblock_rank_full[[set]][i] <- dat$kcv_hblock_para[1]
    fpca_M_full[[set]][i] <- dat$fpca_loss
    fpca_rank_full[[set]][i] <- dat$fpca_rank
    cp_M_full[[set]][i] <- dat$cp_loss
    cp_rank_full[[set]][i] <- dat$cp_rank
    cp_no_M_full[[set]][i] <- dat$cp_no_loss
    cp_no_rank_full[[set]][i] <- dat$cp_no_rank
    mfpca_M_full[[set]][i] <- dat$mfpca_loss
    mfpca_rank_full[[set]][i] <- dat$mfpca_rank
    
  }
}

methods_full <- list(Oracle = oracle_M_full, SmoothHOOI = kcv_M_full,hblock=kcv_hblock_M_full, FPCA = fpca_M_full, CP_Ortho = cp_M_full, CP = cp_no_M_full, MFPCA = mfpca_M_full)

data_list_full <- list()

for (method in names(methods_full)) {
  for (i in seq_along(Setting_codes)) {
    data_list_full[[length(data_list_full) + 1]] <- data.frame(
      Method = method,
      panel = panel_labels[i],
      Setting = Setting_labels[i],
      y_value = methods_full[[method]][[Setting_codes[i]]]
    )
  }
}

data_full <- do.call(rbind, data_list_full)

data <- data.frame(data_full)

data$Method <- factor(data$Method, levels=c("Oracle", "SmoothHOOI","hblock", "FPCA",  "MFPCA", "CP_Ortho", "CP"))

panel_limits <- data %>%
  group_by(panel) %>%
  summarise(ymin = min(y_value, na.rm = TRUE),
            ymax = max(y_value, na.rm = TRUE))

data <- data %>%
  left_join(panel_limits, by = "panel")


ggplot(data, aes(x = Method, y = y_value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  labs(y = "Loss of M", fill = "Method", title = "") +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(hjust = 0.5, size = 15),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )
#ggsave("./full_lossM_boxplot_complete.pdf", width=12, height=10)

data_subset <- data.frame(data_full) 
data_subset <- subset(data_subset, Method %in% c("Oracle", "SmoothHOOI","hblock", "FPCA",  "MFPCA", "CP"))

data_subset$Method <- factor(data_subset$Method, levels=c("Oracle", "SmoothHOOI","hblock", "FPCA",  "MFPCA", "CP"))

panel_limits <- data_subset %>%
  group_by(panel) %>%
  summarise(ymin = min(y_value, na.rm = TRUE),
            ymax = max(y_value, na.rm = TRUE))


data_subset <- data_subset %>%
  left_join(panel_limits, by = "panel") 

ggplot(data_subset, aes(x = Method, y = y_value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 0.5) +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  labs(y = "Loss of M", fill = "Method", title = "") +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_text(hjust = 0.5, size = 15),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )
#ggsave("./full_ar1_lossM_noOrtho.pdf", width=8, height=6)

# table() rank for all the methods in c("Oracle", "SmoothHOOI", "FPCA",  "MFPCA", "CP_Ortho", "CP")
rank_tables <- list()
for (method in c("oracle_rank_full", "kcv_rank_full","kcv_hblock_rank_full", "fpca_rank_full", "mfpca_rank_full", "cp_rank_full", "cp_no_rank_full")) {
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
  filter(method %in% c("kcv_rank_full", "kcv_hblock_rank_full", "fpca_rank_full")) 

rank_data_list_full$method <- dplyr::recode(rank_data_list_full$method,
                                            kcv_rank_full = "SmoothHOOI",
                                            kcv_hblock_rank_full = "hblock",
                                            fpca_rank_full = "FPCA")

rank_data_list_full$method <- factor(rank_data_list_full$method, levels=c("SmoothHOOI","hblock", "FPCA"))

rank_data_list_full$Rank <- factor(rank_data_list_full$Rank,
                                   levels = c("2","3","4"))


rank_data_list_full

sum(rank_data_list_full$Count[which(rank_data_list_full$method=="SmoothHOOI" & rank_data_list_full$Rank=="3")])

sum(rank_data_list_full$Count[which(rank_data_list_full$method=="FPCA" & rank_data_list_full$Rank=="3")])

custom_colors <- c("2" = "#FED9A6", "3" = "#B3CDE3", "4" = "#CCEBC5")

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

ggplot(rank_data_list_full, aes(x = method, y = Count, fill = Rank)) +
  scale_fill_manual(values = custom_colors, name = "r1 values") +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + 
  #facet_nested(rows = NULL, cols = vars(panel, factor(method)), scales = "free_x", strip = custom_strips, space="free_x") +
  labs(x = "Method", y = "Counts", fill = "r1 values", 
       title = "") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle=30, hjust = 1, vjust = 1, size=12),
        axis.title = element_text(hjust=0.5, size=12),
        axis.text.y = element_text(hjust=0.5, size=12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12),       
        legend.title = element_text(size = 12)) +
  guides(fill = guide_legend(reverse = FALSE)) 

#ggsave("./rankselect_ar1.pdf", dpi=600, width = 5, height=5)


