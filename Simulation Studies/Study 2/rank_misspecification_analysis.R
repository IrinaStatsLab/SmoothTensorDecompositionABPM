library(ggplot2)
library(ggh4x)
library(dplyr)

Setting <- c("m0", "m10", "m20", "m50", "struc", "p50", "p500", "n0", "n05", "n15", "n2")
Setting_codes <- c("m0", "m10", "m20", "m50", "struc", "p50", "m20", "p500", "n0", "n05", "m20", "n15", "n2")
Setting_labels <- c("0%", "10%", "20%", "50%", "structured", 
                    "n=50", "n=200", "n=500", 
                    "Noise: 0", "Noise: 0.5",
                    "Noise: 1", "Noise: 1.5", "Noise: 2")

panel_labels <- c(rep("Missingness", 5), rep("Sample size", 3), rep("Noise level", 5))


rda_files <- list.files("~/Documents/GitHub/TensorABPMSmooth/Simulation_Studies_R1_Final/rank_misspecification_ev", pattern = "\\.Rda$", full.names = TRUE)
for (f in rda_files) {
  load(f)
}

rank_results_list <- list(m0 = miss0, m10 = miss10,
                          m20 = miss20, m50 = miss50,
                          struc = missstruc,
                          p50 = p50, p500 = p500,
                          n0 = noise0, 
                          n05 = noise05, n15 = noise15,
                          n2 = noise2)

rank_data_list <- list()

for (set in Setting) {
  
  res_set <- rank_results_list[[set]]
  
  for (i in seq_along(res_set)) {
    
    dat <- res_set[[i]]
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(3,2)",
      loss_M_full = dat$kcv_true_loss$loss_M,
      loss_M_truncated = dat$kcv_true_loss$loss_M,
      loss_L = dat$kcv_true_loss$loss_L,
      loss_R = dat$kcv_true_loss$loss_R,
      lambda = dat$lambda_true
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(4,2)",
      loss_M_full = dat$kcv_42_lossM$loss_M,
      loss_M_truncated = dat$kcv_42_loss_truncated$loss_M,
      loss_L = dat$kcv_42_loss_truncated$loss_L,
      loss_R = dat$kcv_42_loss_truncated$loss_R,
      lambda = dat$lambda_42
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(6,2)",
      loss_M_full = dat$kcv_62_lossM$loss_M,
      loss_M_truncated = dat$kcv_62_loss_truncated$loss_M,
      loss_L = dat$kcv_62_loss_truncated$loss_L,
      loss_R = dat$kcv_62_loss_truncated$loss_R,
      lambda = dat$lambda_62
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(3,3)",
      loss_M_full = dat$kcv_33_lossM$loss_M,
      loss_M_truncated = dat$kcv_33_loss_truncated$loss_M,
      loss_L = dat$kcv_33_loss_truncated$loss_L,
      loss_R = dat$kcv_33_loss_truncated$loss_R,
      lambda = dat$lambda_33
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(4,3)",
      loss_M_full = dat$kcv_43_lossM$loss_M,
      loss_M_truncated = dat$kcv_43_loss_truncated$loss_M,
      loss_L = dat$kcv_43_loss_truncated$loss_L,
      loss_R = dat$kcv_43_loss_truncated$loss_R,
      lambda = dat$lambda_43
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(6,3)",
      loss_M_full = dat$kcv_63_lossM$loss_M,
      loss_M_truncated = dat$kcv_63_loss_truncated$loss_M,
      loss_L = dat$kcv_63_loss_truncated$loss_L,
      loss_R = dat$kcv_63_loss_truncated$loss_R,
      lambda = dat$lambda_63
    )
  }
}

rank_data_raw <- bind_rows(rank_data_list)

############################################################
## Add panel information
############################################################

setting_map <- data.frame(
  raw_setting = Setting_codes,
  panel = panel_labels,
  Setting = Setting_labels
)

# Important:
# m20 appears three times in setting_map, so the same baseline
# simulation results are deliberately reused in three panels.
rank_data <- rank_data_raw %>%
  left_join(
    setting_map,
    by = "raw_setting",
    relationship = "many-to-many"
  )


############################################################
## Factor ordering
############################################################

rank_data$rank <- factor(
  rank_data$rank,
  levels = c(
    "(3,2)",
    "(4,2)",
    "(6,2)",
    "(3,3)",
    "(4,3)",
    "(6,3)"
  )
)

rank_data$Setting <- factor(
  rank_data$Setting,
  levels = Setting_labels
)

rank_data$panel <- factor(
  rank_data$panel,
  levels = c(
    "Missingness",
    "Sample size",
    "Noise level"
  )
)

rank_colors <- c(
  "(3,2)" = "#B79F00",  
  "(4,2)" = "#E69F00",
  "(6,2)" = "#D55E00",
  
  "(3,3)" = "#73C8D3",
  "(4,3)" = "#3690B5",
  "(6,3)" = "#25618B"
)

## Loss of L for lambda selected for each fitted rank

plot_L <- ggplot(
  rank_data,
  aes(
    x = Setting,
    y = loss_L,
    fill = rank
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 1,
    outlier.size = 0.5
  ) +
  scale_fill_manual(values = rank_colors) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_nested(
    ~ panel,
    scales = "free_x",
    space = "free_x",
    independent = "none"
  ) +
  guides(
    fill = guide_legend(
      byrow = TRUE,
      nrow = 1
    )
  ) +
  labs(
    x = "",
    # y = "Loss of L after truncation (optimal)",
    y = expression("Loss of L after truncation: " * lambda * " selected for each fitted rank"),
    fill = "Fitted rank"
  ) +
  #theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 13
    ),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

plot_L
#ggsave("./rank_miss_lossL_trunc_ylim.pdf", width=12, height=6)

## Selected lambda under each fitted rank

ranks <- c("(3,2)", "(4,2)", "(6,2)", "(3,3)", "(4,3)", "(6,3)")
lambda_cols <- c(
  "lambda_true", "lambda_42", "lambda_62",
  "lambda_33", "lambda_43", "lambda_63"
)


lambda_data_raw <- imap_dfr(rank_results_list[Setting], \(res_set, set) {
  imap_dfr(res_set, \(dat, i) {
    tibble(
      raw_setting = set,
      replicate = i,
      rank = ranks,
      lambda = unlist(dat[lambda_cols])
    )
  })
})

lambda_data <- lambda_data_raw %>%
  inner_join(setting_map, by = "raw_setting") %>%
  mutate(
    rank = factor(rank, levels = ranks),
    panel = factor(
      panel,
      levels = c("Missingness", "Sample size", "Noise level")
    ),
    Setting = factor(
      Setting,
      levels = setting_map$Setting
    )
  )


p_lambda <- ggplot(
  lambda_data,
  aes(Setting, lambda, fill = rank)
) +
  geom_boxplot(
    position = position_dodge(0.8),
    width = 1,
    outlier.size = 0.7
  ) +
  scale_fill_manual(values = rank_colors) +
  facet_nested(
    ~ panel,
    scales = "free_x",
    space = "free_x",
    independent = "none"
  ) +
  labs(
    x = NULL,
    y = expression("CV-selected " * lambda),
    fill = "Fitted rank"
  ) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 1)) +
  #theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    strip.text = element_text(size = 13),
    strip.background = element_rect(fill = "grey85"),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
    #panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)
  )

p_lambda

#ggsave("./rank_miss_selected_lambda.pdf", width=12, height=6)


sepL_m20_list <- list()

EV_names <- c(
  "(3,2)" = "EV_true",
  "(4,2)" = "EV_42",
  "(6,2)" = "EV_62",
  "(3,3)" = "EV_33",
  "(4,3)" = "EV_43",
  "(6,3)" = "EV_63"
)

for (i in seq_along(miss20)) {
  
  dat <- miss20[[i]]
  
  for (r in names(EV_names)) {
    
    ev_obj <- dat[[EV_names[r]]]
    x <- ev_obj$sep_EV_L
    
    sepL_m20_list[[length(sepL_m20_list) + 1]] <- data.frame(
      replicate = i,
      rank = r,
      component = paste0("L", seq_along(x)),
      sep_EV_L = as.numeric(x)
    )
  }
}

sepL_m20 <- bind_rows(sepL_m20_list)

sepL_m20$rank <- factor(
  sepL_m20$rank,
  levels = c(
    "(3,2)", "(4,2)", "(6,2)",
    "(3,3)", "(4,3)", "(6,3)"
  )
)

sepL_m20$component <- factor(
  sepL_m20$component,
  levels = paste0("L", 1:6)
)

ggplot(sepL_m20, aes(x = component, y = 100 * sep_EV_L, fill = rank)) +
  scale_fill_manual(values = rank_colors) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 1,
    outlier.size = 0.6
  ) +
  labs(
    x = "",
    y = "Incremental Explained Variability in L (%)",
    fill = "Fitted rank"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  #theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

#ggsave("./rank_miss_L_EV.pdf", width=10, height=5)

sepR_m20_list <- list()

EV_names <- c(
  "(3,2)" = "EV_true",
  "(4,2)" = "EV_42",
  "(6,2)" = "EV_62",
  "(3,3)" = "EV_33",
  "(4,3)" = "EV_43",
  "(6,3)" = "EV_63"
)

for (i in seq_along(miss20)) {
  
  dat <- miss20[[i]]
  
  for (r in names(EV_names)) {
    
    ev_obj <- dat[[EV_names[r]]]
    x <- ev_obj$sep_EV_R
    
    sepR_m20_list[[length(sepR_m20_list) + 1]] <- data.frame(
      replicate = i,
      rank = r,
      component = paste0("R", seq_along(x)),
      sep_EV_R = as.numeric(x)
    )
  }
}

sepR_m20 <- bind_rows(sepR_m20_list)

sepR_m20$rank <- factor(
  sepR_m20$rank,
  levels = c(
    "(3,2)", "(4,2)", "(6,2)",
    "(3,3)", "(4,3)", "(6,3)"
  )
)

sepR_m20$component <- factor(
  sepR_m20$component,
  levels = c("R1", "R2", "R3")
)

ggplot(
  sepR_m20,
  aes(
    x = component,
    y = 100 * sep_EV_R,
    fill = rank
  )
) +
  scale_fill_manual(values = rank_colors) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 1,
    outlier.size = 0.6
  ) +
  labs(
    x = "",
    y = "Incremental Explained Variability in R (%)",
    fill = "Fitted rank"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  coord_cartesian(ylim = c(0, 50)) +
  #theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

#ggsave("./rank_miss_R_EV.pdf", width=10, height=5)

# Loss of L with lambda selected under true rank

rda_files <- list.files("./rank_misspecification_samelambda", pattern = "\\.Rda$", full.names = TRUE)
for (f in rda_files) {
  load(f)
}

rank_results_list <- list(m0 = miss0, m10 = miss10,
                          m20 = miss20, m50 = miss50,
                          struc = missstruc,
                          p50 = p50, p500 = p500,
                          n0 = noise0, 
                          n05 = noise05, n15 = noise15,
                          n2 = noise2)

rank_data_list <- list()

for (set in Setting) {
  
  res_set <- rank_results_list[[set]]
  
  for (i in seq_along(res_set)) {
    
    dat <- res_set[[i]]
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(3,2)",
      loss_M_full = dat$kcv_true_loss$loss_M,
      loss_M_truncated = dat$kcv_true_loss$loss_M,
      loss_L = dat$kcv_true_loss$loss_L,
      loss_R = dat$kcv_true_loss$loss_R,
      lambda = dat$lambda_true
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(4,2)",
      loss_M_full = dat$kcv_42_lossM$loss_M,
      loss_M_truncated = dat$kcv_42_loss_truncated$loss_M,
      loss_L = dat$kcv_42_loss_truncated$loss_L,
      loss_R = dat$kcv_42_loss_truncated$loss_R
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(6,2)",
      loss_M_full = dat$kcv_62_lossM$loss_M,
      loss_M_truncated = dat$kcv_62_loss_truncated$loss_M,
      loss_L = dat$kcv_62_loss_truncated$loss_L,
      loss_R = dat$kcv_62_loss_truncated$loss_R
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(3,3)",
      loss_M_full = dat$kcv_33_lossM$loss_M,
      loss_M_truncated = dat$kcv_33_loss_truncated$loss_M,
      loss_L = dat$kcv_33_loss_truncated$loss_L,
      loss_R = dat$kcv_33_loss_truncated$loss_R
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(4,3)",
      loss_M_full = dat$kcv_43_lossM$loss_M,
      loss_M_truncated = dat$kcv_43_loss_truncated$loss_M,
      loss_L = dat$kcv_43_loss_truncated$loss_L,
      loss_R = dat$kcv_43_loss_truncated$loss_R
    )
    
    rank_data_list[[length(rank_data_list) + 1]] <- data.frame(
      raw_setting = set,
      replicate = i,
      rank = "(6,3)",
      loss_M_full = dat$kcv_63_lossM$loss_M,
      loss_M_truncated = dat$kcv_63_loss_truncated$loss_M,
      loss_L = dat$kcv_63_loss_truncated$loss_L,
      loss_R = dat$kcv_63_loss_truncated$loss_R
    )
  }
}

rank_data_raw <- bind_rows(rank_data_list)

rank_data <- rank_data_raw %>%
  left_join(
    setting_map,
    by = "raw_setting",
    relationship = "many-to-many"
  )


rank_data$rank <- factor(
  rank_data$rank,
  levels = c(
    "(3,2)",
    "(4,2)",
    "(6,2)",
    "(3,3)",
    "(4,3)",
    "(6,3)"
  )
)

rank_data$Setting <- factor(
  rank_data$Setting,
  levels = Setting_labels
)

rank_data$panel <- factor(
  rank_data$panel,
  levels = c(
    "Missingness",
    "Sample size",
    "Noise level"
  )
)


plot_L <- ggplot(
  rank_data,
  aes(
    x = Setting,
    y = loss_L,
    fill = rank
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 1,
    outlier.size = 0.5
  ) +
  scale_fill_manual(values = rank_colors) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_nested(
    ~ panel,
    scales = "free_x",
    space = "free_x",
    independent = "none"
  ) +
  guides(
    fill = guide_legend(
      byrow = TRUE,
      nrow = 1
    )
  ) +
  labs(
    x = "",
    # y = "Loss of L after truncation (optimal)",
    y = expression("Loss of L after truncation: " * lambda * " selected under true rank"),
    fill = "Fitted rank"
  ) +
  #theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 13
    ),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

plot_L
#ggsave("./rank_miss_samelambda_lossL_trunc_ylim.pdf", width=12, height=6)









