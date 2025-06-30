load("./Application/AlgorithmResIdent.Rda") # L_tilde, R_tilde, G_tilde
load("./Application/filtered_ABPM_new.Rda") # filtered_abpm_new

library(ggplot2)
library(dplyr)

g31 <- G_tilde[3,1,]

# Match g31 scores with sleep data 
wake_data <- filtered_abpm_new
ids <- unique(wake_data$ID)
for (i in 1:207){
  id <- ids[i]
  wake_data$g31[wake_data$ID==id] <- g31[i]
}

# Check whether there are multiple patterns within the same hour
id_clean <- wake_data %>%
  group_by(ID, DATE, HOUR) %>%
  summarise(n_val = n_distinct(WAKE), .groups = "drop") %>%
  group_by(ID) %>%
  summarise(conflict = any(n_val > 1), .groups = "drop") %>%
  filter(!conflict) 

# Remove subjects who have more multiple patterns within the same hour
wake_clean <- wake_data[wake_data$ID %in% id_clean$ID, ]

# Check subjects who only have WAKE == 1 and no WAKE == 0
wake_only <- wake_clean %>%
  group_by(ID) %>%
  summarise(has_sleep = any(WAKE == 0), .groups = "drop") %>%
  filter(!has_sleep)
length(wake_only$ID)

# Remove subjects who do not have sleep times
wake_clean <- wake_clean %>% 
  filter(WAKE==0) 

#pdf("./Application/wake_plots.pdf")

# Plot the sleep time for all subjects in wake_reordered
for (i in 1:length(unique(wake_clean$ID))){
  id <- unique(wake_clean$ID)[i]
  print(ggplot(filtered_abpm_new[filtered_abpm_new$ID==id,], aes(x=DATE_TIME, y=WAKE)) + 
          geom_point()+
          labs(title=id))
}
#dev.off()

# Remove subjects who have more than one-day record of sleep time (ID==70217, 70356, 70501, 70567, 70610)
wake_clean <- wake_clean %>%
  filter(!(ID %in% c(70217, 70356, 70501, 70567, 70610)))

# Keep only the sleep time, and reorder the time to put 12am in the middle
wake_reordered <- wake_clean %>% 
  group_by(ID) %>%
  mutate(hour_rescale=(HOUR - 12) %% 24) %>%
  arrange((HOUR - 12) %% 24, .by_group = TRUE)

# Summarize sleep start time and end time
wake_summary <- wake_reordered %>%
  group_by(ID, g31) %>%
  summarise(
    start_hour = min(hour_rescale),
    end_hour = max(hour_rescale),
    mean_hour = mean(c(min(hour_rescale), max(hour_rescale))),
    .groups = "drop"
  ) %>% 
  arrange(g31)

breaks <- seq(min(wake_summary$g31), max(wake_summary$g31), length.out=11)

min(wake_summary$g31)
max(wake_summary$g31)

# Create bins for g31 scores
wake_summary <- wake_summary %>%
  mutate(g31_bin = cut(g31, breaks = breaks, include.lowest = TRUE, labels = FALSE))

# Calculate mean g31 and mean hour for each bin
wake_10_summary <- wake_summary %>%
  group_by(g31_bin) %>%
  summarise(
    mean_g31 = mean(g31),
    mean_10 = mean(mean_hour),
    .groups = "drop"
  ) %>% 
  arrange(mean_10)

ggplot() +
  geom_segment(data=wake_summary, aes(x = start_hour, xend = end_hour, y = g31, yend = g31, color=g31)) +
  geom_segment(data=wake_summary,aes(x = start_hour, xend = start_hour, 
                                     y = g31 - 0.1, yend = g31 + 0.1, color = g31)) +
  geom_segment(data=wake_summary,aes(x = end_hour, xend = end_hour, 
                                     y = g31 - 0.1, yend = g31 + 0.1, color = g31)) +
  geom_point(data = wake_10_summary, aes(x = mean_10, y = mean_g31), color = "red", size = 2) +
  labs(x = "Hour", y = "G score", color = "G score") +
  scale_x_continuous(#limits=c(0,24), 
    breaks = seq(0,24, by=1), labels = c(seq(12,23, by=1), seq(0,12, by=1))) +
  theme_minimal()+
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size=15))

# ggsave("./Application/g31_vs_sleeptime.pdf", dpi=600, width=7, height=5)


