library(tidyr)
library(ggplot2)
library(dplyr)
library(pheatmap)

# Plasmid cluster gene presence absence 

cluster1 <- read.csv(file='results/plasmid/cluster_gene_pw/complete_contigs/presence_absence_matrix_cluster_1.csv')
cluster2 <- read.csv(file='results/plasmid/cluster_gene_pw/complete_contigs/presence_absence_matrix_cluster_2.csv')


# Reshape the data into long format
df_long <- cluster1 %>%
  pivot_longer(cols = -c(Key, Description), names_to = "Sample", values_to = "Presence")

# Convert "+" to 1 and others to 0 for visualization
df_long <- df_long %>%
  mutate(Presence = ifelse(Presence == "+", 1, 0))

# Calculate the total presence for each Description and reorder
description_order <- df_long %>%
  group_by(Description) %>%
  summarize(total_presence = sum(Presence, na.rm = TRUE)) %>%
  arrange(desc(total_presence)) %>%
  pull(Description)

df_long <- df_long %>%
  mutate(Description = factor(Description, levels = description_order))

# Create the heatmap
heatmap_plot_1 <- ggplot(df_long, aes(y = Description, x = Sample, fill = factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "blue"), name = "Presence") +
  labs(title = "Plasmid Cluster 1",
       x = "Genome",
       y = "OLG Function") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.y = element_text(hjust = 1, size = 6))

# Print the heatmap
heatmap_plot_1
# Save
ggsave("results/plasmid/cluster_gene_pw/cluster_1_complete_pam.png", plot = heatmap_plot_1, width = 12, height = 40, dpi = 300)




#
# Plasmid cluster 2
# 

# Reshape the data into long format
df_long <- cluster2 %>%
  pivot_longer(cols = -c(Key, Description), names_to = "Sample", values_to = "Presence")

# Convert "+" to 1 and others to 0 for visualization
df_long <- df_long %>%
  mutate(Presence = ifelse(Presence == "+", 1, 0))

# Calculate the total presence for each Description and reorder
description_order <- df_long %>%
  group_by(Description) %>%
  summarize(total_presence = sum(Presence, na.rm = TRUE)) %>%
  arrange(desc(total_presence)) %>%
  pull(Description)

df_long <- df_long %>%
  mutate(Description = factor(Description, levels = description_order))

# Create the heatmap
heatmap_plot_2 <- ggplot(df_long, aes(y = Description, x = Sample, fill = factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "blue"), name = "Presence") +
  labs(title = "Plasmid Cluster 2",
       x = "Genome",
       y = "OLG Function") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1, size = 6))

# Print the heatmap
heatmap_plot_2

ggsave("results/plasmid/cluster_gene_pw/cluster_2_complete_pam.png", plot = heatmap_plot_2, width = 12, height = 20, dpi = 300)

