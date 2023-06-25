scores_df <- read.csv("all_200_rescaled_vf.csv")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
scores_df <- scores_df %>% mutate(rank=dense_rank(desc(-AUC_mean)))

scores_df$median_group <- ifelse(scores_df$AUC_mean > median(scores_df$AUC_mean), "High", "Low")

PRISM <- readRDS("./targeted_survival/files/PSet_PRISM.rds")
PRISM <- updateObject(PRISM)
treatment <- PRISM@treatment[,colnames(PRISM@treatment) %in% c("moa","treatmentid", "cid", "smiles", "phase","indication")]

find_closest_match <- function(drug, treatment_table, threshold = 3) {
  treatmentid <- treatment_table$treatmentid
  distances <- adist(drug, treatmentid)
  
  closest_index <- which.min(distances)
  min_distance <- min(distances)
  
  if (min_distance <= threshold) {
    return(treatment_table[closest_index, ])
  } else {
    return(data.frame(moa = NA, treatmentid = NA, cid = NA, smiles = NA, phase = NA, indication = NA))
  }
}

result <- scores_df %>%
  rowwise() %>%
  mutate(closest_treatment = list(find_closest_match(drugs, treatment))) %>%
  unnest(cols = closest_treatment) %>%
  select(-treatmentid) 

##edge outgoing
png("dyn_plots/edge_outgoing_bar.png", units="in", width=4, height=6, res=300, type = "cairo")
ggplot(scores_df, aes(x = median_group, y = edge_outgoing, fill = median_group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    hide.ns = TRUE,
    bracket.size = 0.5,
    label.y = max(scores_df$edge_outgoing) * 0.95,
    label.x = 1.5
  ) +
  labs(
    x = "AAC Score",
    y = "Edge Outgoing Vector Field Score"
  ) +
  theme_bw() +
  theme(
    text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Low" = "#1b9e77", "High" = "#7570b3")) +
  scale_x_discrete(labels = c("High (n = 35)","Low (n = 35)"))
dev.off()


##core incoming
png("dyn_plots/core_incoming_bar.png", units="in", width=4, height=6, res=300, type = "cairo")
ggplot(scores_df, aes(x = median_group, y = core_incoming, fill = median_group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    hide.ns = TRUE,
    bracket.size = 0.5,
    label.y = max(scores_df$edge_outgoing) * 0.95,
    label.x = 1.5
  ) +
  labs(
    x = "AAC Score",
    y = "Core Incoming Vector Field Score"
  ) +
  theme_bw() +
  theme(
    text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Low" = "#1b9e77", "High" = "#7570b3")) +
  scale_x_discrete(labels = c("High (n = 35)","Low (n = 35)"))
dev.off()

###drug class analysis

#filter for only classes with multiple drugs
filtered_result <- result %>%
  group_by(moa) %>%
  filter(n() > 1) %>%
  ungroup()

filtered_result <- filtered_result[!is.na(filtered_result$moa),]

png("dyn_plots/edge_outgoing_drug_class.png", units="in", width=8, height=6, res=300, type = "cairo")
ggplot(filtered_result, aes(x = moa, y = edge_outgoing, fill = moa)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    hide.ns = TRUE,
    bracket.size = 0.5,
    label.y = max(scores_df$edge_outgoing) * 0.95,
    label.x = 1.5
  ) +
  labs(
    x = "Drug MOA",
    y = "Edge Outgoing Vector Field Score"
  ) +
  theme_bw() +
  theme(
    text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
dev.off()

png("dyn_plots/core_incoming_drug_class.png", units="in", width=8, height=6, res=300, type = "cairo")

ggplot(filtered_result, aes(x = moa, y = core_incoming, fill = moa)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    hide.ns = TRUE,
    bracket.size = 0.5,
    label.y = max(scores_df$core_incoming) * 0.95,
    label.x = 1.5
  ) +
  labs(
    x = "Drug MOA",
    y = "Core Incoming Vector Field Score"
  ) +
  theme_bw() +
  theme(
    text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
dev.off()

write.csv(result, file = "Dyanmo_Vf_Pertrub.csv")
