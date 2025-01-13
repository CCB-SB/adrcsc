library(data.table)
library(Seurat)
library(ggplot2)
library(stats)
library(dplyr)
library(cowplot)
library(aplot)
source("scripts/helper.R")
source("scripts/Functions_Fig2.R")


pbmc <- readRDS("data/CompleteObjectAnnotated_onlyFirstVisit.rds")
meta.data <- pbmc@meta.data


meta.data$celltype_cluster[meta.data$celltype_cluster == "CD4+ Proliferating T cell"] <-
  "CD4+ Proliferating"
meta.data <-
  meta.data[!meta.data$celltype_cluster %in% c("HSPC / Progenitor TO_DELETE", "Red Blood Cells"), ]


counts_per_sample <- meta.data %>%
  group_by(SCMD) %>%
  mutate(count_pat = n()) %>%
  group_by(SCMD, celltype_cluster, Diagnosis, Sex) %>%
  mutate(per = 100 * n() / count_pat)


counts_per_sample <-
  counts_per_sample[, colnames(counts_per_sample) %in% c("SCMD", "celltype_cluster", "Diagnosis", "Sex", "per")]
counts_per_sample <-
  counts_per_sample[!duplicated(counts_per_sample), ]
counts_per_sample$p <- 1

for (d in unique(counts_per_sample$Diagnosis[counts_per_sample$Diagnosis != "Healthy Control"])) {
  for (ct in unique(counts_per_sample$celltype_cluster)) {
    for (s in unique(counts_per_sample$Sex)) {
      g1 <-
        counts_per_sample$per[counts_per_sample$Diagnosis == d &
                                counts_per_sample$celltype_cluster == ct &
                                counts_per_sample$Sex == s]
      g2 <-
        counts_per_sample$per[counts_per_sample$Diagnosis == "Healthy Control" &
                                counts_per_sample$celltype_cluster == ct &
                                counts_per_sample$Sex == s]
      if (length(g1) < 3 | length(g2) < 3)
        next
      t <- t.test(g1, g2, paired = F)
      counts_per_sample$p[counts_per_sample$Diagnosis == d &
                            counts_per_sample$celltype_cluster == ct &
                            counts_per_sample$Sex == s] <- t$p.value
    }
  }
}

counts_per_sample <-
  counts_per_sample[, colnames(counts_per_sample) %in% c("celltype_cluster", "Diagnosis", "Sex", "p")]
counts_per_sample <-
  counts_per_sample[!duplicated(counts_per_sample), ]

counts_per_sample$p_adj <- 1
counts_per_sample$p_adj[counts_per_sample$p != 1] <-
  p.adjust(counts_per_sample$p[counts_per_sample$p != 1], method = "BH")

c <-
  seq(nrow(meta.data[, colnames(meta.data) %in% c("celltype_cluster", "Diagnosis", "Sex")]))
counts <-
  aggregate(c ~ ., data = meta.data[, colnames(meta.data) %in% c("celltype_cluster", "Diagnosis", "Sex")], FUN =
              length)
counts$CellAbbr <-
  labels_cell[counts$celltype_cluster]#, " cell", "")
counts$CellAbbr <-
  factor(counts$CellAbbr, levels = labels_cell[celltype_cluster_order])
counts$celltype_cluster <-
  factor(counts$celltype_cluster, levels = celltype_cluster_order)
counts <- counts[counts$Diagnosis != "", ]

counts$c <-
  sapply(1:length(counts$celltype_cluster), function(i)
    (counts$c[i] / sum(counts$c[counts$Sex == counts$Sex[i] &
                                  counts$Diagnosis == counts$Diagnosis[i]])))
counts$c <-
  sapply(1:length(counts$celltype_cluster), function(i)
    ((counts$c[i]) / (counts$c[counts$celltype_cluster == counts$celltype_cluster[i] &
                                 counts$Sex == counts$Sex[i] &
                                 counts$Diagnosis == "Healthy Control"]))[1])
counts <- counts[counts$Diagnosis != "Healthy Control", ]

counts <-
  merge(counts,
        counts_per_sample,
        by = c("Sex", "Diagnosis", "celltype_cluster"))


counts$literature <- F
counts$literature[counts$Diagnosis == "Parkinson's Disease only" &
                    counts$celltype_cluster %in% unique(counts$celltype_cluster)[grepl("CD8", unique(counts$celltype_cluster)) |
                                                                                   grepl("CD4", unique(counts$celltype_cluster)) |
                                                                                   grepl("Th1", unique(counts$celltype_cluster)) |
                                                                                   grepl("B", unique(counts$celltype_cluster))]] <- T
counts$literature[counts$Diagnosis == "Alzheimer's disease" &
                    counts$celltype_cluster %in% unique(counts$celltype_cluster)[grepl("Monocyte", unique(counts$celltype_cluster)) |
                                                                                   grepl("B", unique(counts$celltype_cluster)) |
                                                                                   grepl("CD8", unique(counts$celltype_cluster))]] <- T
counts$literature[counts$Diagnosis == "Alzheimer's disease" &
                    counts$celltype_cluster %in% c("Treg CD4+ cell", "CD4+ Tfh cell", "Th1/Th2/Th17 cell")] <-
  T

plot <-
  ggplot(counts[counts$Sex == "male", ], aes(x = celltype_cluster, y = name2short[Diagnosis])) +
  geom_tile(
    aes(
      fill = c,
      color = ifelse(p_adj < 0.05, "black", "white")
    ),
    width = 0.9,
    height = 0.9,
    linewidth = .2
  ) +
  geom_point(data = counts[counts$Sex == "male" &
                             counts$literature, ],
             aes(x = celltype_cluster, y = name2short[Diagnosis]),
             size = .5) +
  scale_color_manual(values = c("black", "white")) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 1,
    limits = c(0, 3.8)
  ) +
  theme_adrc() + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) + facet_grid(Sex ~ .) +
  scale_x_discrete(name = "", labels = labels_cell[meta.data$celltype_cluster]) +
  scale_y_discrete(name = "") + labs(fill = 'Rel. Change\nof the\n Cell-type\nproportion') + guides(color = "none")

plot2 <-
  ggplot(counts[counts$Sex == "female", ], aes(x = celltype_cluster, y = name2short[Diagnosis])) +
  geom_tile(
    aes(
      fill = c,
      color = ifelse(p_adj < 0.05, "black", "white")
    ),
    width = 0.9,
    height = 0.9,
    linewidth = .2
  ) +  guides(color = "none") +
  geom_point(data = counts[counts$Sex == "female" &
                             counts$literature, ],
             aes(x = celltype_cluster, y = name2short[Diagnosis]),
             size = .5) +
  scale_color_manual(values = c("black", "white")) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 1,
    limits = c(0, 3.8)
  ) +
  theme_adrc() + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) + facet_grid(Sex ~ .) +
  scale_x_discrete(name = "", labels = labels_cell[meta.data$celltype_cluster]) +
  scale_y_discrete(name = "") + labs(fill = 'Rel. Change\nof the\n Cell-type\nproportion') +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

labels = ggplot(counts, aes(x = celltype_cluster, y = 1, fill = celltype_cluster)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill = "none")

p <-
  (plot2 %>% insert_top(labels, height = .25)) %>% insert_bottom(plot, height = 1.0)

ggsave("figures/figure_2_c.svg",
       p,
       width = 5,
       height = 3)

write.csv(counts, "SourceData/Figure 2c.csv")
