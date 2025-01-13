library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(aplot)
library(ComplexHeatmap)
library(viridis)
library(ggplotify)
library(ggh4x)

source("scripts/helper.R")



degs_male <-
  read.csv(
    "results_male/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )
degs_male <-
  degs_male[degs_male$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"), ]
degs_male$Sex = "male"
degs_female <-
  read.csv(
    "results_female/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )
degs_female$Sex = "female"
degs_female <-
  degs_female[degs_female$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"), ]

degs_male_t <- degs_male[(degs_male$gene %in% degs_female$gene) , ]
degs_female_t <-
  degs_female[(degs_female$gene %in% degs_male$gene), ]

degs <-
  degs_male_t %>% right_join(degs_female_t, by = c("gene", "cluster_id", "contrast"))

Correlation <- NA
for (cellt in unique(degs$cluster_id)) {
  for (cont in unique(degs$contrast[degs$cluster_id == cellt])) {
    tryCatch({
      corr <-
        cor.test(degs$logFC.x[degs$contrast == cont &
                                degs$cluster_id == cellt],
                 degs$logFC.y[degs$contrast == cont &
                                degs$cluster_id == cellt],
                 method = c("pearson", "kendall", "spearman"))
      res <-
        data.frame(
          Value = corr$estimate,
          p_val = corr$p.value,
          Celltype = cellt,
          Comparison = cont
        )
      if (any(is.na(Correlation))) {
        Correlation <- res
      } else {
        Correlation <- rbind(Correlation, res)
      }
    },
    error = function(cond) {
      print(cond)
    })
  }
}
Correlation$p_val_adj <- p.adjust(Correlation$p_val, method = "BH")

# Figure 3 b)


degs_male <-
  degs_male[degs_male$p_adj.loc < 0.05 & abs(degs_male$logFC) > 0.5, ]
degs_female <-
  degs_female[degs_female$p_adj.loc < 0.05 &
                abs(degs_female$logFC) > 0.5, ]
degs_dereg <-
  full_join(degs_male[, colnames(degs_male) %in% c("gene", "cluster_id", "contrast")],
            degs_female[, colnames(degs_female) %in% c("gene", "cluster_id", "contrast")])
c <-
  seq(nrow(degs_dereg[, colnames(degs_dereg) %in% c("cluster_id", "contrast")]))
counts <-
  aggregate(c ~ ., data = degs_dereg[, colnames(degs_dereg) %in% c("cluster_id", "contrast")], FUN =
              length)
counts <-
  counts[!counts$cluster_id %in% c("HSPC / Progenitor TO_DELETE", "Red Blood Cells"), ]
counts$Celltype <- counts$cluster_id
counts$Comparison <- counts$contrast


Correlation$Comparison <-
  factor(Correlation$Comparison,
         levels = c("PDvsHC", "MCIvsHC", "ADvsHC"))
Correlation$Celltype[Correlation$Celltype == "CD4+ Proliferating T cell"] <-
  "CD4+ Proliferating"
Correlation <-
  Correlation[!Correlation$Celltype %in% c("HSPC / Progenitor TO_DELETE", "Red Blood Cells"), ]
Correlation$CellAbbr <- labels_cell[Correlation$Celltype]
Correlation$CellAbbr <-
  factor(Correlation$CellAbbr, levels = labels_cell[celltype_cluster_order])
Correlation$Celltype <-
  factor(Correlation$Celltype, levels = celltype_cluster_order)


data <- merge(Correlation, counts, by = c("Celltype", "Comparison"))

data$p_val[data$p_val == 0] <- min(data$p_val[data$p_val > 0])

plot <-
  ggplot(data, aes(
    x = c,
    y = Value,
    size = -log10(p_val),
    fill = Comparison
  )) +
  geom_point(color = "black", pch = 21) +
  theme_adrc() +
  scale_size(range = c(0, 2)) +
  facet_grid(. ~ Comparison) +
  scale_fill_manual(values = color_v, guide = "none") +
  geom_text_repel(
    data = subset(
      data,
      p_val < 0.05 &
        labels_cell[as.character(Celltype)] %in% c("Naive CD4+ T", "Prolif. CD8+ T", "CD56- NK")
    ),
    aes(label = labels_cell[as.character(Celltype)]),
    color = "black",
    size = 2
  ) +
  labs(size = "-log10\nadj.\np-value", Comparison = "") +
  geom_hline(yintercept = 0, linetype = 'dashed')  +
  xlab("Number of DEGs") + ylab("Correlation") +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    aspect.ratio = 1.3,
    plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
    panel.grid = element_line(colour = "gray92"),
    panel.border = element_rect(color = "black")
  )

figure_b <- as.ggplot(fill_title(plot, color_v))

ggsave(
  "figures/figure_3_b.svg",
  figure_b,
  height = 4,
  width = 9,
  unit = "cm"
)

write.csv(data, "SourceData/figure_3_bc.csv")

## Figure 3 c)



labels = ggplot(Correlation, aes(x = Celltype, y = 1, fill = Celltype)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill = "none")

plot <-
  ggplot(Correlation, aes(x = Celltype, y = Comparison, fill = Value)) +
  geom_tile() +
  geom_point(aes(
    size = ifelse(
      Correlation$p_val_adj < 0.05,
      ifelse(Correlation$Value > 0, "pos", "neg"),
      "no_dot"
    ),
    shape = ifelse(
      Correlation$p_val_adj < 0.05,
      ifelse(Correlation$Value > 0, "pos", "neg"),
      "no_dot"
    )
  )) +
  scale_size_manual(values = c(pos = 0.5, neg = 0.5, no_dot = NA),
                    guide = "none") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    limits = c(-1, 1)
  ) +
  scale_shape_manual(values = c(pos = 19, neg = 19, no_dot = NA),
                     guide = "none") +
  theme_adrc() + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.7,
    hjust = 1
  )) +
  scale_x_discrete(name = "", labels = labels_cell[as.character(Correlation$Celltype)]) +
  scale_y_discrete(name = "") +
  labs(fill = 'Correlation') +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  theme(legend.key.size = unit(3, 'mm'))

figure_c <- plot <- plot %>% insert_top(labels, height = .3)
ggsave(
  "figures/figure_3_c.svg",
  plot,
  width = 9,
  height = 4,
  unit = "cm"
)
