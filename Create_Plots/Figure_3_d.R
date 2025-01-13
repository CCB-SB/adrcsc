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


## Figure 3 d)

degs_male <-
  read.csv(
    "results_male/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )
degs_male <-
  degs_male[degs_male$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"), ]
degs_female <-
  read.csv(
    "results_female/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )
degs_male$Sex = "male"
degs_female$Sex = "female"
degs_female <-
  degs_female[degs_female$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"), ]
dge <-
  merge(degs_male, degs_female, by = c("gene", "cluster_id", "contrast"))

dge$color <- "none"
dge$color[dge$p_adj.loc.x < 0.05] <- "male"
dge$color[dge$p_adj.loc.y < 0.05] <- "female"
dge$color[dge$p_adj.loc.y < 0.05  &
            dge$p_adj.loc.x < 0.05] <- "male and female"

m <-
  dge[dge$cluster_id %in% c("CD4+ Tcm/Tscm T cell", "NKT-like cell"), ]

m$cl_abbr <- m$cluster_id

p <- ggplot(m, aes(x = logFC.x, y = logFC.y)) +
  geom_point(data = m[m$color == "none", ],
             aes(x = logFC.x, y = logFC.y, color = color),
             alpha = 0.7) +
  geom_point(data = m[m$color != "none", ],
             aes(x = logFC.x, y = logFC.y, color = color),
             alpha = 0.7) +
  facet_grid(cols = vars(contrast), rows = vars(cluster_id)) +  xlim(-3.5, 3.5) + ylim(-3.5, 3.5) +
  scale_fill_gradient2(low = "#d6604d", mid = "white", high = "#4393c3") +
  scale_colour_manual(
    name = "significant in",
    values = c(
      "male and female" = "black",
      "female" = "#9C489A",
      "none" = "lightgray",
      "male" = "#5D7322",
      "gray" = "gray"
    )
  ) +
  stat_cor(
    method = "pearson",
    label.x = -2,
    label.y = 3,
    size = 2.5
  ) +
  geom_smooth(method = lm, color = "black") +
  geom_hline(yintercept = 0,
             color = "lightgrey",
             linetype = "dashed") +
  geom_vline(xintercept = 0,
             color = "lightgrey",
             linetype = "dashed") +
  geom_hline(yintercept = 0.5,
             color = "lightgrey",
             linetype = "dotted") +
  geom_hline(yintercept = -0.5,
             color = "lightgrey",
             linetype = "dotted") +
  geom_vline(xintercept = 0.5,
             color = "lightgrey",
             linetype = "dotted") +
  geom_vline(xintercept = -0.5,
             color = "lightgrey",
             linetype = "dotted") + xlab("log2 Fold-change in male") + ylab("log2 Fold-change in female") +
  theme_adrc() + theme(aspect.ratio = 1, legend.position = "bottom")

figure_d <-  fill_title(p, color_v)
ggsave(
  "figures/figure_3_d.svg",
  as.ggplot(figure_d),
  height = 7,
  width = 9,
  unit = "cm"
)

write.csv(m, "SourceData/figure 3d.csv")
