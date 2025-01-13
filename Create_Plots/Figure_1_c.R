library(Seurat)
library(ggplot2)
library(data.table)

source("scripts/helper.R")

adrc <- readRDS("data/CompleteObjectAnnotated_final.rds")

############## Figure 1d)
adrc = SetIdent(adrc, value = "L2")
figure1_d = DimPlot(adrc, reduction = "umap", label = F) +
  scale_color_manual(values = color_v) +
  NoLegend() +
  scale_x_continuous(expand = expand_scale(mult = c(0.025, 0.025))) +
  theme_void() +
  theme(aspect.ratio = 1) +
  theme(legend.position = "none")

write.csv(figure1_d$data, "SourceData/Figure1c.csv")

ggsave("figures/figure_1_c.pdf",
       figure1_d,
       width = 5,
       height = 5)

figure1_d = DimPlot(adrc, reduction = "umap", label = T) +
  scale_color_manual(values = color_v) +
  NoLegend() +
  scale_x_continuous(expand = expand_scale(mult = c(0.025, 0.025))) +
  theme_void() +
  theme(aspect.ratio = 1) +
  theme(legend.position = "none")

ggsave("figures/figure_1_c_label.pdf",
       figure1_d,
       width = 5,
       height = 5)
