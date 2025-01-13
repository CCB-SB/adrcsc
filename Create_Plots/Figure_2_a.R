library(ggplot2)
library(cowplot)
library(data.table)
library(Seurat)
library(grid)


source("scripts/helper.R")
source("scripts/Functions_Fig2.R")

pbmc <- readRDS("data/CompleteObjectAnnotated_onlyFirstVisit.rds")

pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]
pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels = diagnosis_order)

umap_data = data.table(
  pbmc@reductions$umap@cell.embeddings,
  biogroup = pbmc$biogroup_short,
  celltype = pbmc$celltype,
  celltype_cluster = pbmc$celltype_cluster,
  sample = pbmc$Sample,
  Sex = pbmc$Sex
)
setnames(umap_data, c("UMAP_1", "UMAP_2"), c("V1", "V2"))

p_mci_f = plot_density_difference(umap_data[biogroup == "MCI" &
                                              Sex == "female"], umap_data[biogroup == "HC" &
                                                                            Sex == "female"]) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_ad_f = plot_density_difference(umap_data[biogroup == "AD" &
                                             Sex == "female"], umap_data[biogroup == "HC" &
                                                                           Sex == "female"]) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_pd_f = plot_density_difference(umap_data[biogroup == "PD" &
                                             Sex == "female"], umap_data[biogroup == "HC" &
                                                                           Sex == "female"]) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_mci_m = plot_density_difference(umap_data[biogroup == "MCI" &
                                              Sex == "male"], umap_data[biogroup == "HC" &
                                                                          Sex == "male"]) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_ad_m = plot_density_difference(umap_data[biogroup == "AD" &
                                             Sex == "male"], umap_data[biogroup == "HC" &
                                                                         Sex == "male"]) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_pd_m = plot_density_difference(umap_data[biogroup == "PD" &
                                             Sex == "male"], umap_data[biogroup == "HC" &
                                                                         Sex == "male"]) +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

legend <- cowplot::get_legend(p_pd_m)
plot <-
  plot_grid(p_ad_f, p_mci_f, p_pd_f, p_ad_m, p_mci_m, p_pd_m, nrow = 2)


ggsave(
  "figures/figure_2a.svg",
  plot,
  width = 9,
  height = 6,
  unit = "cm"
)
ggsave(
  "figures/figure_2a.png",
  plot,
  width = 9,
  height = 6,
  unit = "cm"
)

pdf("figures/figure_2_a_legend.pdf")
grid.newpage()
grid.draw(legend)
dev.off()

write.csv(p_ad_f$data, "SourceData/Figure 2a_ad_f.csv")
write.csv(p_ad_m$data, "SourceData/Figure 2a_ad_m.csv")
write.csv(p_pd_f$data, "SourceData/Figure 2a_pd_f.csv")
write.csv(p_pd_m$data, "SourceData/Figure 2a_pd_m.csv")
write.csv(p_mci_f$data, "SourceData/Figure 2a_mci_f.csv")
write.csv(p_mci_m$data, "SourceData/Figure 2a_mci_m.csv")
