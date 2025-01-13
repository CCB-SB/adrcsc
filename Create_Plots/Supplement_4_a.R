library(ggplot2)
library(cowplot)
library(data.table)
library(ggsci)
library(ggforce)
library(readxl)
library(ggraph)
library(igraph)
library(gghalves)
library(Seurat)
library(reshape2)
library(ggsignif)
library(ggrastr)
library(UpSetR)
library(ComplexUpset)
library(pbapply)
library(ggpubr)
library(grid)
library(viridis)
library(aplot)
library(ggplotify)

source("scripts/helper.R")


pbmc = readRDS("data/CompleteObjectAnnotated_onlyFirstVisit.rds")
tbl = fread("results_all/metadata.filtered.csv")
tbl[Diagnosis == "Parkinson's Disease only", Diagnosis := "Parkinson's Disease"]

name2short = c(
  "All" = "All",
  "Neurodegeneration" = "",
  "Cognitive Impairment" = "",
  "Healthy Control" = "HC",
  "Parkinson's Disease" = "PD",
  "Parkinson's Disease only" = "PD",
  "Alzheimer's disease" = "AD",
  "Parkinson's Disease with MCI" = "PD-MCI",
  "Mild Cognitive Impairment" = "MCI"
)

pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]
pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels = diagnosis_order)

umap_embedding_df = as.data.table(pbmc@reductions$umap@cell.embeddings)
umap_embedding_df$biogroup = factor(pbmc$biogroup_short, levels = diagnosis_order)
umap_embedding_df$Sample = pbmc$Sample

density_umap_per_biogroup = ggplot(umap_embedding_df, aes(x = UMAP_1, y =
                                                            UMAP_2)) +
  geom_point_rast(size = 0.05, alpha = 0.25) +
  stat_density_2d(aes(fill = stat(nlevel)),
                  geom = "polygon",
                  n = 200,
                  size = 0.5) +
  facet_grid(. ~ biogroup) + scale_fill_viridis_c() +
  theme_adrc() + xlab("") + ylab("") +
  theme(legend.position = "none") + theme(
    legend.position = "none",
    #strip.text=element_text(size=7),
    #text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    aspect.ratio = 1,
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

figure1_f <-
  as.ggplot(fill_title(density_umap_per_biogroup, color_v))

save_plot(
  "figures/figure4_a.pdf",
  fill_title(density_umap_per_biogroup, color_v),
  base_height = 50,
  base_width = 120,
  units = "mm"
)
