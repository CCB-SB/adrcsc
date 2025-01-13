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
library(tidyverse)
library(rstatix)

source("ADRC_theme.R")
source("scripts/helper.R")
source("scripts/Functions_Fig2.R")

pbmc <- readRDS("data/CompleteObjectAnnotated_onlyFirstVisit.rds")

pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]
pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels = diagnosis_order)


pbmc_fem <- subset(pbmc, subset = Sex == "female")
cells_per_sample_perc_df_female <-
  get_cells_per_sample_perc(pbmc_fem)

pbmc_male <- subset(pbmc, subset = Sex == "male")
cells_per_sample_perc_df_male <-
  get_cells_per_sample_perc(pbmc_male)

p1 <-
  as.ggplot(
    getCelltypeProportion_Boxplot_scaled(cells_per_sample_perc_df_male, c("CD4+ T cell"), y_lim = 1.20)
  )
p2 <-
  as.ggplot(
    getCelltypeProportion_Boxplot_scaled(
      cells_per_sample_perc_df_male,
      c("Monocyte", "Dendritic cell", "Megacaryocytes", "Plasma cell"),
      y_lim = 0.3
    )
  )
p3 <-
  as.ggplot(
    getCelltypeProportion_Boxplot_scaled(
      cells_per_sample_perc_df_male,
      c("CD8+ T cell", "NK cell", "B cell", "NKT-like cell"),
      y_lim = 0.6
    )
  )
p <-
  p1 %>% insert_right(p2, width = 3) %>% insert_right(p3, width = 3)
save_plot(
  "figures/suppl_4_a_male.svg",
  p,
  base_height = 40,
  base_width = 220,
  units = "mm"
)

# Significant differences in males
p_sign <-
  cells_per_sample_perc_df_male %>% group_by(L2) %>% pairwise_t_test(value ~ Diagnosis, p.adjust.method = "BH")
p_sign[p_sign$p.adj < 0.05 & p_sign$group1 == "HC",]


p1 <-
  as.ggplot(
    getCelltypeProportion_Boxplot_scaled(cells_per_sample_perc_df_female, c("CD4+ T cell"), y_lim = 1.20)
  )
p2 <-
  as.ggplot(
    getCelltypeProportion_Boxplot_scaled(
      cells_per_sample_perc_df_female,
      c("Monocyte", "Dendritic cell", "Megacaryocytes", "Plasma cell"),
      y_lim = 0.3
    )
  )
p3 <-
  as.ggplot(
    getCelltypeProportion_Boxplot_scaled(
      cells_per_sample_perc_df_female,
      c("CD8+ T cell", "NK cell", "B cell", "NKT-like cell"),
      y_lim = 0.6
    )
  )
p <-
  p1 %>% insert_right(p2, width = 3) %>% insert_right(p3, width = 3)
save_plot(
  "figures/suppl_4_a_female.svg",
  p,
  base_height = 40,
  base_width = 220,
  units = "mm"
)
# Significant differences in females
p_sign <-
  cells_per_sample_perc_df_female %>% group_by(L2) %>% pairwise_t_test(value ~ Diagnosis, p.adjust.method = "BH")
p_sign[p_sign$p.adj < 0.05 & p_sign$group1 == "HC",]
