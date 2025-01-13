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


degs <-
  read.csv(
    "results_male/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )

degs <- degs[degs$contrast %in% c("PDvsHC", "ADvsHC", "MCIvsHC"), ]
degs <- degs[degs$p_adj.loc < 0.05 & abs(degs$logFC) > 0.5, ]

degs_upset <-
  list(
    PDvsHC = degs$gene[degs$contrast == "PDvsHC"],
    ADvsHC = degs$gene[degs$contrast == "ADvsHC"],
    MCIvsHC = degs$gene[degs$contrast == "MCIvsHC"]
  )

m1 = make_comb_mat(degs_upset)


top_annotation = upset_top_annotation(
  m1,
  axis_param = list(gp = gpar(fontsize = 6)),
  gp = gpar(fontsize = 8, fill = "black"),
  annotation_name_gp = grid::gpar(fontsize = 8),
  height = unit(1.5, "cm")
)
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)

plot <-
  UpSet(
    m1,
    comb_order = order(comb_size(m1)),
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8)
  )
figure_e1 <- as.ggplot(plot)

ggsave("figures/figure_3_e.svg",
       figure_e1 ,
       height = 2,
       width = 3)


degs <-
  read.csv(
    "results_female/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )

degs <- degs[degs$contrast %in% c("PDvsHC", "ADvsHC", "MCIvsHC"), ]
degs <- degs[degs$p_adj.loc < 0.05 & abs(degs$logFC) > 0.5, ]

degs_upset <-
  list(
    PDvsHC = degs$gene[degs$contrast == "PDvsHC"],
    ADvsHC = degs$gene[degs$contrast == "ADvsHC"],
    MCIvsHC = degs$gene[degs$contrast == "MCIvsHC"]
  )

m1 = make_comb_mat(degs_upset)

top_annotation = upset_top_annotation(
  m1,
  axis_param = list(gp = gpar(fontsize = 6)),
  gp = gpar(fontsize = 8, fill = "black"),
  annotation_name_gp = grid::gpar(fontsize = 8),
  height = unit(1.5, "cm")
)
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)

plot <-
  UpSet(
    m1,
    comb_order = order(comb_size(m1)),
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8)
  )
figure_e2 <- as.ggplot(plot)
ggsave("figures/figure_3_f.svg",
       figure_e2 ,
       height = 2,
       width = 3)
