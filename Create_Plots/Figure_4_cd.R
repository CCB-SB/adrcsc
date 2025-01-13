library(circlize)
library(ComplexHeatmap)
library(data.table)
library(ggplotify)
library(ggplot2)
source("scripts/helper.R")

pathways_female <-
  read.csv("results_female/pathway_analysis_results_all.csv")
pathways_female$sex <- "female"
pathways_male <-
  read.csv("results_male/pathway_analysis_results_All.csv")
pathways_male$sex <- "male"

pw_collected <- rbind(pathways_female, pathways_male)



pw_collected_f <- pw_collected[pw_collected$sex == "female", ]
pw_upset <-
  list(
    PDvsHC = pw_collected_f$X.Name[pw_collected_f$contrast == "PDvsHC"],
    ADvsHC = pw_collected_f$X.Name[pw_collected_f$contrast == "ADvsHC"],
    MCIvsHC = pw_collected_f$X.Name[pw_collected_f$contrast == "MCIvsHC"]
  )

m1 = make_comb_mat(pw_upset)

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
ggsave(
  "figures/figure_4_d.svg",
  figure_e2 ,
  height = 4,
  width = 7,
  unit = "cm"
)


pw_collected_m <- pw_collected[pw_collected$sex == "male", ]
pw_upset <-
  list(
    PDvsHC = pw_collected_m$X.Name[pw_collected_m$contrast == "PDvsHC"],
    ADvsHC = pw_collected_m$X.Name[pw_collected_m$contrast == "ADvsHC"],
    MCIvsHC = pw_collected_m$X.Name[pw_collected_m$contrast == "MCIvsHC"]
  )

m1 = make_comb_mat(pw_upset)

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
ggsave(
  "figures/figure_4_c.svg",
  figure_e2  ,
  height = 4,
  width = 7,
  unit = "cm"
)
