library(circlize)
library(data.table)
library(ComplexHeatmap)

source("scripts/helper.R")

df_to_matrix <- function(df, x_col, y_col, data_col) {
  data <-
    sapply(unique(df[[x_col]]), function(x)
      sapply(unique(df[[y_col]]), function(y) {
        ret <-
          as.numeric(as.character(df[[data_col]][df[[x_col]] == x &
                                                   df[[y_col]] == y][1]))
        if (is.na(ret)[1]) {
          ret = 0
        }
        ret
      }))
  
  return(data)
}

plot_heatmap <- function(dtp, lim = 20) {
  df <- df_to_matrix(dtp, "X.Name", "cluster_id", "color")
  mat <- t(as.matrix(df))
  mat <- mat[rowSums(mat != 0)  > 3 , colSums(mat != 0) > 4]
  
  colors = colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))
  
  plot <-
    ComplexHeatmap::Heatmap(
      t(mat),
      name = "mat",
      row_names_side = "left",
      col = colors,
      column_names_rot = 45,
      column_km = 3,
      row_km = 3,
      show_row_dend = FALSE,
      show_column_dend = FALSE
    )
  return(plot)
}

pathways_female <-
  read.csv("results_female/pathway_analysis_results_all.csv")
pathways_female$sex <- "female"
pathways_male <-
  read.csv("results_male/pathway_analysis_results_All.csv")
pathways_male$sex <- "male"

pw_collected <- rbind(pathways_female, pathways_male)

##### Cluster the pathways and get most frequent per group

pw_subset <-
  pw_collected[pw_collected$contrast == "ADvsHC" &
                 pw_collected$sex == "female", ]
pw_upset <-
  lapply(unique(pw_subset$cluster_id), function(cluster)
    pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)
plot <-
  UpSet(m1,
        right_annotation = right_annotation,
        comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_b.svg",
    width = 8,
    height = 5)
plot
dev.off()

pw_subset <-
  pw_collected[pw_collected$contrast == "ADvsHC" &
                 pw_collected$sex == "male", ]
pw_upset <-
  lapply(unique(pw_subset$cluster_id), function(cluster)
    pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)
plot <-
  UpSet(m1,
        right_annotation = right_annotation,
        comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_a.svg",
    width = 8,
    height = 5)
plot
dev.off()

pw_subset <-
  pw_collected[pw_collected$contrast == "PDvsHC" &
                 pw_collected$sex == "female", ]
pw_upset <-
  lapply(unique(pw_subset$cluster_id), function(cluster)
    pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)
plot <-
  UpSet(m1,
        right_annotation = right_annotation,
        comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_d.svg",
    width = 8,
    height = 5)
plot
dev.off()

pw_subset <-
  pw_collected[pw_collected$contrast == "PDvsHC" &
                 pw_collected$sex == "male", ]
pw_upset <-
  lapply(unique(pw_subset$cluster_id), function(cluster)
    pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)
plot <-
  UpSet(m1,
        right_annotation = right_annotation,
        comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_c.svg",
    width = 8,
    height = 5)
plot
dev.off()
pw_subset <-
  pw_collected[pw_collected$contrast == "MCIvsHC" &
                 pw_collected$sex == "female", ]
pw_upset <-
  lapply(unique(pw_subset$cluster_id), function(cluster)
    pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)
plot <-
  UpSet(m1,
        right_annotation = right_annotation,
        comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_f.svg",
    width = 8,
    height = 5)
plot
dev.off()

pw_subset <-
  pw_collected[pw_collected$contrast == "MCIvsHC" &
                 pw_collected$sex == "male", ]
pw_upset <-
  lapply(unique(pw_subset$cluster_id), function(cluster)
    pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(
  m1,
  gp = gpar(fill = color_v[rownames(m1)], fontsize = 2),
  annotation_name_gp = grid::gpar(fontsize = 8)
)
plot <-
  UpSet(m1,
        right_annotation = right_annotation,
        comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_e.svg",
    width = 8,
    height = 5)
plot
dev.off()
