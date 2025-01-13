
source("ADRC_theme.R")

colors = fread("data/colors.csv", strip.white = F)
color_v = colors$Color
names(color_v) = colors$ID

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) = names(color)
  col
}


diagnosis_order = c("HC", "MCI", "AD", "PD-MCI", "PD")
celltype_order = c("CD4+ T-Helper Cell","Naive CD4+ T cell","CD4+ Memory T cell",
                   "Gamma delta T cell","Mucosal associated invariant T cell",
                   "CD8+ Memory T cell","Naive CD8+ T cell","Proliferating CD8+ T cell",
                   "Memory B cell","Double negative B cell","Naive B cell",
                   "B cell","Transitional B cell","Conventional Dendritic cell",
                   "CD16+ Monocyte","Megacaryocytes", "Proliferating Monocyte",
                   "CD14+ Monocyte","Plasmacytoid Dendritic cell","HSPC",
                   "Red Blood Cells","NKT-like cell","Proliferating CD4+ T cell",
                   "Plasmablasts","CD56-Dim NK cell", "CD56-Bright NK cell")
                   
CellCluster_order <- c("CD4+ T cell", "CD8+ T cell", "NK cell", "NKT-like cell", "Gamma delta T cell", "Mucosal associated invariant T cell", "B cell", "Plasma cell", "Monocyte", "Megacaryocytes", "Dendritic cell", "HSPC/Progenitors","Proliferating T cell","Red Blood Cells")


celltype_cluster_order = c("Naive CD4+ T cell","Treg CD4+ cell","CD4+ Tfh cell","Th1/Th2/Th17 cell","CD4+ Tcm/Tscm T cell",
"CD4+ Tcm T cell","CD4+ Tem cell","CD4+ Proliferating","Naive CD8+ T cell","CD8+ T cell","CD8+ Tcm T cell","CD8+ Tem T cell",
"CD8+ Tte T cell","Proliferating T cell","NKT-like cell","CD56-Dim, CD16 NK cell","CD56-Bright NK cell","Naive B cell",
"C-Memory B cell","M-Memory B cell","Transitional B cell","DN1 B cell","DN B cell","Plasmablasts","pDC cell","cDC1","cDC2",
"CD14 Monocytes","CD16 Monocytes","Proliferating Monocytes","Megacaryocytes","HSPC / Progenitor TO_DELETE","Red Blood Cells")

name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease only"="PD",
               "Parkinson's Disease"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")


fill_title = function(p, palette){
  g <- ggplot_gtable(ggplot_build(p))
  
  strips <- which(grepl('strip-', g$layout$name))
  
  for (i in seq_along(strips)) {
    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette[g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$label]
    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white"
  }
  return(g)
}


sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}

