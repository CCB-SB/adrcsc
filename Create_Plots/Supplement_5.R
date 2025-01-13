library(dplyr)
library(ggplot2)
library(ggsankey)
library(ggrepel)
library(ggplotify)
library(lsa)
library(ComplexHeatmap)

degs_female <- read.csv("/local/s8frgran/ADRC/Upload/Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep = "\t")
degs_female$Sex <- "female"
degs_male <- read.csv("/local/s8frgran/ADRC/Upload/Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep = "\t")
degs_male$Sex <- "male"

degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "ADvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "female",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")

getSimilarity <- function(data, cell_column){
  degs_matrix <- NA
  for (celltype in unique(data[[cell_column]])){
    print(celltype)
    degs_cellt <- data[data[[cell_column]] == celltype,colnames(data) %in% c("logFC", "gene")]
    
    degs_cellt$rank <- NA
    degs_cellt$rank[order((degs_cellt$logFC), decreasing = T)] <- 1:nrow(degs_cellt)
    
    degs_cellt <- degs_cellt[,colnames(degs_cellt) %in% c("rank", "gene")]
    if (is.na(is.na(degs_matrix)[1])| is.na(degs_matrix)[1]){
      degs_matrix <- degs_cellt
      colnames(degs_matrix)[colnames(degs_matrix) == "rank"] <- celltype
    } else {
      degs_matrix <- merge(degs_matrix, degs_cellt, by = "gene", all = T)
      colnames(degs_matrix)[colnames(degs_matrix) == "rank"] <- celltype
    }
  }

  degs_matrix[is.na(degs_matrix)] <- 0
  degs_matrix <- as.matrix(degs_matrix[,-1])
  similarity <- cosine(degs_matrix)
  similarity_without_mid <- similarity
  similarity_without_mid[similarity_without_mid == 1] <- 0
  similarity <- similarity[abs(apply(similarity_without_mid, 1, max, na.rm=TRUE))>0.5,abs(apply(similarity_without_mid, 1, max, na.rm=TRUE))>0.5]
  
  return(similarity)
  
}



colors = fread("colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

library(circlize)
color_scale = colorRamp2(seq(0, 1, length = 3), c("blue", "#EEEEEE", "red"))

degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "ADvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "female",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")

similarity <- getSimilarity(degs_collected, "Celltype_sex")

column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
plot1 <- Heatmap(similarity, name = "mat", 
                col = color_scale, 
                top_annotation = column_ha, right_annotation = row_ha, 
                show_column_names =FALSE, show_row_names =FALSE, 
                column_title = NULL, 
                row_title = NULL, 
                row_km = 5, column_km = 5, 
                show_column_dend = F, show_row_dend = F, 
                show_heatmap_legend = F)
ggsave("figures/Supplement_5_AD_f.svg", as.ggplot(plot1), width = 4, height = 4)
write.csv(similarity, "SourceData/Supplement_5_AD_f.csv")

degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "ADvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "male",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")
similarity <- getSimilarity(degs_collected, "Celltype_sex")

column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
plot2 <- Heatmap(similarity, name = "mat", 
                col = color_scale, 
                top_annotation = column_ha, right_annotation = row_ha, 
                show_column_names =FALSE, show_row_names =FALSE, 
                column_title = NULL, 
                row_title = NULL, 
                row_km = 5, column_km = 5, 
                show_column_dend = F, show_row_dend = F, 
                show_heatmap_legend = F)
ggsave("figures/Supplement_5_AD_m.svg", as.ggplot(plot2), width = 4, height = 4)
write.csv(similarity, "SourceData/Supplement_5_AD_m.csv")
# 
# degs_collected <- rbind(degs_female, degs_male)
# degs_collected <- degs_collected[degs_collected$contrast == "ADvsHC",]
# degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")
# similarity <- getSimilarity(degs_collected, "Celltype_sex")
# similarity <- similarity[str_replace_all(colnames(similarity), ".*~", "") == "female", str_replace_all(rownames(similarity), ".*~", "") == "male"]
# 
# column_ha = HeatmapAnnotation(CellGroup = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(2, "cm"), col = color_v), show_legend = FALSE)
# row_ha = rowAnnotation(CellGroup = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(2, "cm"), col = color_v), show_legend = FALSE)
# plot <- Heatmap(similarity, name = "mat", 
#                 col = color_scale, 
#                 top_annotation = column_ha, right_annotation = row_ha, 
#                 show_column_names =FALSE, show_row_names =FALSE, 
#                 row_km = 5, column_km = 5, 
#                 show_column_dend = F, show_row_dend = F, 
#                 show_heatmap_legend = T)
# ggsave("figures/Supplement_5_AD.svg", as.ggplot(plot), width = 3, height = 3)
# 

degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "PDvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "female",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")

similarity <- getSimilarity(degs_collected, "Celltype_sex")

column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
plot3 <- Heatmap(similarity, name = "mat", 
                col = color_scale, 
                top_annotation = column_ha, right_annotation = row_ha, 
                show_column_names =FALSE, show_row_names =FALSE, 
                column_title = NULL, 
                row_title = NULL, 
                row_km = 5, column_km = 5, 
                show_column_dend = F, show_row_dend = F, 
                show_heatmap_legend = F)
ggsave("figures/Supplement_5_PD_f.svg", as.ggplot(plot3), width = 4, height = 4)
write.csv(similarity, "SourceData/Supplement_5_PD_f.csv")

degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "PDvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "male",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")
similarity <- getSimilarity(degs_collected, "Celltype_sex")

column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
plot4 <- Heatmap(similarity, name = "mat", 
                col = color_scale, 
                top_annotation = column_ha, right_annotation = row_ha, 
                show_column_names =FALSE, show_row_names =FALSE, 
                column_title = NULL, 
                row_title = NULL, 
                row_km = 5, column_km = 5, 
                show_column_dend = F, show_row_dend = F, 
                show_heatmap_legend = F)
ggsave("figures/Supplement_5_PD_m.svg", as.ggplot(plot4), width = 4, height = 4)
write.csv(similarity, "SourceData/Supplement_5_PD_m.csv")
# 
# 
# degs_collected <- rbind(degs_female, degs_male)
# degs_collected <- degs_collected[degs_collected$contrast == "PDvsHC",]
# degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")
# similarity <- getSimilarity(degs_collected, "Celltype_sex")
# similarity <- similarity[str_replace_all(colnames(similarity), ".*~", "") == "female", str_replace_all(rownames(similarity), ".*~", "") == "male"]
# 
# column_ha = HeatmapAnnotation(CellGroup = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
# row_ha = rowAnnotation(CellGroup = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
# plot <- Heatmap(similarity, name = "mat", 
#                 col = color_scale, 
#                 top_annotation = column_ha, right_annotation = row_ha, 
#                 show_column_names =FALSE, show_row_names =FALSE, 
#                 row_km = 5, column_km = 5, 
#                 show_column_dend = F, show_row_dend = F, 
#                 show_heatmap_legend = T)
# ggsave("figures/Heatmal_Celltype_Similarity_PD_mvsf.svg", as.ggplot(plot), width = 4, height = 4)



degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "MCIvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "female",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")

similarity <- getSimilarity(degs_collected, "Celltype_sex")

column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = T)
row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = T)
plot5 <- Heatmap(similarity, name = "mat", 
                col = color_scale, 
                top_annotation = column_ha, right_annotation = row_ha, 
                show_column_names =F, show_row_names =F, 
                column_title = NULL, 
                row_title = NULL, 
                row_km = 5, column_km = 5, 
                show_column_dend = F, show_row_dend = F, 
                show_heatmap_legend = F)
ggsave("figures/Supplement_5_MCI_f.svg", as.ggplot(plot5), width = 5, height = 4)
write.csv(similarity, "SourceData/Supplement_5_MCI_f.csv")

degs_collected <- rbind(degs_female, degs_male)
degs_collected <- degs_collected[degs_collected$contrast == "MCIvsHC",]
degs_collected <- degs_collected[degs_collected$Sex == "male",]
degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")
similarity <- getSimilarity(degs_collected, "Celltype_sex")

column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
plot6 <- Heatmap(similarity, name = "mat", 
                col = color_scale, 
                top_annotation = column_ha, right_annotation = row_ha, 
                show_column_names =FALSE, show_row_names =FALSE, 
                column_title = NULL, 
                row_title = NULL, 
                row_km = 5, column_km = 5, 
                show_column_dend = F, show_row_dend = F, 
                show_heatmap_legend = F)

ggsave("figures/Supplement_5_MCI_m.svg", as.ggplot(plot6), width = 4, height = 4)
write.csv(similarity, "SourceData/Supplement_5_MCI_m.csv")
lgd = (Legend(col_fun = color_scale, title = "cosine\nsimilarity"))


plot <- ggarrange(as.ggplot(plot1) + as.ggplot(plot3) + as.ggplot(plot5),
          as.ggplot(plot2) + as.ggplot(plot4) + as.ggplot(plot6), nrow = 2)
ggsave("figures/Supplement_5.svg", plot, width = 12, height = 8)
svg("figures/Legend_Supplement_5.svg", width = 1, height = 1.5)
draw(lgd)
dev.off()
# 
# degs_collected <- rbind(degs_female, degs_male)
# degs_collected <- degs_collected[degs_collected$contrast == "MCIvsHC",]
# degs_collected$Celltype_sex <- paste(degs_collected$cluster_id, degs_collected$Sex, sep = "~")
# similarity <- getSimilarity(degs_collected, "Celltype_sex")
# similarity <- similarity[str_replace_all(colnames(similarity), ".*~", "") == "female", str_replace_all(rownames(similarity), ".*~", "") == "male"]
# 
# column_ha = HeatmapAnnotation(Type = anno_simple(str_replace_all(colnames(similarity), "~.*", ""), height = unit(1, "cm"), col = color_v), show_legend = FALSE)
# row_ha = rowAnnotation(Type = anno_simple(str_replace_all(rownames(similarity), "~.*", ""), width = unit(1, "cm"), col = color_v), show_legend = FALSE)
# plot <- Heatmap(similarity, name = "mat", 
#                 col = color_scale, 
#                 top_annotation = column_ha, right_annotation = row_ha, 
#                 show_column_names =FALSE, show_row_names =FALSE, 
#                 column_title = NULL, 
#                 row_title = NULL, 
#                 row_km = 5, column_km = 5, 
#                 show_column_dend = F, show_row_dend = F, 
#                 show_heatmap_legend = T)
# ggsave("figures/Heatmal_Celltype_Similarity_MCI_mvsf.svg", as.ggplot(plot), width = 4, height = 4)



