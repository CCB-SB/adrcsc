library(circlize)
source("ADRC_theme.R")
source("src/helper.R")

df_to_matrix <- function(df, x_col, y_col, data_col){
  data <- sapply(unique(df[[x_col]]), function(x) sapply(unique(df[[y_col]]), function(y) {
    ret <- as.numeric(as.character(df[[data_col]][df[[x_col]]==x & df[[y_col]]==y][1]))
    if(is.na(ret)[1]){ret = 0}
    ret
  }
  ))
  
  return(data)
}

plot_heatmap <- function(dtp, lim = 20){
  df <- df_to_matrix(dtp, "X.Name", "cluster_id", "color")
  mat <- t(as.matrix(df))
  mat <- mat[rowSums(mat != 0)  > 3 , colSums(mat != 0) > 4]
  
  colors = colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027")) # black, red, green, blue

  plot <- ComplexHeatmap::Heatmap(t(mat), name = "mat",  row_names_side = "left",col = colors, 
                                  column_names_rot = 45,
                                  #left_annotation = rowAnnotation(celltype = colnames(mat),show_legend = c( FALSE), annotation_name_side = "top"),
                                  column_km = 3, row_km =3,
                                  show_row_dend = FALSE, show_column_dend = FALSE)
  return(plot)
}

pathways_female <- read.csv("../Upload/Pipeline_Female/results/pathway_analysis_results_all.csv")
pathways_female$sex <- "female"
pathways_male <- read.csv("../Upload/Pipeline_Male/results/pathway_analysis_results_All.csv")
pathways_male$sex <- "male"

pw_collected <- rbind(pathways_female, pathways_male)

# dtp <- pw_collected[pw_collected$contrast == "ADvsHC" & pw_collected$sex == "male",]
# dtp$Regulation_direction[dtp$Regulation_direction == 0] <- -1
# dtp$P.value[dtp$P.value< 10^-10] <- 10^-10
# dtp$color <- -log10(dtp$P.value) * as.numeric(dtp$Regulation_direction)
# plot <- plot_heatmap(dtp, 10)
# plot
# 
# dtp <- pw_collected[pw_collected$contrast == "ADvsHC" & pw_collected$sex == "female",]
# dtp$Regulation_direction[dtp$Regulation_direction == 0] <- -1
# dtp$P.value[dtp$P.value< 10^-10] <- 10^-10
# dtp$color <- -log10(dtp$P.value) * as.numeric(dtp$Regulation_direction)
# plot <- plot_heatmap(dtp, 10)
# plot

##############################

# plot_pathways_comparison <- function(comp, pw_count, label = T){
#   pw_count_gr1 <- pw_count[pw_count$contrast == comp & pw_count$sex == "female",]
#   pw_count_gr2 <- pw_count[pw_count$contrast == comp & pw_count$sex == "male",]
#   
#   data_to_plot <- merge(pw_count_gr1, pw_count_gr2, by = "X.Name", all = T)
#   data_to_plot$PW_count.x[is.na(data_to_plot$PW_count.x)] <- 0
#   data_to_plot$PW_count.y[is.na(data_to_plot$PW_count.y)] <- 0
#   data_to_plot$in_group <- "female and male"
#   data_to_plot$in_group[(data_to_plot$PW_count.x==0)] <- "only male"
#   data_to_plot$in_group[(data_to_plot$PW_count.y == 0)] <- "only female"
#   counts <- data_to_plot %>% group_by(PW_count.x, PW_count.y, in_group) %>% summarise(num_pathways = n())
#   counts$label <- sapply(1:length(counts$PW_count.x), function(i) paste(data_to_plot$X.Name[data_to_plot$PW_count.x == counts$PW_count.x[i] & data_to_plot$PW_count.y == counts$PW_count.y[i]],collapse = '\n'))
#   counts$label[counts$PW_count.x<2 & counts$PW_count.y<7] <- ""
#   counts$label[is.na(counts$label)] <- ""
#   counts$num_pathways[counts$num_pathways > 20] <- 20
#   if (!label) {counts$label <- ""}
#   plot <- ggplot(counts, aes(x = PW_count.x, y = PW_count.y, size = num_pathways, label = label, color = in_group)) + 
#     geom_point()+ 
#     geom_text_repel(size = 1)+ 
#     geom_abline(slope = 1, intercept = 0, linetype=3)+ 
#     theme_adrc()+ 
#     theme(aspect.ratio = 1)+ 
#     xlab("Number occurences\nin female")+ 
#     ylab("Number occurences\nin male")+
#     scale_size(range = c(.1,1), breaks = c(1,2, 5,10, 20) )  
#     #scale_x_continuous(limits = c(0, max(counts$PW_count.x, counts$PW_count.y)))+ 
#     #scale_y_continuous(limits = c(0, max(counts$PW_count.x, counts$PW_count.y)))
#   return(plot)
# }

# pw_count <- pw_collected %>% group_by(X.Name, contrast, sex) %>% summarise(PW_count = n())
# 
# colors <- c("#5D7322", "#9C489A", "black")
# names(colors) <- c("only male", "only female", "female and male")
# # plot1 <- plot_pathways_comparison("ADvsHC", pw_count, label = F)+ scale_color_manual(values = colors) + 
# #   labs(color="Significant in", size = "#Pathways")+ 
# #   theme(plot.margin = margin(0,1,0,0, unit = "cm"), legend.position = "none")
# # 
# # ggsave("figures/Comparison_Pathways_MF_AD.svg", plot1, width = 5, height = 5, unit = "cm")
# # 
# # plot2 <- plot_pathways_comparison("PDvsHC", pw_count, label = F)+ scale_color_manual(values = colors) + facet_wrap(.~contrast)+
# #   labs(color="Significant in", size = "#Pathways")+ 
# #   theme(plot.margin = margin(0,1,0,0, unit = "cm"), legend.position = "none")
# # ggsave("figures/Comparison_Pathways_MF_PD.svg", plot2, width = 5, height = 5, unit = "cm")
# # 
# # plot3 <- plot_pathways_comparison("MCIvsHC", pw_count, label = F)+ 
# #   scale_color_manual(values = colors) + labs(color="Significant in", size = "#Pathways")+
# #   theme(plot.margin = margin(0,1,0,0, unit = "cm"))
# # ggsave("figures/Comparison_Pathways_MF_MCI.svg", plot3, width = 5, height = 5, unit = "cm")
# # 
# # plot1 <- plot1 + theme(legend.position = "none")
# # plot2 <- plot2 + theme(legend.position = "none")
# # #plot <- ggarrange(plot1, plot2, plot3,  nrow = 1, width = c(1,1,2))
# # plot <- plot1+plot2+plot3
# # ggsave("figures/Comparison_Pathways_MF.svg", plot, width = 18, height = 6, unit = "cm")
# # 
# 
# # pw_count_gr1 <- pw_count[pw_count$sex == "female",]
# # pw_count_gr2 <- pw_count[pw_count$sex == "male",]
# # 
# # data_to_plot <- merge(pw_count_gr1, pw_count_gr2, by = c("X.Name", "contrast"), all = T)
# # data_to_plot$PW_count.x[is.na(data_to_plot$PW_count.x)] <- 0
# # data_to_plot$PW_count.y[is.na(data_to_plot$PW_count.y)] <- 0
# # data_to_plot$in_group <- "female and male"
# # data_to_plot$in_group[(data_to_plot$PW_count.x==0)] <- "only male"
# # data_to_plot$in_group[(data_to_plot$PW_count.y == 0)] <- "only female"
# # counts <- data_to_plot %>% group_by(PW_count.x, PW_count.y, in_group, contrast) %>% summarise(num_pathways = n())
# # counts$label <- sapply(1:length(counts$PW_count.x), function(i) paste(data_to_plot$X.Name[data_to_plot$PW_count.x == counts$PW_count.x[i] & data_to_plot$PW_count.y == counts$PW_count.y[i]],collapse = '\n'))
# # counts$label[counts$PW_count.x<2 & counts$PW_count.y<7] <- ""
# # counts$label[is.na(counts$label)] <- ""
# # counts$num_pathways[counts$num_pathways > 20] <- 20
# # counts$label <- ""
# # plot <- ggplot(counts, aes(x = PW_count.x, y = PW_count.y, size = num_pathways, label = label, color = in_group)) + 
# #   facet_grid2(.~contrast, scales = "free", independent = "y")+
# #   scale_color_manual(values = colors) + labs(color="Significant in", size = "#Occurences")+
# #   geom_point()+ 
# #   geom_text_repel(size = 1)+ 
# #   geom_abline(slope = 1, intercept = 0, linetype=3)+ 
# #   theme_adrc()+ 
# #   theme(aspect.ratio = 1)+ 
# #   xlab("Number occurences\nin female")+ 
# #   ylab("Number occurences\nin male")+
# #   scale_size(range = c(.1,1), breaks = c(1,2, 5,10, 20) )  
# # #scale_x_continuous(limits = c(0, max(counts$PW_count.x, counts$PW_count.y)))+ 
# # #scale_y_continuous(limits = c(0, max(counts$PW_count.x, counts$PW_count.y)))
# # plot <-  fill_title(plot, color_v)
# # ggsave("figures/Comparison_Pathways_MF.svg", plot, width = 18, height = 7, unit = "cm")
# 
# 
# 
# pw_tmp <- pw_count[pw_count$contrast == "ADvsHC"& pw_count$sex == "female",]
# ad_pw_f <- pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
# pw_tmp <- pw_count[pw_count$contrast == "ADvsHC"& pw_count$sex == "male",]
# ad_pw_m <- pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
# pw_tmp <- pw_count[pw_count$contrast == "PDvsHC"& pw_count$sex == "female",]
# pd_pw_f <- pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
# pw_tmp <- pw_count[pw_count$contrast == "PDvsHC"& pw_count$sex == "male",]
# pd_pw_m <- pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
# pw_tmp <- pw_count[pw_count$contrast == "MCIvsHC"& pw_count$sex == "female",]
# mci_pw_f <- pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
# pw_tmp <- pw_count[pw_count$contrast == "MCIvsHC"& pw_count$sex == "male",]
# mci_pw_m <- pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
# 
# pw_col <- c(ad_pw_f, ad_pw_m, pd_pw_f, pd_pw_m, mci_pw_f, mci_pw_m)
# df_pw <- data.frame(X.Name = pw_col, contrast = c(rep("ADvsHC",10), rep("PDvsHC", 10), rep("MCIvsHC",10)), sex = rep(c(rep("female",5), rep("male", 5)), 3))
# df_pw$amongTop5 <- T
# df <- pw_collected[pw_collected$X.Name %in% df_pw$X.Name,]
# summary <- df %>% group_by(X.Name, contrast, sex, Regulation_direction) %>% summarise(n = n())
# merged  <- merge(df_pw, summary, all = T, by = c("X.Name", "contrast", "sex"))
# merged$amongTop5[is.na(merged$amongTop5)] <- F
# 
# merged$X.Name <- unlist(lapply(as.character(merged$X.Name), function(string) {
#   n = unlist(gregexpr(pattern =' ', string))
# 
#   if (any(n>50)) {
#     pos <- n[n>50][1]
#     substr(string, pos, pos) <- "\n"
#   }
# 
#   if (any(n>25)) {
#     pos <- n[n>25][1]
#     substr(string, pos, pos) <- "\n"
#   }
#   
#   string
# }))
# merged$n_orig <- merged$n
# #merged$n[merged$n > 10] <- 10
# merged$n <- ifelse(merged$Regulation_direction, merged$n, -merged$n)
# merged$n[merged$n == 0] <- NA
# plot <- ggplot(merged, aes(x= as.character(Regulation_direction), y = X.Name,  fill = sign(n)* log10(abs(n)),label = ifelse(abs(merged$n_orig)>10, merged$n_orig, "")) ) + 
#   geom_point(shape = ifelse(merged$Regulation_direction==1,95,43 ))  +
#   geom_tile()  +
#   scale_x_discrete(labels = c("Enr.", "Depl."),breaks = c("1", "0"))+
#                    facet_nested(.~contrast+ sex) + 
#   scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", midpoint = 0, name = "#Occurences", breaks = c(-log10(5),-log10(2), 0, log10(2), log10(5), log10(10), log10(20), log10(30)), labels =  c(5,2, 0,2, 5,10, 20, 30))+
#   theme_adrc()+xlab("")+ ylab("")+ theme(axis.text.y = element_text(lineheight = 0.7), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# plot <-  fill_title(plot, color_v)
# ggsave("figures/Details_Pathways.svg", plot, width = 12, height = 9, unit = "cm")
# 

# 
# pw_collected_f <- pw_collected[pw_collected$sex == "female",]
# pw_upset <- list(PDvsHC = pw_collected_f$X.Name[pw_collected_f$contrast == "PDvsHC"], ADvsHC = pw_collected_f$X.Name[pw_collected_f$contrast == "ADvsHC"],  MCIvsHC = pw_collected_f$X.Name[pw_collected_f$contrast == "MCIvsHC"])
# 
# m1 = make_comb_mat(pw_upset)
# 
# top_annotation = upset_top_annotation(m1,
#                                       axis_param = list(gp = gpar(fontsize = 6)), 
#                                       gp = gpar(fontsize = 8, fill = "black"), annotation_name_gp = grid::gpar(fontsize = 8), height = unit(1.5, "cm")
# )
# right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
# 
# plot <- UpSet(m1, comb_order = order(comb_size(m1)),   top_annotation= top_annotation, right_annotation = right_annotation, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8))
# figure_e2 <- as.ggplot(plot)
# ggsave("figures/UpsetPathways_f.svg",figure_e2 , height = 4, width = 7, unit = "cm")
# 
# 
# pw_collected_m <- pw_collected[pw_collected$sex == "male",]
# pw_upset <- list(PDvsHC = pw_collected_m$X.Name[pw_collected_m$contrast == "PDvsHC"], ADvsHC = pw_collected_m$X.Name[pw_collected_m$contrast == "ADvsHC"],  MCIvsHC = pw_collected_m$X.Name[pw_collected_m$contrast == "MCIvsHC"])
# 
# m1 = make_comb_mat(pw_upset)
# 
# top_annotation = upset_top_annotation(m1,
#                                       axis_param = list(gp = gpar(fontsize = 6)), 
#                                       gp = gpar(fontsize = 8, fill = "black"), annotation_name_gp = grid::gpar(fontsize = 8), height = unit(1.5, "cm")
# )
# right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
# 
# plot <- UpSet(m1, comb_order = order(comb_size(m1)),   top_annotation= top_annotation, right_annotation = right_annotation, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8))
# figure_e2 <- as.ggplot(plot)
# ggsave("figures/UpsetPathways_m.svg",figure_e2  , height = 4, width = 7, unit = "cm")

##### Cluster the pathways and get most frequent per group

pw_subset <- pw_collected[pw_collected$contrast == "ADvsHC" & pw_collected$sex == "female",]
pw_upset <- lapply(unique(pw_subset$cluster_id), function(cluster) pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
plot <- UpSet(m1, right_annotation = right_annotation, comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_b.svg", width = 8, height = 5)
plot
dev.off()

pw_subset <- pw_collected[pw_collected$contrast == "ADvsHC" & pw_collected$sex == "male",]
pw_upset <- lapply(unique(pw_subset$cluster_id), function(cluster) pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
plot <- UpSet(m1, right_annotation = right_annotation, comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_a.svg", width = 8, height = 5)
plot
dev.off()

pw_subset <- pw_collected[pw_collected$contrast == "PDvsHC" & pw_collected$sex == "female",]
pw_upset <- lapply(unique(pw_subset$cluster_id), function(cluster) pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
plot <- UpSet(m1, right_annotation = right_annotation, comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_d.svg", width = 8, height = 5)
plot
dev.off()

pw_subset <- pw_collected[pw_collected$contrast == "PDvsHC" & pw_collected$sex == "male",]
pw_upset <- lapply(unique(pw_subset$cluster_id), function(cluster) pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
plot <- UpSet(m1, right_annotation = right_annotation, comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_c.svg", width = 8, height = 5)
plot
dev.off()
pw_subset <- pw_collected[pw_collected$contrast == "MCIvsHC" & pw_collected$sex == "female",]
pw_upset <- lapply(unique(pw_subset$cluster_id), function(cluster) pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
plot <- UpSet(m1, right_annotation = right_annotation, comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_f.svg", width = 8, height = 5)
plot
dev.off()

pw_subset <- pw_collected[pw_collected$contrast == "MCIvsHC" & pw_collected$sex == "male",]
pw_upset <- lapply(unique(pw_subset$cluster_id), function(cluster) pw_subset$X.Name[pw_subset$cluster_id == cluster])
names(pw_upset) <- unique(pw_subset$cluster_id)
m1 = make_comb_mat(pw_upset)
m1 <- m1[comb_degree(m1) > 1 & comb_size(m1) > 1]
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))
plot <- UpSet(m1, right_annotation = right_annotation, comb_order = order(comb_size(m1)))
svg("figures/Supplement_7_e.svg", width = 8, height = 5)
plot
dev.off()


