library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggplotify)
library(data.table)
source("ADRC_theme.R")
source("scripts/helper.R")


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




degs_male <- read.csv("../Upload/Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_male$gender <- "male"
degs_female <- read.csv("../Upload/Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_female$gender <- "female"
degs <- rbind(degs_female, degs_male)



unique(degs$contrast)
degs_ad <- degs[degs$contrast %in% c("ADvsHC"),]#, "PDMCIvsHC", "PDvsMCI", "PDvsHC", "ADvsPD"),]
genes <- as.data.frame(degs_ad$gene)

Alzheimer_kegg <- read.csv("external_data/AlzheimerGenes.csv", header=FALSE)

degs_kegg_ad <- degs_ad[gsub("^MT-", "", degs_ad$gene) %in% Alzheimer_kegg$V1 & degs_ad$contrast == "ADvsHC",]

freq <- as.data.frame(table(degs_kegg_ad$gene[degs_kegg_ad$p_adj.loc<0.1]))
use <- freq$Var1[freq$Freq >1]
degs_kegg_ad <- degs_kegg_ad[(degs_kegg_ad$gene %in% use) ,]

Cellt_to_use <- unique(degs_kegg_ad$cluster_id[degs_kegg_ad$p_adj.loc<0.05])
degs_kegg_ad <- degs_kegg_ad[ (degs_kegg_ad$cluster_id %in% Cellt_to_use),]

degs_kegg_ad$logFC[degs_kegg_ad$logFC>1.5] <- 1.5
degs_kegg_ad$logFC[degs_kegg_ad$logFC< (-1.5)] <- -1.5


female = sprintf(intToUtf8(9792))
male = intToUtf8(9794)

label_names <- list(
  'male'=male,
  'female'=female
)

sex_labeller <- function(variable,value){
  variable$gender <- unlist(lapply(variable$gender, function(x) ifelse(x=="female", female, male)))
  print((labels_cell[variable$cluster_id]))
  variable$cluster_id <- unlist(lapply(variable$cluster_id, function(x) as.character(labels_cell[x])))
  print(variable)                   
  return(variable)
}

p <- ggplot() +   geom_point(data=degs_kegg_ad, aes(x=cluster_id , y=gene,  size=-log10(p_adj.loc), fill=logFC, color = p_adj.loc<0.05 ), alpha = 0.8, shape = 21) +
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-1.5,1.5))+ 
  scale_color_manual(values=c("white","black"), guide = "none")+
  facet_nested(contrast~cluster_id+ gender, scales = "free_x",  labeller=sex_labeller,  remove_labels = "x")  + 
  scale_x_discrete(name = "", labels = element_blank()) +
  xlab("") + ylab("") + 
  theme_adrc() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(size = "log10 adj.\np-value")


color_v2 <- color_v
names(color_v2) <- labels_cell[names(color_v)]

female = sprintf(intToUtf8(9792))
male = intToUtf8(9794)

label_names <- c(color_v["Male"],color_v["Female"])
names(label_names) <- c(male, female)

color <- c(color_v, color_v2, label_names)
figure_g <- as.ggplot(fill_title(p, color))
save_plot("figures/supplement_6_a.svg", figure_g, base_height =5, base_width=7)#,device = cairo_pdf) 



degs_male <- read.csv("../Upload/Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_male$gender <- "male"
degs_female <- read.csv("../Upload/Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_female$gender <- "female"
degs <- rbind(degs_female, degs_male)


unique(degs$contrast)
degs_pd <- degs[degs$contrast %in% c("PDvsHC"),]#, "PDMCIvsHC", "PDvsMCI", "PDvsHC", "ADvsPD"),]
genes <- as.data.frame(degs_pd$gene)

Parkinson_kegg <- read.csv("external_data/ParkinsonGenes.csv", header=FALSE)

degs_kegg_pd <- degs_pd[gsub("^MT-", "", degs_pd$gene) %in% Parkinson_kegg$V1 & degs_pd$contrast == "PDvsHC",]

freq <- as.data.frame(table(degs_kegg_pd$gene[degs_kegg_pd$p_adj.loc<0.1]))
use <- freq$Var1[freq$Freq >1]
degs_kegg_pd <- degs_kegg_pd[(degs_kegg_pd$gene %in% use) ,]

Cellt_to_use <- unique(degs_kegg_pd$cluster_id[degs_kegg_pd$p_adj.loc<0.05])
degs_kegg_pd <- degs_kegg_pd[ (degs_kegg_pd$cluster_id %in% Cellt_to_use),]

degs_kegg_pd$logFC[degs_kegg_pd$logFC>1.5] <- 1.5
degs_kegg_pd$logFC[degs_kegg_pd$logFC< (-1.5)] <- -1.5


female = sprintf(intToUtf8(9792))
male = intToUtf8(9794)

label_names <- list(
  'male'=male,
  'female'=female
)

sex_labeller <- function(variable,value){
  variable$gender <- unlist(lapply(variable$gender, function(x) ifelse(x=="female", female, male)))
  print((labels_cell[variable$cluster_id]))
  variable$cluster_id <- unlist(lapply(variable$cluster_id, function(x) as.character(labels_cell[x])))
  print(variable)                   
  return(variable)
}

p <- ggplot() +   geom_point(data=degs_kegg_pd, aes(x=cluster_id , y=gene,  size=-log10(p_adj.loc), fill=logFC, color = p_adj.loc<0.05 ), alpha = 0.8, shape = 21) +
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-1.5,1.5))+ 
  scale_color_manual(values=c("white","black"), guide = "none")+
  facet_nested(contrast~cluster_id+ gender, scales = "free_x",  labeller=sex_labeller,  remove_labels = "x")  + 
  scale_x_discrete(name = "", labels = element_blank()) +
  xlab("") + ylab("") + 
  theme_adrc() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(size = "log10 adj.\np-value")

figure_g <- as.ggplot(fill_title(p, color))
save_plot("figures/supplement_6_b.svg", figure_g, base_height =5, base_width=7)#, device = cairo_pdf) 








