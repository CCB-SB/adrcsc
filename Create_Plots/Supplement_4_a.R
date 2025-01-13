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
#library(svglite)
source("ADRC_theme.R")
source("Pipeline_All/src/helper.R")

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



pbmc = readRDS("/local/s8frgran/ADRC/CompleteObjectAnnotated_NewAnnotation_onlyFirstVisit.rds")
tbl = fread("Pipeline_All/results/metadata.filtered.csv")
tbl[Diagnosis == "Parkinson's Disease only", Diagnosis:="Parkinson's Disease"]
print("here1")
colors = fread("Pipeline_All/data/colors.csv", strip.white = F)
color_v = colors$Color
names(color_v) = colors$ID

name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease"="PD", "Parkinson's Disease only"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")

pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]
pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels=diagnosis_order)

print("here2")



umap_embedding_df = as.data.table(pbmc@reductions$umap@cell.embeddings)
umap_embedding_df$biogroup = factor(pbmc$biogroup_short, levels=diagnosis_order)
umap_embedding_df$Sample = pbmc$Sample

print("here3")

density_umap_per_biogroup = ggplot(umap_embedding_df, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point_rast(size=0.05, alpha=0.25) + 
  stat_density_2d(aes(fill=stat(nlevel)), geom="polygon", n=200, size=0.5) +
  facet_grid(. ~ biogroup) + scale_fill_viridis_c() +
  theme_adrc() + xlab("") + ylab("") +
  theme(legend.position = "none") + theme(legend.position="none", 
                                          #strip.text=element_text(size=7), 
                                          #text = element_text(size=7), 
                                          axis.text.x = element_blank(), axis.text.y = element_blank(),aspect.ratio = 1, axis.line=element_blank(), axis.ticks=element_blank())
print("here3.5")
figure1_f <- as.ggplot(fill_title(density_umap_per_biogroup, color_v))
print("here4")
#pdf()
#figure1_f
#dev.off()
save_plot("figures/figure4_a.pdf", fill_title(density_umap_per_biogroup, color_v), base_height = 50, base_width=120, units="mm")

