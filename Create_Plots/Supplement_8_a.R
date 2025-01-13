library(ggplot2)
library(cowplot)
source("src/helper.R")
source("ADRC_theme.R")
library(aplot)
library(ggpubr)
library(viridis)

path <- paste("BrainData/Brain_PFC_Colltected.csv")
degs_brain <- read.csv(path, sep = ",")
degs_brain$sex <- degs_brain$Sex
degs_brain$region <- "Prefrontal Cortex"
degs_brain$tissue <- "Rosmap"
degs_rosmap <- degs_brain[, colnames(degs_brain) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id", "sex", "region")]


path <- paste("BrainData/human_cortex_sex_condition_marker.csv")
degs_brain <- read.csv(path, sep = ",")
degs_brain <- degs_brain[degs_brain$contrast == "AD-CT",]
degs_brain$cluster_id <- degs_brain$cell_type

#degs_brain$tissue <- "brain"
degs_brain$logFC <- degs_brain$logFC
degs_brain$p_adj.loc <- degs_brain$FDR
degs_brain$gene <- degs_brain$gene
degs_brain$region <- "Cortex"
degs_brain$tissue <- "Zebra"
degs_zebra <- degs_brain[, colnames(degs_brain) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id", "sex", "region")]


read_pbmc_input <- function(path){
  degs_pbmc <- read.csv(path, sep = "\t")
  
  degs_pbmc <- degs_pbmc[degs_pbmc$contrast == "ADvsHC" ,] 
  degs_pbmc$tissue <- "pbmc"
  degs_pbmc <- degs_pbmc[, colnames(degs_pbmc) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
  return(degs_pbmc)
}

plotComparison <- function(brain_data, pbmc_data, rosmap_data, celltypes, celltype_to_order){
  
  degs <- rbind(brain_data, pbmc_data, rosmap_data)
  # print(degs[degs$gene %in% degs$gene[degs$tissue == "brain"] & degs$tissue == "pbmc",])
  # degs <- degs[degs$gene %in% degs$gene[degs$tissue == "brain"],]
  # degs <- degs[degs$gene %in% degs$gene[!degs$tissue %in% c("brain", "pbmc")],]
  # degs <- degs[degs$gene %in% degs$gene[degs$tissue == "pbmc"],]
  # degs <- degs[degs$cluster_id %in% celltypes | degs$tissue == "brain",]
  #degs <- degs[degs$gene %in% degs$gene[degs$cluster_id %in% celltypes],]
  degs <- degs[!is.na(degs$p_adj.loc),]
  gene_list_to_order <- unique(degs$gene)
  fc_to_order <- lapply(gene_list_to_order, function(gene) {
    if (gene %in% degs$genes[degs$cluster_id %in% celltype_to_order]) {degs$logFC[degs$cluster_id %in% celltype_to_order & degs$gene == gene]} else 
    {degs$logFC[degs$gene == gene][1]}
  })
  levels <- unique(degs$gene[degs$cluster_id %in% celltype_to_order][order(degs$logFC[degs$cluster_id %in% celltype_to_order], decreasing = T)])
  degs$gene <- factor(degs$gene, levels = levels)
  
  degs$tissue <- factor(degs$tissue, levels = c("pbmc", "Zebra", unique(degs$tissue[!degs$tissue %in% c("Zebra", "pbmc")])))
  
  degs$logFC[degs$logFC>0.5] <- 0.5
  degs$logFC[degs$logFC<(-0.5)] <- (-0.5)
  degs$p_adj.loc[degs$p_adj.loc == 0] <- min(degs$p_adj.loc[degs$p_adj.loc != 0])
  #print(head(degs))
  degs$p_adj.loc_norm <- NA
  degs$p_adj.loc_norm[degs$tissue == "pbmc"] <- -log(degs$p_adj.loc[degs$tissue == "pbmc"]) / max(-log(degs$p_adj.loc[degs$tissue == "pbmc"]))
  degs$p_adj.loc_norm[degs$tissue == "Zebra"] <- -log(degs$p_adj.loc[degs$tissue == "Zebra"]) / max(-log(degs$p_adj.loc[degs$tissue == "Zebra"]))
  degs$p_adj.loc_norm[degs$tissue == "Rosmap"] <- -log(degs$p_adj.loc[degs$tissue == "Rosmap"]) / max(-log(degs$p_adj.loc[degs$tissue == "Rosmap"]))
  return(degs)
}

input_female_pbmc <- "../Upload/Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"
input_male_pbmc <- "../Upload/Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"

cell_types <-c("Astrocyte 2","Microglia","Oligodendrocyte 1","OPC", "Inhibitory neuron 1", "Excitatory neuron 4")

cellt_male <- c("CD56-Dim, CD16 NK cell", "Naive CD4+ T cell", "Naive CD8+ T cell", "Transitional B cell") #unique(degs_pbmc_male$cluster_id) #c("Naive CD4+ T cells","Naive CD8+ T cell")
cellt_female <- cellt_male

celltype_to_order <- c("CD4+ T cells")

degs_pbmc_male <- read_pbmc_input(input_male_pbmc)
degs_pbmc_female <- read_pbmc_input(input_female_pbmc)
degs_pbmc_male$region <- "PBMC"
degs_pbmc_female$region <- "PBMC"
degs_pbmc_male$sex <- "M"
degs_pbmc_female$sex <- "F"
celltype_to_order <- "Naive CD8+ T cell"


cell_types <- unique(degs_brain$cluster_id) #c("Astrocyte 2","Microglia","Oligodendrocyte 1","OPC", "Inhibitory neuron 1", "Excitatory neuron 4")

cellt_male <- unique(degs_pbmc_male$cluster_id)  # c("CD56-Dim, CD16 NK cell", "Naive CD4+ T cell", "Naive CD8+ T cell", "Transitional B cell") 
cellt_female <- cellt_male

## Figure b
degs_zebra_male <- degs_zebra[degs_zebra$sex == "M",]
degs_rosmap_male <- degs_rosmap[degs_rosmap$sex == "M",]
#degs_brain_male <-  degs_brain_male[, colnames(degs_brain_male) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
#degs_brain_male <- degs_brain_male[degs_brain_male$cluster_id %in% cell_types,]

#degs_pbmc_male <- read_pbmc_input(input_male_pbmc)
#degs_pbmc_male <- degs_pbmc_male[degs_pbmc_male$cluster_id %in% cellt_male,]

male_genes <- degs_pbmc_male$gene[degs_pbmc_male$p_adj.loc < 0.05 & 
                                    #degs_pbmc_male$cluster_id %in% cellt_male & 
                                    degs_pbmc_male$gene %in% degs_zebra_male$gene &
                                    degs_pbmc_male$gene %in% degs_rosmap_male$gene] #read.csv("Alzheimer_male", header = F)$V1

plot_male_degs <- plotComparison(degs_zebra_male[degs_zebra_male$gene %in% male_genes,],
                                 degs_rosmap_male[degs_rosmap_male$gene %in% male_genes,],
                                 degs_pbmc_male[degs_pbmc_male$gene %in% male_genes,], cellt_male, cellt_male)
# plot_male_degs$logFC[plot_male_degs$logFC > 0.5] <- 0.5
# plot_male_degs$logFC[plot_male_degs$logFC < (- 0.5)] <- (-0.5)

plot_male <- ggplot() + 
  geom_tile(data = plot_male_degs, aes(y= cluster_id, x=gene, fill = logFC)) + 
  facet_grid(tissue~., scales = "free", space = "free")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-0.5,0.5))+ 
  theme_adrc()+ theme(axis.title.x=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())+ ylab("Celltypes")+ xlab("Genes")


cell_types <- unique(degs_brain$cluster_id) #c("Astrocyte 2","Microglia","Oligodendrocyte 1","OPC", "Inhibitory neuron 1", "Excitatory neuron 4")

cellt_male <- unique(degs_pbmc_male$cluster_id)  # c("CD56-Dim, CD16 NK cell", "Naive CD4+ T cell", "Naive CD8+ T cell", "Transitional B cell") 
cellt_female <- cellt_male

## Figure b
degs_zebra_female <- degs_zebra[degs_zebra$sex == "F",]
degs_rosmap_female <- degs_rosmap[degs_rosmap$sex == "F",]

female_genes <- degs_pbmc_female$gene[degs_pbmc_female$p_adj.loc < 0.05 & 
                                        #degs_pbmc_male$cluster_id %in% cellt_male & 
                                        degs_pbmc_female$gene %in% degs_zebra_female$gene &
                                        degs_pbmc_female$gene %in% degs_rosmap_female$gene] #read.csv("Alzheimer_male", header = F)$V1

plot_female_degs <- plotComparison(degs_zebra_female[degs_zebra_female$gene %in% female_genes,],
                                   degs_rosmap_female[degs_rosmap_female$gene %in% female_genes,],
                                   degs_pbmc_female[degs_pbmc_female$gene %in% female_genes,], cellt_female, cellt_female)
# plot_male_degs$logFC[plot_male_degs$logFC > 0.5] <- 0.5
# plot_male_degs$logFC[plot_male_degs$logFC < (- 0.5)] <- (-0.5)


plot_female <- ggplot() + 
  geom_tile(data = plot_female_degs, aes(y= cluster_id, x=gene, fill = logFC)) + 
  facet_grid(tissue~., scales = "free", space = "free")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-0.5,0.5))+ 
  theme_adrc()+ theme(axis.title.x=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())+ ylab("Celltypes")+ xlab("Genes")


plot_male  <- plot_male + theme(legend.position = "none")
plot_female  <- plot_female + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

plot <- plot_male+ plot_female

save_plot("figures/Supplement_8_a.svg", plot, base_width = 150, base_height = 50, unit = "mm")