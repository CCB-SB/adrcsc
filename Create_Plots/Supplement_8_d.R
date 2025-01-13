get_signaling_genes <- function(ds1, ds2){
  
  signal_collected <- NA
  for (pw in unique(CellChatDB.human$interaction$pathway_name)){
    ligand_data <- ds1[ds1$gene %in% CellChatDB.human$interaction$ligand[CellChatDB.human$interaction$pathway_name == pw],]
    receptor_data <- ds2[ds2$gene %in% CellChatDB.human$interaction$receptor[CellChatDB.human$interaction$pathway_name == pw],]
    if (length(ligand_data$gene)<1 | length(receptor_data$gene)<1) next
    #print(CellChatDB.human$interaction$evidence[CellChatDB.human$interaction$pathway_name == pw][1])
    #print(unique(ligand_data))
    #print(unique(receptor_data$cell_type))
    data_pw <- data.frame(role = c(rep("ligand", length(ligand_data$gene)), rep("receptor", length(receptor_data$gene))), tissue = c(rep("pbmc", length(ligand_data$gene)), rep("brain", length(receptor_data$gene))),gene = c(ligand_data$gene, receptor_data$gene), FC = c(ligand_data$logFC, receptor_data$logFC), p = c(ligand_data$p_adj.loc, receptor_data$FDR), cluster =  c(ligand_data$cluster_id, receptor_data$cell_type))
    data_pw$pathway <- pw
    if(is.na(signal_collected)[1]){signal_collected <- data_pw} else{signal_collected <- rbind(signal_collected, data_pw)}
    
    receptor_data <- ds1[ds1$gene %in% CellChatDB.human$interaction$receptor[CellChatDB.human$interaction$pathway_name == pw],]
    ligand_data <- ds2[ds2$gene %in% CellChatDB.human$interaction$ligand[CellChatDB.human$interaction$pathway_name == pw],]
    if (length(ligand_data$gene)<1 | length(receptor_data$gene)<1) next
    print(pw)
    #print(unique(ligand_data))
    #print(unique(receptor_data$cell_type))
    data_pw <- data.frame(role = c(rep("ligand", length(ligand_data$gene)), rep("receptor", length(receptor_data$gene))), tissue = c(rep("brain", length(ligand_data$gene)), rep("pbmc", length(receptor_data$gene))),gene = c(ligand_data$gene, receptor_data$gene), FC = c(ligand_data$logFC, receptor_data$logFC), p = c(ligand_data$FDR, receptor_data$p_adj.loc), cluster =  c(ligand_data$cell_type, receptor_data$cluster_id))
    data_pw$pathway <- pw
    if(is.na(signal_collected)[1]){signal_collected <- data_pw} else{signal_collected <- rbind(signal_collected, data_pw)}
  }
  return(signal_collected)
}

plot_signals <- function(signal_collected){
  
  ggplot(signal_collected, aes(x = cluster, y = gene), color = "black")+
    geom_point(data = signal_collected[signal_collected$FC> 0,], aes(x = cluster, y = gene, fill = FC, color = ifelse(p<0.05, "black", "white")), shape = 24)+
    geom_point(data = signal_collected[signal_collected$FC< 0,], aes(x = cluster, y = gene, fill = FC, color = ifelse(p<0.05, "black", "white")), shape = 25)+
    #geom_point(size = 5)+ 
    facet_grid(role ~ tissue, scales = "free", space= "free")+ 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+ 
    scale_color_manual(values = c("black", "white")) + theme_classic()+ xlab("")+ylab("")+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")+ guides(color = "none")+ 
    labs(fill = "log2\nFold-change")
  
}

brain <- read.csv("BrainData/human_cortex_sex_condition_marker.csv")
brain <- brain[brain$contrast == "AD-CT",]
brain <- brain[brain$FDR < 0.05,]
pbmc_f <- read.csv("../Upload/Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.p_filtered.csv", sep = "\t")
pbmc_f$sex <- "F"
pbmc_m <- read.csv("../Upload/Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.p_filtered.csv", sep = "\t")
pbmc_m$sex <- "M"
pbmc <- rbind(pbmc_f, pbmc_m)
pbmc <- pbmc[pbmc$contrast == "ADvsHC",]

load("../Upload/external_data/CellChatDB.human.rda")
#CellChatDB



signal_collected_f <- get_signaling_genes(pbmc[pbmc$sex == "F",], brain[brain$sex == "F",])
plot2 = plot_signals(signal_collected_f[signal_collected_f$pathway == "CCL",])+ theme(legend.position = "right")
ggsave("figures/Signaling_Female_CCL_filtered.svg", plot2, width = 9, height = 7, unit = "cm")

# 
# brain <- read.csv("BrainData/human_cortex_sex_condition_marker.csv")
# pbmc_f <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep = "\t")
# pbmc_f$sex <- "F"
# pbmc_m <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep = "\t")
# pbmc_m$sex <- "M"
# pbmc <- rbind(pbmc_f, pbmc_m)
# brain <- brain[brain$contrast == "AD-CT",]
# pbmc <- pbmc[pbmc$contrast == "ADvsHC",]
# 
# load("CellChatDB.human.rda")
# CellChatDB <- CellChatDB.human
# 
# 
# signal_collected <- get_signaling_genes(pbmc, brain)
# plot1 = plot_signals(signal_collected[signal_collected$pathway == "MHC-I",])
# ggsave("figures/Signaling_All_MCI-I.svg", plot1)
# plot2 = plot_signals(signal_collected[signal_collected$pathway == "CCL",])
# ggsave("figures/Signaling_All_CCl.svg", plot2)
# 
# signal_collected_m <- get_signaling_genes(pbmc[pbmc$sex == "M",], brain[brain$sex == "M",])
# plot1 = plot_signals(signal_collected_m[signal_collected_m$pathway == "CCL",])
# ggsave("figures/Signaling_Male_CCL.svg", plot1)
# plot2 = plot_signals(signal_collected_m[signal_collected_m$pathway == "CD45",])
# ggsave("figures/Signaling_Male_CD45.svg", plot2)
# 
# signal_collected_f <- get_signaling_genes(pbmc[pbmc$sex == "F",], brain[brain$sex == "F",])
# plot1 = plot_signals(signal_collected_f[signal_collected_f$pathway == "MHC-I",])
# ggsave("figures/Signaling_Female_MHC-I.svg", plot1)
# plot2 = plot_signals(signal_collected_f[signal_collected_f$pathway == "CCL",])
# ggsave("figures/Signaling_Female_CCL.svg", plot2)
