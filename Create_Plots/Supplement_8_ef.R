library(ggplot2)
library(cowplot)
library(aplot)
library(ggpubr)
library(viridis)
library(data.table)
source("scripts/helper.R")

input_male_pbmc <-
  "results_male/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"

bulk_male <- read.csv("ROSMAP_bulk/DEGs_male.csv")
path <- paste("BrainData/human_cortex_sex_condition_marker.csv")
degs_brain <- read.csv(path, sep = ",")
degs_brain <- degs_brain[degs_brain$contrast == "AD-CT",]
degs_brain$cluster_id <- degs_brain$cell_type

degs_brain$logFC <- degs_brain$logFC
degs_brain$p_adj.loc <- degs_brain$FDR
degs_brain$gene <- degs_brain$gene
degs_brain$region <- "Cortex"
degs_brain$tissue <- "Zebra"


degs_pbmc_male <- read.csv(input_male_pbmc, sep = "\t")
degs_pbmc_male <-
  degs_pbmc_male[degs_pbmc_male$contrast == "ADvsHC" ,]
degs_pbmc_male$tissue <- "pbmc"

a <- c(degs_pbmc_male$gene_id)
b <- c(bulk_male$X)

intersect_genelist <- unique(a[a %in% b])
intersect_genelist_gene <-
  unlist(lapply(intersect_genelist, function(gene)
    degs_pbmc_male$gene[degs_pbmc_male$gene_id == gene][1]))

df_to_plot <- data.frame(gene = intersect_genelist,
                         bulk = unlist(lapply(intersect_genelist, function(g)
                           ifelse(g %in% bulk_male$gene, bulk_male$logFC[bulk_male$X == g], NA))),
                         pbmc = unlist(lapply(intersect_genelist, function(g)
                           ifelse(g %in% degs_pbmc_male$gene_id[degs_pbmc_male$cluster_id == "CD56-Dim, CD16 NK cell"], degs_pbmc_male$logFC[degs_pbmc_male$cluster_id == "CD56-Dim, CD16 NK cell" &
                                                                                                                                               degs_pbmc_male$gene_id == g], NA))))

plot1 <-
  ggscatter(
    df_to_plot,
    x = "bulk",
    y = "pbmc",
    conf.int = T,
    add = "reg.line"
  ) +
  stat_cor(method = "pearson",
           label.x = 0,
           label.y = 1.5) + scale_y_continuous(limits = c(-2, 2))

degs_brain_male <- degs_brain[degs_brain$sex == "M",]
degs_brain_male <-
  degs_brain_male[, colnames(degs_brain_male) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]

df_to_plot2 <- data.frame(gene = intersect_genelist,
                          bulk = unlist(lapply(intersect_genelist, function(g)
                            ifelse(g %in% bulk_male$gene, bulk_male$logFC[bulk_male$X == g], NA))),
                          brain = unlist(lapply(intersect_genelist_gene, function(g)
                            ifelse(g %in% degs_brain_male$gene, degs_brain_male$logFC[degs_brain_male$gene == g][1], NA))))

plot2 <-
  ggscatter(
    df_to_plot2,
    x = "bulk",
    y = "brain",
    conf.int = T,
    add = "reg.line"
  ) +
  stat_cor(method = "pearson",
           label.x = 0,
           label.y = 1.5) + scale_y_continuous(limits = c(-2, 2))

plot <- plot1 + plot2
ggsave("figures/Supplement_8f.svg",
       plot,
       width = 7,
       height = 3)
rownames(df_to_plot) <- df_to_plot$gene
rownames(df_to_plot2) <- df_to_plot2$gene
write.csv(merge(df_to_plot, df_to_plot2, by = c("gene", "bulk")),
          "SourceData/Supplement_8f.csv")
