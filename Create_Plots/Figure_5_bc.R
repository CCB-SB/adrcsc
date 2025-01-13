library(ggplot2)
library(cowplot)
library(aplot)
library(ggpubr)
library(data.table)

source("scripts/helper.R")

read_pbmc_input <- function(path) {
  degs_pbmc <- read.csv(path, sep = "\t")
  
  degs_pbmc <- degs_pbmc[degs_pbmc$contrast == "ADvsHC" , ]
  degs_pbmc$tissue <- "pbmc"
  degs_pbmc <-
    degs_pbmc[, colnames(degs_pbmc) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
  return(degs_pbmc)
}

plotComparison <- function(brain_data, pbmc_data) {
  degs <- rbind(brain_data, pbmc_data)
  degs <- degs[degs$gene %in% degs$gene[degs$tissue == "brain"], ]
  degs <- degs[degs$gene %in% degs$gene[degs$tissue == "pbmc"], ]
  
  freq_cellt <-
    as.data.frame(table(degs$cluster_id[degs$tissue == "pbmc"]))
  celltype_to_order <-
    freq_cellt$Var1[freq_cellt$Freq == max(freq_cellt$Freq)]
  
  gene_list_to_order <- unique(degs$gene)
  fc_to_order <- lapply(gene_list_to_order, function(gene) {
    if (gene %in% degs$genes[degs$cluster_id %in% celltype_to_order]) {
      degs$logFC[degs$cluster_id %in% celltype_to_order &
                   degs$gene == gene]
    } else
    {
      degs$logFC[degs$gene == gene][1]
    }
  })
  levels <-
    unique(degs$gene[degs$cluster_id %in% celltype_to_order][order(degs$logFC[degs$cluster_id %in% celltype_to_order], decreasing = T)])
  degs$gene <- factor(degs$gene, levels = levels)
  
  degs$tissue <- factor(degs$tissue, levels = c("pbmc", "brain"))
  
  degs$logFC[degs$logFC > 0.5] <- 0.5
  degs$logFC[degs$logFC < (-0.5)] <- (-0.5)
  degs$p_adj.loc[degs$p_adj.loc == 0] <-
    min(degs$p_adj.loc[degs$p_adj.loc != 0])
  
  degs$p_adj.loc_norm <- NA
  degs$p_adj.loc_norm[degs$tissue == "pbmc"] <-
    -log(degs$p_adj.loc[degs$tissue == "pbmc"]) / max(-log(degs$p_adj.loc[degs$tissue == "pbmc"]))
  degs$p_adj.loc_norm[degs$tissue == "brain"] <-
    -log(degs$p_adj.loc[degs$tissue == "brain"]) / max(-log(degs$p_adj.loc[degs$tissue == "brain"]))
  
  return(degs)
}


path <- paste("BrainData/human_cortex_sex_condition_marker.csv")
degs_brain <- read.csv(path, sep = ",")
degs_brain <- degs_brain[degs_brain$contrast == "AD-CT", ]
degs_brain$cluster_id <- degs_brain$cell_type

degs_brain$tissue <- "brain"
degs_brain$logFC <- degs_brain$logFC
degs_brain$p_adj.loc <- degs_brain$FDR
degs_brain$gene <- degs_brain$gene

degs_brain <-
  degs_brain[, colnames(degs_brain) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id", "sex")]


input_female_pbmc <-
  "results_female/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"
input_male_pbmc <-
  "results_male/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"

cell_types <-
  c(
    "Astrocyte 2",
    "Microglia",
    "Oligodendrocyte 1",
    "OPC",
    "Inhibitory neuron 1",
    "Excitatory neuron 4"
  )


################################################################################


degs_pbmc_female <- read_pbmc_input(input_female_pbmc)
degs_brain_female <- degs_brain[degs_brain$sex == "F", ]
degs_brain_female <-
  degs_brain_female[, colnames(degs_brain_female) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
degs_brain_female <-
  degs_brain_female[degs_brain_female$cluster_id %in% cell_types, ]

unique(degs_pbmc_female$cluster_id[degs_pbmc_female$p_adj.loc < 0.05])
cellt_female <- unique(degs_pbmc_female$cluster_id)
female_genes <-
  unique(degs_pbmc_female$gene[degs_pbmc_female$p_adj.loc < 0.05  &
                                 abs(degs_pbmc_female$logFC) > 0.6])
female_genes <-
  degs_brain_female$gene[degs_brain_female$gene %in% female_genes &
                           !is.na(degs_brain_female$p_adj.loc) &
                           degs_brain_female$p_adj.loc < 0.05]

use_clusters <-
  unique(degs_pbmc_female$cluster_id[degs_pbmc_female$p_adj.loc < 0.05 &
                                       degs_pbmc_female$gene %in% female_genes])

plot_female_degs <-
  plotComparison(degs_brain_female[degs_brain_female$gene %in% female_genes, ], degs_pbmc_female[degs_pbmc_female$gene %in% female_genes &
                                                                                                   degs_pbmc_female$cluster_id %in% use_clusters, ])

plot_female <- ggplot() +
  geom_point(
    data = plot_female_degs,
    aes(
      y = cluster_id,
      x = gene,
      fill = logFC,
      size = p_adj.loc_norm,
      colour = p_adj.loc < 0.05
    ),
    shape = 21
  ) +
  facet_grid(tissue ~ ., scales = "free", space = "free") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    limits = c(-0.5, 0.5)
  ) +
  scale_colour_manual(values = c("white", "black"), guide = "none") +
  theme_adrc() + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  scale_size(range = c(0, 3)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot_female
ggsave(
  "figures/figure_5_c.svg",
  plot_female,
  width = 10,
  height = 5,
  unit = "cm"
)

################################################################################

degs_pbmc_male <- read_pbmc_input(input_male_pbmc)
degs_brain_male <- degs_brain[degs_brain$sex == "M", ]
degs_brain_male <-
  degs_brain_male[, colnames(degs_brain_male) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
degs_brain_male <-
  degs_brain_male[degs_brain_male$cluster_id %in% cell_types, ]

male_genes <-
  unique(degs_pbmc_male$gene[degs_pbmc_male$p_adj.loc < 0.05 &
                               abs(degs_pbmc_male$logFC) > 0.6])
male_genes <-
  degs_brain_male$gene[degs_brain_male$gene %in% male_genes &
                         !is.na(degs_brain_male$p_adj.loc) &
                         degs_brain_male$p_adj.loc < 0.05]
use_clusters <-
  unique(degs_pbmc_male$cluster_id[degs_pbmc_male$p_adj.loc < 0.05 &
                                     degs_pbmc_male$gene %in% male_genes])
use_clusters <- use_clusters[use_clusters != "Red Blood Cells"]
male_genes <-
  male_genes[male_genes %in% degs_pbmc_male$gene[degs_pbmc_male$p_adj.loc < 0.05 &
                                                   degs_pbmc_male$cluster_id %in% use_clusters]]
cellt_male <- unique(degs_pbmc_male$cluster_id)

plot_male_degs <-
  plotComparison(degs_brain_male[degs_brain_male$gene %in% male_genes, ], degs_pbmc_male[degs_pbmc_male$gene %in% male_genes &
                                                                                           degs_pbmc_male$cluster_id %in% use_clusters, ])

plot_male <- ggplot() +
  geom_point(
    data = plot_male_degs,
    aes(
      y = cluster_id,
      x = gene,
      fill = logFC,
      size = p_adj.loc_norm,
      colour = p_adj.loc < 0.05
    ),
    shape = 21
  ) +
  facet_grid(tissue ~ ., scales = "free", space = "free") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    limits = c(-0.5, 0.5)
  ) +
  scale_colour_manual(values = c("white", "black"), guide = "none") +
  theme_adrc() + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) + scale_size(range = c(0, 3)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


ggsave(
  "figures/figure_5_b.svg",
  plot_male,
  width = 12.5,
  height = 5,
  unit = "cm"
)