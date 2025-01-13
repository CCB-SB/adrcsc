library(ggplot2)
library(cowplot)
library(aplot)
library(ggpubr)
library(viridis)
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

plotComparison <-
  function(brain_data,
           pbmc_data,
           rosmap_data,
           celltypes,
           celltype_to_order) {
    degs <- rbind(brain_data, pbmc_data, rosmap_data)
    
    degs <- degs[!is.na(degs$p_adj.loc), ]
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
    
    degs$tissue <-
      factor(degs$tissue, levels = c("pbmc", "Zebra", unique(degs$tissue[!degs$tissue %in% c("Zebra", "pbmc")])))
    
    degs$logFC[degs$logFC > 0.5] <- 0.5
    degs$logFC[degs$logFC < (-0.5)] <- (-0.5)
    degs$p_adj.loc[degs$p_adj.loc == 0] <-
      min(degs$p_adj.loc[degs$p_adj.loc != 0])
    
    degs$p_adj.loc_norm <- NA
    degs$p_adj.loc_norm[degs$tissue == "pbmc"] <-
      -log(degs$p_adj.loc[degs$tissue == "pbmc"]) / max(-log(degs$p_adj.loc[degs$tissue == "pbmc"]))
    degs$p_adj.loc_norm[degs$tissue == "Zebra"] <-
      -log(degs$p_adj.loc[degs$tissue == "Zebra"]) / max(-log(degs$p_adj.loc[degs$tissue == "Zebra"]))
    degs$p_adj.loc_norm[degs$tissue == "Rosmap"] <-
      -log(degs$p_adj.loc[degs$tissue == "Rosmap"]) / max(-log(degs$p_adj.loc[degs$tissue == "Rosmap"]))
    return(degs)
  }


path <- paste("BrainData/human_cortex_sex_condition_marker.csv")
degs_brain <- read.csv(path, sep = ",")
degs_brain <- degs_brain[degs_brain$contrast == "AD-CT", ]
degs_brain$cluster_id <- degs_brain$cell_type


degs_brain$logFC <- degs_brain$logFC
degs_brain$p_adj.loc <- degs_brain$FDR
degs_brain$gene <- degs_brain$gene
degs_brain$region <- "Cortex"
degs_brain$tissue <- "Zebra"
degs_zebra <-
  degs_brain[, colnames(degs_brain) %in% c("tissue",
                                           "logFC",
                                           "p_adj.loc",
                                           "gene",
                                           "cluster_id",
                                           "sex",
                                           "region")]


input_female_pbmc <-
  "../Upload/Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"
input_male_pbmc <-
  "../Upload/Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv"


degs_pbmc_male <- read_pbmc_input(input_male_pbmc)
degs_pbmc_female <- read_pbmc_input(input_female_pbmc)
degs_pbmc_male$region <- "PBMC"
degs_pbmc_female$region <- "PBMC"
degs_pbmc_male$sex <- "M"
degs_pbmc_female$sex <- "F"
celltype_to_order <- "Naive CD8+ T cell"


path <- paste("BrainData/Brain_Rosmap_Colltected.csv")
degs_brain <- read.csv(path, sep = ",")
degs_brain$sex <- degs_brain$Sex

degs_brain$tissue <- "Rosmap"
degs_rosmap <-
  degs_brain[, colnames(degs_brain) %in% c("tissue",
                                           "logFC",
                                           "p_adj.loc",
                                           "gene",
                                           "cluster_id",
                                           "sex",
                                           "region")]



degs_pbmc_male <-
  degs_pbmc_male[abs(degs_pbmc_male$logFC) > .5 &
                   degs_pbmc_male$p_adj.loc < 0.05, ]
degs_pbmc_female <-
  degs_pbmc_female[abs(degs_pbmc_female$logFC) > .5 &
                     degs_pbmc_female$p_adj.loc < 0.05, ]
degs_brain_zebra <-
  degs_zebra[abs(degs_zebra$logFC) > .5 &
               degs_zebra$p_adj.loc < 0.05, ]
degs_collected_rosmap <-
  degs_rosmap[abs(degs_rosmap$logFC) > .5 &
                degs_rosmap$p_adj.loc < 0.05, ]

male_genes <- intersect(
  degs_pbmc_male$gene,
  intersect(degs_brain_zebra$gene[degs_brain_zebra$sex == "M"],
            degs_collected_rosmap$gene[degs_collected_rosmap$sex == "M"])
)

female_genes <- intersect(
  degs_pbmc_female$gene,
  intersect(degs_brain_zebra$gene[degs_brain_zebra$sex == "F"],
            degs_collected_rosmap$gene[degs_collected_rosmap$sex == "F"])
)


gene_info <- read.csv("data/Gene_Information_Brain.csv")

gene_info[gene_info$gene %in% female_genes, ]
gene_info[gene_info$gene %in% male_genes, ]


degs_pbmc_male <-
  degs_pbmc_male[abs(degs_pbmc_male$logFC) > .5 &
                   degs_pbmc_male$p_adj.loc < 0.05 &
                   degs_pbmc_male$gene %in% male_genes, ]
degs_pbmc_female <-
  degs_pbmc_female[abs(degs_pbmc_female$logFC) > .5 &
                     degs_pbmc_female$p_adj.loc < 0.05 &
                     degs_pbmc_female$gene %in% female_genes, ]
degs_brain_zebra_male <-
  degs_zebra[abs(degs_zebra$logFC) > .5 &
               degs_zebra$p_adj.loc < 0.05 &
               degs_zebra$gene %in% male_genes & degs_zebra$sex == "M", ]
degs_brain_zebra_female <-
  degs_zebra[abs(degs_zebra$logFC) > .5 &
               degs_zebra$p_adj.loc < 0.05 &
               degs_zebra$gene %in% female_genes & degs_zebra$sex == "F", ]
degs_collected_rosmap_male <-
  degs_rosmap[abs(degs_rosmap$logFC) > .5 &
                degs_rosmap$p_adj.loc < 0.05 &
                degs_rosmap$gene %in% male_genes & degs_rosmap$sex == "M", ]
degs_collected_rosmap_female <-
  degs_rosmap[abs(degs_rosmap$logFC) > .5 &
                degs_rosmap$p_adj.loc < 0.05 &
                degs_rosmap$gene %in% female_genes & degs_rosmap$sex == "F", ]


degs_pbmc_male$sign <-
  sapply(degs_pbmc_male$gene , function(g)
    ifelse(any(degs_pbmc_male$logFC[degs_pbmc_male$gene == g] > 0), ifelse(any(
      degs_pbmc_male$logFC[degs_pbmc_male$gene == g] < 0
    ), "+/-", "+"), "-"))
degs_brain_zebra_male$sign <-
  sapply(degs_brain_zebra_male$gene , function(g)
    ifelse(
      any(degs_brain_zebra_male$logFC[degs_brain_zebra_male$gene == g] > 0),
      ifelse(any(degs_brain_zebra_male$logFC[degs_brain_zebra_male$gene == g] <
                   0), "+/-", "+"),
      "-"
    ))

degs_collected_rosmap_male$sign <-
  sapply(1:length(degs_collected_rosmap_male$gene) , function(i) {
    g <- degs_collected_rosmap_male$gene[i]
    r <- degs_collected_rosmap_male$region[i]
    ifelse(any(degs_collected_rosmap_male$logFC[degs_collected_rosmap_male$gene == g &
                                                  degs_collected_rosmap_male$region == r] > 0),
           ifelse(any(degs_collected_rosmap_male$logFC[degs_collected_rosmap_male$gene == g &
                                                         degs_collected_rosmap_male$region == r] < 0), "+/-", "+"),
           "-")
  })

data_to_plot <-
  rbind(degs_pbmc_male,
        degs_brain_zebra_male,
        degs_collected_rosmap_male)
data_to_plot$category <-
  sapply(data_to_plot$gene, function(g)
    gene_info$Category[gene_info$gene == g][1])
unique(data_to_plot$gene[is.na(data_to_plot$category)])


data_to_plot$category <-
  factor(
    data_to_plot$category,
    levels =  c(
      "AD-general",
      "Cognition",
      "Astrocytes",
      "Astrocytic Ferroptosis",
      "Microglia",
      "Ribosome",
      "Mitochondrium",
      "Plaques",
      "BBB and Immune",
      "Immune response",
      "Inflammation",
      "Not reported"
    )
  )

plot <-
  ggplot(data_to_plot, aes(x = region, y = gene)) + geom_point(aes(color = sign, shape = sign), size = 6) +
  facet_grid(category ~ tissue, scales = "free", space = "free") + theme_adrc() +
  theme(strip.text.y =  element_text(angle = 0)) +
  scale_color_manual(
    name = "Deregulation",
    values = c("#d73027", "gray", "#4575b4"),
    breaks = c("+", "+/-", "-"),
    labels = c("Up-regulated", "Up- and \ndown-regulated", "Downregulated")
  ) +
  scale_shape_manual(
    name = "Deregulation",
    values = c('+', "*", '-'),
    breaks = c("+", "+/-", "-"),
    labels = c("Up-regulated", "Up- and \ndown-regulated", "Downregulated")
  )
ggsave("figures/figure_5_d.svg",
       plot,
       width = 5.5,
       height = 5)




degs_pbmc_female$sign <-
  sapply(degs_pbmc_female$gene , function(g)
    ifelse(any(degs_pbmc_female$logFC[degs_pbmc_female$gene == g] > 0), ifelse(
      any(degs_pbmc_female$logFC[degs_pbmc_female$gene == g] < 0), "+/-", "+"
    ), "-"))
degs_brain_zebra_female$sign <-
  sapply(degs_brain_zebra_female$gene , function(g)
    ifelse(
      any(degs_brain_zebra_female$logFC[degs_brain_zebra_female$gene == g] > 0),
      ifelse(any(degs_brain_zebra_female$logFC[degs_brain_zebra_female$gene == g] <
                   0), "+/-", "+"),
      "-"
    ))
degs_collected_rosmap_female$sign <-
  sapply(1:length(degs_collected_rosmap_female$gene) , function(i) {
    g <- degs_collected_rosmap_female$gene[i]
    r <- degs_collected_rosmap_female$region[i]
    
    ifelse(any(degs_collected_rosmap_female$logFC[degs_collected_rosmap_female$gene == g &
                                                    degs_collected_rosmap_female$region == r] > 0),
           ifelse(any(
             degs_collected_rosmap_female$logFC[degs_collected_rosmap_female$gene == g &
                                                  degs_collected_rosmap_female$region == r] < 0
           ), "+/-", "+"),
           "-")
  })

data_to_plot <-
  rbind(degs_pbmc_female,
        degs_brain_zebra_female,
        degs_collected_rosmap_female)
data_to_plot$category <-
  sapply(data_to_plot$gene, function(g)
    gene_info$Category[gene_info$gene == g][1])
unique(data_to_plot$gene[is.na(data_to_plot$category)])


data_to_plot$category <-
  factor(
    data_to_plot$category,
    levels =  c(
      "AD-general",
      "Cognition",
      "Astrocytes",
      "Astrocytic Ferroptosis",
      "Microglia",
      "Microglia/Immune",
      "Ribosome",
      "Mitochondrium",
      "Plaques",
      "BBB and Immune",
      "Immune response",
      "Inflammation",
      "Not reported"
    )
  )

plot <-
  ggplot(data_to_plot, aes(x = region, y = gene)) + geom_point(aes(color = sign, shape = sign), size = 6) +
  facet_grid(category ~ tissue, scales = "free", space = "free") + theme_adrc() +
  theme(strip.text.y =  element_text(angle = 0)) +
  scale_color_manual(
    name = "Deregulation",
    values = c("#d73027", "gray", "#4575b4"),
    breaks = c("+", "+/-", "-"),
    labels = c("Up-regulated", "Up- and \ndown-regulated", "Downregulated")
  ) +
  scale_shape_manual(
    name = "Deregulation",
    values = c('+', "*", '-'),
    breaks = c("+", "+/-", "-"),
    labels = c("Up-regulated", "Up- and \ndown-regulated", "Downregulated")
  )
plot
ggsave("figures/figure_5_e.svg",
       plot,
       width = 5,
       height = 2)