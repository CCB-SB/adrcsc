library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(aplot)
library(ComplexHeatmap)
library(viridis)
library(ggplotify)
library(ggh4x)

source("scripts/helper.R")


color_v2 <- color_v
names(color_v2) <- labels_cell[names(color_v2)]
color_v <- c(color_v, color_v2)

degs_male <-
  read.csv(
    "results_male/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )
degs_male <-
  degs_male[degs_male$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"), ]
degs_male$Sex = "male"
degs_male <-
  degs_male[degs_male$p_adj.loc < 0.05 & abs(degs_male$logFC) > 0.5, ]
degs_female <-
  read.csv(
    "results_female/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv",
    sep = "\t"
  )
degs_female$Sex = "female"
degs_female <-
  degs_female[degs_female$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"), ]
degs_female <-
  degs_female[degs_female$p_adj.loc < 0.05 &
                abs(degs_female$logFC) > 0.5, ]

# Figure 3 a)


degs_number <- rbind(degs_male, degs_female)

c <-
  seq(nrow(degs_number[, colnames(degs_number) %in% c("cluster_id", "contrast", "Sex")]))
counts <-
  aggregate(c ~ ., data = degs_number[, colnames(degs_number) %in% c("cluster_id", "contrast", "Sex")], FUN =
              length,drop=FALSE)

counts$cluster_id[counts$cluster_id == "CD4+ Proliferating T cell"] <-
  "CD4+ Proliferating"
counts <-
  counts[!counts$cluster_id %in% c("HSPC / Progenitor TO_DELETE", "Red Blood Cells"), ]

counts$contrast <-
  factor(counts$contrast, levels = c("PDvsHC", "MCIvsHC", "ADvsHC"))

counts$CellAbbr <- labels_cell[counts$cluster_id]
counts$CellAbbr <-
  factor(counts$CellAbbr, levels = labels_cell[celltype_cluster_order])
counts$cluster_id <-
  factor(counts$cluster_id, levels = celltype_cluster_order)

labels = ggplot(counts, aes(x = CellAbbr, y = 1, fill = CellAbbr)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill = "none")

counts$c[counts$c == 0] <- NA

plot <- ggplot(counts[counts$contrast == "ADvsHC",], aes(x = CellAbbr, y = Sex))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "white", high = "#ff7f00", trans = "log", limits = c(1,1000), breaks = c(1,5,50,500),na.value="white") +
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(contrast~.)+
  scale_x_discrete(name = "", labels = labels_cell[as.character(counts$cluster_id)]) + scale_y_discrete(name = "") + 
  labs(fill='#DEGs ') + theme(legend.position="bottom") + theme(legend.title=element_text(size=8), legend.text=element_text(size=6))

plot2 <- ggplot(counts[counts$contrast == "MCIvsHC",], aes(x = CellAbbr, y = Sex))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "white", high = "#e31a1c", trans = "log", limits = c(1,1000), breaks = c(1,5,50,500),na.value="white") +
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(contrast~.)+
  scale_x_discrete(name = "", labels = labels_cell[as.character(counts$cluster_id)]) + scale_y_discrete(name = "") + 
  labs(fill='#DEGs ')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(legend.position="bottom") + theme(legend.title=element_text(size=8), legend.text=element_text(size=6))

plot3 <- ggplot(counts[counts$contrast == "PDvsHC",], aes(x = CellAbbr, y = Sex))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "white", high = "#6a3d9a", trans = "log", limits = c(1,1000), breaks = c(1,5,50,500),na.value="white") +
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(contrast~.)+
  scale_x_discrete(name = "", labels = labels_cell[as.character(counts$cluster_id)]) + scale_y_discrete(name = "") + labs(fill='#DEGs ') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(legend.position="bottom")+ theme(legend.title=element_text(size=8), legend.text=element_text(size=6))


legend1 <- cowplot::get_legend(plot)
legend2 <- cowplot::get_legend(plot2)
legend3 <- cowplot::get_legend(plot3)
plot <- plot + theme(legend.position="none")
plot2 <- plot2 + theme(legend.position="none")
plot3 <- plot3 + theme(legend.position="none")


p <- (plot3 %>% insert_top(labels, height=.25)) %>% insert_bottom(plot2, height = 1.0) %>% insert_bottom(plot, height = 1.0)


l <- as.ggplot(plot_grid(legend1, legend2,legend3, nrow = 1))
l <- l + theme(plot.margin = unit(c(0, 0, 1, 2), "cm"))
figure_a <- plot_grid(as.ggplot(p), l, nrow = 2, rel_heights = c(1,0.1))

ggsave(
  "figures/figure_3_a.svg",
  figure_a,
  width = 9,
  height = 6.5,
  unit = "cm"
)

write.csv(counts, "SourceData/Figure_3a.csv")
