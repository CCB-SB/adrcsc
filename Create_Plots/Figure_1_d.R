library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(ggforce)

source("scripts/helper.R")

adrc <- readRDS("data/CompleteObjectAnnotated_final.rds")

meta.data <- adrc@meta.data
meta.data$L2[meta.data$L2 == "CD4+ Proliferating"] <-
  "CD4+ Proliferating T cell"
meta.data$L1[meta.data$L1 == "CD4+ Proliferating"] <-
  "CD4+ Proliferating T cell"

l1 <- data.frame(table(meta.data$Manual_Annotation))
l2 <- data.frame(table(meta.data$L1))
l3 <- data.frame(table(meta.data$L2))

l1$level <- 3
l2$level <- 2
l3$level <- 1


list_factor <- c()
l3_order <- l3$Var1[order(l3$Freq, decreasing = F)]

for (ct in l3_order) {
  subs <- unique(meta.data$L1[meta.data$L2 == ct])
  order <- l2$Var1[order(l2$Freq, decreasing = F)]
  l2_order <- order[order %in% subs]
  
  list_factor <- c(list_factor, as.character(ct))
  for (ct2 in l2_order) {
    subs <- unique(meta.data$Manual_Annotation[meta.data$L1 == ct2])
    order <- l1$Var1[order(l1$Freq, decreasing = F)]
    list_factor <- c(list_factor, as.character(ct2))
    list_factor <-
      c(list_factor, as.character(order[order %in% subs]))
  }
}


data_to_plot <- rbind(l1, l3)
data_to_plot$Freq <- data_to_plot$Freq / sum(l1$Freq) * 100

data_to_plot$Var1 <-
  factor(data_to_plot$Var1, levels = unique(list_factor))

test.matrix.wide <- tidyr::spread(data_to_plot, level, Freq) %>%
  tidyr::replace_na(list(`1` = 0, `2` = 0)) %>%
  mutate(y = 100 -  cumsum(`1`), yend = 100 -  cumsum(`1`)) %>%
  add_row(y = 0, yend = 0)

data_to_plot$label <- as.character(data_to_plot$Var1)
data_to_plot$label[data_to_plot$level > 1 &
                     data_to_plot$Freq < 1.5] <- ""

plot <-
  ggplot(data_to_plot, aes(x = level, fill = Var1)) + geom_bar(aes(y = Freq),
                                                               width = 0.2,
                                                               colour = "black",
                                                               stat = "identity") +
  geom_text(
    aes(
      label = gsub("cell", "", gsub("Cells", "",  data_to_plot$label)),
      y = data_to_plot$Freq,
      hjust = 0,
      x = ifelse(data_to_plot$level == 1, 1.2, 3.2)
    ),
    position = position_stack(vjust = 0.5),
    stat = "identity",
    size = 3
  ) +
  scale_fill_manual(values = color_v) +
  geom_segment(data = test.matrix.wide,
               colour = "black",
               aes(
                 x = 1 + 0.2 / 2,
                 xend = 3 - 0.2 / 2,
                 y = y,
                 yend = yend
               )) +
  facet_zoom(
    ylim = c(92.65, 100),
    horizontal = T,
    zoom.size = .7
  ) +
  scale_x_continuous(
    limits = c(0.8, 4.5),
    breaks = c(1, 3),
    labels = c("Broad Annotation", "Fine Annotation")
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") + ylab("Proportion of cells in the dataset (%)")

ggsave(
  "figures/figure_1_d.svg",
  plot,
  width = 10.5,
  height = 7,
  unit = "cm"
)

write.csv(data_to_plot, "SourceData/Figure1d.csv")
