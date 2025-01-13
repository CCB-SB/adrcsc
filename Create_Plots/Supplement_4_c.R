library(ggplot2)
library(data.table)
library(stringr)
library(cowplot)
library(aplot)
source("scripts/helper.R")

data_merged <-
  read.csv("scCODA/broadAnnotation/Output_All.csv")

data_merged <-
  data_merged[grep("Diagnosis", data_merged$Covariate),]
data_merged$Cell.Type <-
  factor(data_merged$Cell.Type, levels = CellCluster_order)


data_merged$significant <- T
data_merged$significant[data_merged$Cell.Type == "Dendritic cell"] <-
  F
data_merged$significant[abs(data_merged$log2.fold.change) < 0.5] <-
  F

data_merged$Disease <-
  str_replace(data_merged$Covariate, ".*\\[T\\.", "")
map <- c("MCI", "PD", "PD-MCI", "AD")
names(map) <-
  c(
    "Mild Cognitive Impairment]",
    "Parkinson's Disease only]",
    "Parkinson's Disease with MCI]",
    "Z.Disase.AD]"
  )
data_merged$Disease_short <- map[data_merged$Disease]
data_merged$Disease_short <-
  factor(data_merged$Disease_short,
         levels = c("AD",  "MCI", "PD", "PD-MCI"))

plot <-
  ggplot(data_merged,
         aes(x = Cell.Type, y = Disease_short, fill = log2.fold.change)) +
  geom_tile() +
  geom_point(data = data_merged[data_merged$significant, ], size = .5) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "log2 \nFold-change"
  ) +
  theme_adrc() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) + xlab("") + ylab("")

labels = ggplot(data_merged, aes(x = Cell.Type, y = 1, fill = Cell.Type)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill = "none")

plot <- (plot %>% insert_top(labels, height = .25))

ggsave(
  "figures/scCODA_Broad_All.svg",
  plot,
  width = 8,
  height = 5,
  unit = "cm"
)




labels = ggplot(data_merged, aes(x = Cell.Type, y = 1, fill = Cell.Type)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() #
legend <- get_legend(labels)
ggsave(
  "figures/BroadCelltypes_Legend.svg",
  legend,
  width = 4,
  height = 9,
  unit = "cm"
)
