library(data.table)
library(stringr)
library(ggplot2)
library(cowplot)
library(aplot)
source("scripts/helper.R")


data_fem <- read.csv("scCODA/broadAnnotation/Output_Female.csv")
data_fem$Sex <- "female"
data_male <- read.csv("scCODA/broadAnnotation/Output_Male.csv")
data_male$Sex <- "male"

data_merged <- rbind(data_fem, data_male)

data_merged <-
  data_merged[grep("Diagnosis", data_merged$Covariate),]
data_merged$Cell.Type <-
  factor(data_merged$Cell.Type, levels = CellCluster_order)

data_merged$significant <- T
data_merged$significant[data_merged$Cell.Type == "CD4+ T cell"] <- F
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
plot1 <-
  ggplot(data_merged[data_merged$Sex == "female" , ],
         aes(x = Cell.Type, y = Disease_short, fill = log2.fold.change)) +
  geom_tile() +
  scale_x_discrete(label = labels_cell) +
  geom_point(data = data_merged[data_merged$significant &
                                  data_merged$Sex == "female", ], size = .5) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "log2\nFold-change"
  ) +
  facet_grid(Sex ~ .) +
  theme_adrc() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) + xlab("") + ylab("")

plot2 <-
  ggplot(data_merged[data_merged$Sex == "male", ],
         aes(x = Cell.Type, y = Disease_short, fill = log2.fold.change)) +
  geom_tile() +
  scale_x_discrete(label = labels_cell) +
  geom_point(data = data_merged[data_merged$significant &
                                  data_merged$Sex == "male", ], size = .5) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "log2\nFold-change"
  ) +
  facet_grid(Sex ~ .) +
  theme_adrc() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) + xlab("") + ylab("")


labels = ggplot(data_merged[data_merged$Sex == "female", ], aes(x = Cell.Type, y =
                                                                  1, fill = Cell.Type)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill = "none")



plot <-
  (plot1 %>% insert_top(labels, height = .25)) %>% insert_bottom(plot2, height = 1.0)

ggsave(
  "figures/scCODA_Broad.svg",
  plot,
  width = 7,
  height = 6,
  unit = "cm"
)


data_fem <- read.csv("scCODA/fineAnnotation/Output_Female.csv")
data_fem$Sex <- "female"
data_male <- read.csv("scCODA/fineAnnotation/Output_Male.csv")
data_male$Sex <- "male"

data_merged <- rbind(data_fem, data_male)

data_merged$Cell.Type[data_merged$Cell.Type == "CD4+ Proliferating T cell"] <-
  "CD4+ Proliferating"
data_merged <-
  data_merged[!data_merged$Cell.Type %in% c("HSPC / Progenitor TO_DELETE", "Red Blood Cells"), ]

data_merged$CellAbbr <-
  labels_cell[as.character(data_merged$Cell.Type)]
data_merged$CellAbbr <-
  factor(data_merged$CellAbbr, levels = labels_cell[celltype_cluster_order])

data_merged <-
  data_merged[grep("Diagnosis", data_merged$Covariate),]
data_merged$Cell.Type <-
  factor(data_merged$Cell.Type, levels = celltype_cluster_order)

data_merged$significant <- T
data_merged$significant[data_merged$Cell.Type == "Treg CD4+ cell"] <-
  F
data_merged$significant[abs(data_merged$log2.fold.change) < 1] <- F
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
plot1 <-
  ggplot(data_merged[data_merged$Sex == "female" , ],
         aes(x = Cell.Type, y = Disease_short, fill = log2.fold.change)) +
  geom_tile() +
  scale_x_discrete(label = labels_cell) +
  geom_point(data = data_merged[data_merged$significant &
                                  data_merged$Sex == "female", ], size = .5) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "log2\nFold-change"
  ) +
  facet_grid(Sex ~ .) +
  theme_adrc() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) + xlab("") + ylab("")

plot2 <-
  ggplot(data_merged[data_merged$Sex == "male", ],
         aes(x = Cell.Type, y = Disease_short, fill = log2.fold.change)) +
  geom_tile() +
  scale_x_discrete(label = labels_cell) +
  geom_point(data = data_merged[data_merged$significant &
                                  data_merged$Sex == "male", ], size = .5) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "log2\nFold-change"
  ) +
  facet_grid(Sex ~ .) +
  theme_adrc() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) + xlab("") + ylab("")


labels = ggplot(data_merged[data_merged$Sex == "female", ], 
                aes(x = Cell.Type, y = 1, fill = Cell.Type)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill = "none")



plot <-
  (plot1 %>% insert_top(labels, height = .25)) %>% insert_bottom(plot2, height = 1.0)

ggsave(
  "figures/scCODA_Fine.svg",
  plot,
  width = 12,
  height = 7,
  unit = "cm"
)
