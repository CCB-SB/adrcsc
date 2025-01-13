library(data.table)
library(cowplot)
library(ggplot2)
library(dplyr)

source("scripts/helper.R")



meta_data <- read.csv("data/metadata.csv")
pbmc <- readRDS("data/CompleteObjectAnnotated_final.rds")
pbmc <-
  subset(pbmc, subset = Sample.ID..From.Stanford. %in% meta_data$Sample.ID..From.Stanford.[meta_data$Visit == 1])

cell_counts_all <-
  pbmc@meta.data %>% group_by(L2, SCMD) %>%  summarise(n = n())

meta <-
  as.data.frame(pbmc@meta.data %>% group_by(SCMD, Sex, Age, Diagnosis, ApoE, Visit))
meta <- meta[, colnames(meta) %in% colnames(meta_data)]
meta <- meta[!duplicated(meta), ]

meta_info <- colnames(meta_data)[17:101]

for (ct in unique(cell_counts_all$L2)) {
  cell_counts <- cell_counts_all[cell_counts_all$L2 == ct, ]
  cell_counts <-
    cell_counts[colnames(cell_counts) %in% c("SCMD", "n")]
  print(cell_counts)
  meta <- merge(meta, cell_counts, by = c("SCMD"), all.x = T)
  colnames(meta)[length(colnames(meta))] <- ct
  meta[[ct]][is.na(meta[[ct]])] <- 0
}

corr_data <- NA
for (d in c(
  "Mild Cognitive Impairment",
  "Alzheimer's disease",
  "Healthy Control",
  "Parkinson's Disease only"
)) {
  for (ct in unique(cell_counts_all$L2)) {
    for (meta_col in unique(meta_info)) {
      if (!is.numeric(meta[[meta_col]][1]))
        next
      corr <-
        cor.test(meta[[ct]][meta$Diagnosis == d], meta[[meta_col]][meta$Diagnosis == d])
      m_info <-  meta[[meta_col]][meta$Diagnosis == d]
      data <-
        data.frame(
          cor = corr$estimate,
          p = corr$p.value,
          Diagnosis = d,
          celltype = ct ,
          meta = meta_col,
          l = length(m_info[!is.na(m_info)])
        )
      if (is.na(corr_data)[1]) {
        corr_data <- data
      } else {
        corr_data <- rbind(corr_data, data)
      }
    }
  }
}
corr_data$p.adj <- p.adjust(corr_data$p)

corr_data$celltype <-
  factor(corr_data$celltype, levels = CellCluster_order)

map <- c("MCI", "PD", "PD-MCI", "AD", "HC")
names(map) <-
  c(
    "Mild Cognitive Impairment",
    "Parkinson's Disease only",
    "Parkinson's Disease with MCI",
    "Alzheimer's disease" ,
    "Healthy Control"
  )
corr_data$Disease_short <- map[corr_data$Diagnosis]
corr_data$Disease_short <-
  factor(corr_data$Disease_short,
         levels = c("AD",  "MCI", "PD", "PD-MCI", "HC"))


plot <-
  ggplot(corr_data[corr_data$meta %in% corr_data$meta[corr_data$p.adj < 0.05], ], aes(x = celltype, y = meta, fill = cor)) + geom_tile() +
  geom_point(data = corr_data[corr_data$p.adj < 0.05, ]) +
  xlab("") + facet_grid(. ~ Diagnosis) +
  ylab("") + theme_adrc() +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "Pearson's correlation"
  ) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))


plot <- fill_title(plot, color_v)
ggsave(
  "figures/Supplement_11_c.svg",
  plot,
  width = 18,
  height = 5,
  unit = "cm"
)

tmp <-
  corr_data[corr_data$meta %in% corr_data$meta[corr_data$p.adj < 0.05], ]
write.csv(tmp, "SourceData/Supplement_11_c.csv")
