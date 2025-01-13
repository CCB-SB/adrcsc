library(circlize)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(ggrepel)
library(data.table)
library(cowplot)

source("scripts/helper.R")

pathways_female <-
  read.csv("results_female/pathway_analysis_results_all.csv")
pathways_female$sex <- "female"
pathways_male <-
  read.csv("results_male/pathway_analysis_results_All.csv")
pathways_male$sex <- "male"

pw_collected <- rbind(pathways_female, pathways_male)

pw_count <-
  pw_collected %>% group_by(X.Name, contrast, sex) %>% summarise(PW_count = n())

colors <- c("#5D7322", "#9C489A", "black")
names(colors) <- c("only male", "only female", "female and male")

pw_count_gr1 <- pw_count[pw_count$sex == "female", ]
pw_count_gr2 <- pw_count[pw_count$sex == "male", ]

data_to_plot <-
  merge(pw_count_gr1,
        pw_count_gr2,
        by = c("X.Name", "contrast"),
        all = T)
data_to_plot$PW_count.x[is.na(data_to_plot$PW_count.x)] <- 0
data_to_plot$PW_count.y[is.na(data_to_plot$PW_count.y)] <- 0
data_to_plot$in_group <- "female and male"
data_to_plot$in_group[(data_to_plot$PW_count.x == 0)] <- "only male"
data_to_plot$in_group[(data_to_plot$PW_count.y == 0)] <-
  "only female"
counts <-
  data_to_plot %>% group_by(PW_count.x, PW_count.y, in_group, contrast) %>% summarise(num_pathways = n())
counts$label <-
  sapply(1:length(counts$PW_count.x), function(i)
    paste(data_to_plot$X.Name[data_to_plot$PW_count.x == counts$PW_count.x[i] &
                                data_to_plot$PW_count.y == counts$PW_count.y[i]], collapse = '\n'))
counts$label[counts$PW_count.x < 2 & counts$PW_count.y < 7] <- ""
counts$label[is.na(counts$label)] <- ""
counts$num_pathways[counts$num_pathways > 20] <- 20
counts$label <- ""

plot <-
  ggplot(
    counts,
    aes(
      x = PW_count.x,
      y = PW_count.y,
      size = num_pathways,
      label = label,
      color = in_group
    )
  ) +
  facet_grid2(. ~ contrast, scales = "free", independent = "y") +
  scale_color_manual(values = colors) + labs(color = "Significant in", size = "#Occurences") +
  geom_point() +
  geom_text_repel(size = 1) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = 3) +
  theme_adrc() +
  theme(aspect.ratio = 1) +
  xlab("Number occurences\nin female") +
  ylab("Number occurences\nin male") +
  scale_size(range = c(.1, 1), breaks = c(1, 2, 5, 10, 20))

plot <-  fill_title(plot, color_v)
ggsave(
  "figures/figure_4_a.svg",
  plot,
  width = 18,
  height = 7,
  unit = "cm"
)

write.csv(counts, "SourceData/Figure_4a.csv")
