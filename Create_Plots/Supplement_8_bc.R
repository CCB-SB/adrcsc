library(ggplot2)
library(cowplot)
library(aplot)
library(ggpubr)
library(viridis)
library(data.table)
source("scripts/helper.R")

pw1 <-
  read.csv("data/Genetrail_Male_BrainBlood/GO_-_Biological_Process.tsv",
           sep = "\t")
pw1$Database <- "GO - BP"
pw2 <-
  read.csv("data/Genetrail_Male_BrainBlood/GO_-_Cellular_Component.tsv",
           sep = "\t")
pw2$Database <- "GO - CC"
pw3 <-
  read.csv("data/Genetrail_Male_BrainBlood/GO_-_Molecular_Function.tsv",
           sep = "\t")
pw3$Database <- "GO - MF"
pw4 <-
  read.csv("data/Genetrail_Male_BrainBlood/KEGG_-_Pathways.tsv",
           sep = "\t")
pw4$Database <- "KEGG"


pw <- rbind(pw1, pw2, pw3, pw4)
pw <- pw[pw$P.value < 0.05, ]

pw_plot_m <- pw[order(pw$P.value, decreasing = F), ]

pw$gene_score <-
  ifelse(pw$Regulation_direction,
         pw$Expected.Score,
         -pw$Expected.Score)
plot <-
  ggplot(pw, aes(
    y = -log10(P.value),
    x = gene_score,
    color = Database,
    label = X.Name
  )) + geom_point() +
  scale_color_manual(values = viridis(4)) + theme_adrc() +
  ylab("-log10 adj. p-value") + xlab("Gene Score")
ggsave(
  "figures/Supplement_8_b.svg",
  plot,
  width = 6,
  height = 5,
  unit = "cm"
)

write.csv(pw, "SourceData/Supplement_8b.csv")

pw1 <-
  read.csv("data/Genetrail_Female_BrainBlood/KEGG_-_Pathways.tsv",
           sep = "\t")
pw1$Database <- "KEGG"

pw <- pw1
pw <- pw[pw$P.value < 0.05, ]
pw_plot_f <- pw[order(pw$P.value, decreasing = F), ]

pw_plot_m$sex <- "Male"
pw_plot_f$sex <- "Female"

pw_plot <- rbind(pw_plot_f, pw_plot_m)

pw_plot$label <-
  unlist(lapply(as.character(pw_plot$X.Name), function(string) {
    n = unlist(gregexpr(pattern = ' ', string))
    if (any(n > 30)) {
      pos <- n[n > 30][1]
      substr(string, pos, pos) <- "\n"
    }
    string
  }))
pw_plot$label <-
  factor(pw_plot$label, levels = pw_plot$label[order(pw_plot$P.value, decreasing = T)])
plot <-
  ggplot(pw_plot, aes(
    x = -log10(P.value),
    y = label,
    fill = Database
  )) +
  geom_bar(stat = "identity") + theme_adrc() + facet_grid(sex ~ ., scales = "free", space = "free") +
  scale_fill_manual(values = viridis(4)) + theme_adrc() +
  xlab("-log10 adj. p-value") + ylab("") + theme(legend.position = "none")

plot <- fill_title(plot, color_v)
ggsave(
  "figures/Supplement_8_c.svg",
  plot,
  width = 9,
  height = 8,
  unit = "cm"
)

write.csv(pw_plot, "SourceData/Supplement_8c.csv")