library(ggplot2)
library(data.table)
library(gghalves)
library(cowplot)
library(ggplotify)

source("scripts/helper.R")

meta_data <- read.csv("data/metadata.csv")

meta_data <- meta_data[,!colnames(meta_data) %in% c("X", "Date.Processed.at.CG")]
meta_data <- meta_data[!duplicated(meta_data),]

meta_data <- as.data.table(meta_data)
name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease"="PD", "Parkinson's Disease only"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")
meta_data[, Diagnosis_short:=name2short[Diagnosis]]
meta_data[, Diagnosis_short:=factor(Diagnosis_short, levels=diagnosis_order)]

pat_tbl = meta_data[meta_data$Visit == 1,]
pat_tbl$Diagnosis_short <- factor(pat_tbl$Diagnosis_short, levels = diagnosis_order)

plot_df2 = pat_tbl[, list(Frequency=.N), by=c("Diagnosis_short", "Sex")][order(Diagnosis_short, -Sex)]
plot_df2 = plot_df2[, list(Sex, Frequency, pos=cumsum(Frequency)-0.5*Frequency), by=Diagnosis_short]

############# Figure 1a)

plot <- ggplot(meta_data, aes(x = Visit, fill = Diagnosis)) + 
  geom_bar()+ 
  theme_classic()+ 
  scale_fill_manual(values = color_v)+
  theme(legend.position = "none")+ 
  xlab("Timepoint")+ ylab("Number of Samples")

ggsave("figures/figure_1_a.svg", plot, width = 3, height = 2)
write.csv(table(meta_data$Visit, meta_data$Diagnosis), "SourceData/Figure1a.csv")

############# Figure 1b)

fig_1_b = ggplot(plot_df2, aes(x=factor(Diagnosis_short, levels = diagnosis_order), y = Frequency)) + 
  geom_bar(aes(fill=Sex), stat="identity", alpha = 0.8) +
  geom_text(aes(label=Frequency, y=pos), size=2) +
  scale_y_continuous(expand = expansion(mult=c(0, 0))) +
  theme_adrc() + 
  scale_fill_manual(values=color_v) + 
  ylab("#Patients") + xlab("") 

legend <- get_legend(fig_1_b)

fig_1_b <- fig_1_b + 
  theme(legend.position="none", legend.direction = "vertical") 

ggsave("figures/figure_1_b_1.svg", fig_1_b, width = 7, height = 3.5, unit = "cm")
ggsave("figures/legend_fig1_b_1.svg", as.ggplot(legend), width = .5, height = .5)

############# Figure 1c)

fig_1_c = ggplot(pat_tbl, aes(x=Diagnosis_short, y = Age)) + 
  geom_half_violin(aes(y=Age), data=pat_tbl[Sex == "female"], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + 
  geom_half_boxplot(aes(y=Age), data=pat_tbl[Sex == "female"], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8) +
  geom_half_violin(aes(y=Age), data=pat_tbl[Sex == "male"], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) +
  geom_half_boxplot(aes(y=Age), data=pat_tbl[Sex == "male"], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8) +
  theme_adrc() + 
  ylab("Age (years)") + xlab("") + 
  theme(legend.position="none")

ggsave("figures/figure_1_b_2.svg", fig_1_c, width = 7, height = 3.5, unit = "cm")
