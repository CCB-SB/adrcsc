library(gghalves)
library(data.table)
library(cowplot)
source("scripts/helper.R")

meta_data <- read.csv("data/metadata.csv")
meta_data <-
  meta_data[meta_data$SCMD %in% meta_data$SCMD[meta_data$Visit != 1], ]
meta_data <-
  meta_data[, colnames(meta_data) %in% c("Sex", "Diagnosis", "SCMD", "Age", "Visit")]
meta_data$Visit <- as.factor(meta_data$Visit)

meta_data <- meta_data[!duplicated(meta_data), ]

plot <- ggplot(meta_data, aes(x = Visit, y = Age)) +
  geom_half_violin(
    aes(y = Age),
    data = meta_data[meta_data$Sex == "female", ],
    side = "l",
    color = darken(color_v["Female"]),
    fill = color_v["Female"],
    alpha = 0.8
  ) +
  geom_half_boxplot(
    aes(y = Age),
    data = meta_data[meta_data$Sex == "female", ],
    side = "l",
    fill = color_v["Female"],
    width = 0.2,
    alpha = 0.8
  ) +
  geom_half_violin(
    aes(y = Age),
    data = meta_data[meta_data$Sex == "male", ],
    side = "r",
    color = darken(color_v["Male"]),
    fill = color_v["Male"],
    alpha = 0.8
  ) +
  geom_half_boxplot(
    aes(y = Age),
    data = meta_data[meta_data$Sex == "male", ],
    side = "r",
    fill = color_v["Male"],
    width = 0.2,
    alpha = 0.8
  ) +
  theme_adrc() + ylab("Age (years)") + xlab("Visit") + theme(legend.position =
                                                               "none")
ggsave(
  "figures/Timeseries_Age.svg",
  plot,
  width = 8,
  height = 4,
  unit = "cm"
)


plot <-
  ggplot(meta_data, aes(x = Visit, y = 1)) + geom_bar(aes(fill = Diagnosis), stat =
                                                        "identity", alpha = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_adrc() + scale_fill_manual(values = color_v) + ylab("#Patients") + xlab("Visit")
ggsave(
  "figures/Timeseries_Diagnosis.svg",
  plot,
  width = 8,
  height = 4,
  unit = "cm"
)

meta_data <- read.csv("data/metadata.csv")
meta_data <-
  meta_data[meta_data$SCMD %in% meta_data$SCMD[meta_data$Visit != 1], ]
dates <- as.Date(meta_data$Date.of.blood.draw)
df <-
  merge(meta_data[meta_data$Visit == 1, ], meta_data[meta_data$Visit == 2, ], by = "SCMD", all = T)
df <-
  merge(df, meta_data[meta_data$Visit == 3, ], by = "SCMD", all = T)

data_to_plot <-
  data.frame(
    Comparison = c(
      rep("Visit 1 to Visit 2", length(df$SCMD)),
      rep("Visit 2 to Visit 3", length(df$SCMD))
    ),
    Samples = rep(df$SCMD, 2) ,
    dist = c(
      as.Date(df$Date.of.blood.draw.y) - as.Date(df$Date.of.blood.draw.x),
      as.Date(df$Date.of.blood.draw) - as.Date(df$Date.of.blood.draw.y)
    )
  )

data_to_plot$Diagnosis <-
  unlist(lapply(data_to_plot$Samples, function(sample)
    meta_data$Diagnosis[meta_data$SCMD == sample &
                          meta_data$Visit == 1][1]))
plot <-
  ggplot(data_to_plot,
         aes(
           x = Comparison,
           y = dist,
           color = Diagnosis,
           fill = Diagnosis
         )) +
  geom_boxplot() + ylab("Days between visits") +
  geom_hline(yintercept = 365,
             linetype = "dashed",
             color = "gray") +
  geom_hline(yintercept = 365 * 2,
             linetype = "dashed",
             color = "gray") +
  geom_hline(yintercept = 365 * 3,
             linetype = "dashed",
             color = "gray") + theme_adrc() + xlab("") +
  scale_color_manual(values = color_v) +
  scale_fill_manual(values = color_v)
ggsave(
  "figures/TimeBetweenVisits.svg",
  plot,
  width = 10,
  height = 5,
  unit = "cm"
)
