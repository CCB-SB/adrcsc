library(circlize)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(cowplot)
source("scripts/helper.R")

df_to_matrix <- function(df, x_col, y_col, data_col) {
  data <-
    sapply(unique(df[[x_col]]), function(x)
      sapply(unique(df[[y_col]]), function(y) {
        ret <-
          as.numeric(as.character(df[[data_col]][df[[x_col]] == x &
                                                   df[[y_col]] == y][1]))
        if (is.na(ret)[1]) {
          ret = 0
        }
        ret
      }))
  
  return(data)
}

plot_heatmap <- function(dtp, lim = 20) {
  df <- df_to_matrix(dtp, "X.Name", "cluster_id", "color")
  mat <- t(as.matrix(df))
  mat <- mat[rowSums(mat != 0)  > 3 , colSums(mat != 0) > 4]
  
  colors = colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027")) 
  
  plot <-
    ComplexHeatmap::Heatmap(
      t(mat),
      name = "mat",
      row_names_side = "left",
      col = colors,
      column_names_rot = 45,
      column_km = 3,
      row_km = 3,
      show_row_dend = FALSE,
      show_column_dend = FALSE
    )
  return(plot)
}

pathways_female <-
  read.csv("results_female/pathway_analysis_results_all.csv")
pathways_female$sex <- "female"
pathways_male <-
  read.csv("results_male/pathway_analysis_results_All.csv")
pathways_male$sex <- "male"

pw_collected <- rbind(pathways_female, pathways_male)

pw_count <-
  pw_collected %>% group_by(X.Name, contrast, sex) %>% summarise(PW_count = n())


pw_tmp <-
  pw_count[pw_count$contrast == "ADvsHC" & pw_count$sex == "female", ]
ad_pw_f <-
  pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
pw_tmp <-
  pw_count[pw_count$contrast == "ADvsHC" & pw_count$sex == "male", ]
ad_pw_m <-
  pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
pw_tmp <-
  pw_count[pw_count$contrast == "PDvsHC" & pw_count$sex == "female", ]
pd_pw_f <-
  pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
pw_tmp <-
  pw_count[pw_count$contrast == "PDvsHC" & pw_count$sex == "male", ]
pd_pw_m <-
  pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
pw_tmp <-
  pw_count[pw_count$contrast == "MCIvsHC" & pw_count$sex == "female", ]
mci_pw_f <-
  pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]
pw_tmp <-
  pw_count[pw_count$contrast == "MCIvsHC" & pw_count$sex == "male", ]
mci_pw_m <-
  pw_tmp$X.Name[order(pw_tmp$PW_count, decreasing = T)][1:5]

pw_col <- c(ad_pw_f, ad_pw_m, pd_pw_f, pd_pw_m, mci_pw_f, mci_pw_m)
df_pw <-
  data.frame(
    X.Name = pw_col,
    contrast = c(rep("ADvsHC", 10), rep("PDvsHC", 10), rep("MCIvsHC", 10)),
    sex = rep(c(rep("female", 5), rep("male", 5)), 3)
  )
df_pw$amongTop5 <- T
df <- pw_collected[pw_collected$X.Name %in% df_pw$X.Name, ]
summary <-
  df %>% group_by(X.Name, contrast, sex, Regulation_direction) %>% summarise(n = n())
merged  <-
  merge(df_pw,
        summary,
        all = T,
        by = c("X.Name", "contrast", "sex"))
merged$amongTop5[is.na(merged$amongTop5)] <- F

merged$X.Name <-
  unlist(lapply(as.character(merged$X.Name), function(string) {
    n = unlist(gregexpr(pattern = ' ', string))
    
    if (any(n > 50)) {
      pos <- n[n > 50][1]
      substr(string, pos, pos) <- "\n"
    }
    
    if (any(n > 25)) {
      pos <- n[n > 25][1]
      substr(string, pos, pos) <- "\n"
    }
    
    string
  }))

merged$n_orig <- merged$n

merged$n <- ifelse(merged$Regulation_direction, merged$n,-merged$n)
merged$n[merged$n == 0] <- NA
plot <-
  ggplot(merged,
         aes(
           x = as.character(Regulation_direction),
           y = X.Name,
           fill = sign(n) * log10(abs(n)),
           label = ifelse(abs(merged$n_orig) > 10, merged$n_orig, "")
         )) +
  geom_point(shape = ifelse(merged$Regulation_direction == 1, 95, 43))  +
  geom_tile()  +
  scale_x_discrete(labels = c("Enr.", "Depl."),
                   breaks = c("1", "0")) +
  facet_nested(. ~ contrast + sex) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "#Occurences",
    breaks = c(
      -log10(5),
      -log10(2),
      0,
      log10(2),
      log10(5),
      log10(10),
      log10(20),
      log10(30)
    ),
    labels =  c(5, 2, 0, 2, 5, 10, 20, 30)
  ) +
  theme_adrc() + xlab("") + ylab("") + theme(
    axis.text.y = element_text(lineheight = 0.7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
plot <-  fill_title(plot, color_v)
ggsave(
  "figures/figure_4_b.svg",
  plot,
  width = 12,
  height = 9,
  unit = "cm"
)

write.csv(merged, "SourceData/figure_4_b.csv")
