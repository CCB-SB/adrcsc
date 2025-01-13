df <- read.csv(file = snakemake@input[[1]])

df$x <- ifelse(df$x == "Alzheimer's disease", "AD", "HC")

write.csv(df$x, file = snakemake@output[[1]])