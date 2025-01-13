library(ggrepel)

supp_data <- read.csv("Literature.csv")
this_data = data.frame(Paper = "This paper", DOI = "", method = "sc", total.nb.of.samples= 290, nb.of.patients = 0, nb.of.controls = 0, nb.of.cells = 909000, Disease = "AD, PD, MCI", first.author..short. = "This Dataset", year = "")
supp_data <- rbind(supp_data, this_data)
supp_data$nb.of.cells[supp_data$nb.of.cells == 0] <- 1
supp_data$Disease <- factor(supp_data$Disease, levels = c("AD","PD","AD, PD", "AD, MCI","AD, PD, MCI"))
colors <- c("#ff7f00", "#6a3d9a", "#1f78b4", "#a6cee3", "black")
plot <- ggplot(supp_data, aes(x = total.nb.of.samples, y =nb.of.cells, color = Disease, label = paste(first.author..short., year, sep = ", ")))+
  geom_point()+ 
  geom_text_repel(size = 2, angle = 45, hjust = -0.1)+ 
  scale_y_continuous(limits = c(1, 10000000),trans = "log", breaks = c(100, 1000, 10000, 100000, 1000000), labels = c("100", "1k", "10k", "100k", "1m"))+ 
  scale_x_continuous(limits = c(-80, 1300), breaks = c(0,500, 1000))+
  geom_hline(yintercept = 10)+ 
  geom_vline(xintercept = 1240)+ 
  theme_adrc()+ 
  xlab("Number of samples")+ 
  ylab("Number of cells")+ 
  annotate("text",x = 1290, y = 5000, label = "Single-cell", angle = 270, color = "black")+ 
  annotate("text",x = 1290, y = 2, label = "Bulk", angle = 270, color = "black")+ 
  scale_color_manual(values = colors)
ggsave("figures/Literature_Overvew.svg", plot, width = 5, 
