get_cells_per_sample_perc <- function(seur_obj){
  
  CpS = table(seur_obj$Sample, seur_obj$L2)
  CpS_perc = prop.table(CpS, margin=1)
  CpS_perc_df = melt(CpS_perc, varnames = c("Sample", "L2"))
  CpS_perc_df$Diagnosis = factor(seur_obj$biogroup_short[match(CpS_perc_df$Sample, seur_obj$Sample)], levels=diagnosis_order)
  
  return (CpS_perc_df)
}


getCelltypeProportion_Boxplot_scaled <- function(input_data, celltypes, y_lim = 1){
  
  input_data$L2 <- factor(input_data$L2, levels =  celltypes)
  input_data <- input_data[input_data$value < (y_lim-0.1),]
  
  p_1 = ggplot(input_data[input_data$L2 %in% celltypes,], aes(x=Diagnosis, y=value, fill=Diagnosis, color=Diagnosis)) +
    geom_boxplot(outlier.size=0.15, lwd=0.2) + 
    facet_wrap(~L2, nrow=1) +
    scale_fill_manual(values=color_v) +
    scale_color_manual(values=darken(color_v)) +
    scale_y_continuous(labels=scales::percent_format()) + 
    xlab("") + ylab("PBMCs")+
    theme_adrc() + 
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    ylim(0,y_lim)
  
  celltype_proportions_p = fill_title(p_1, color_v)
  
  return(celltype_proportions_p)
}




plot_density_difference = function(data1, data2) {
  xrng = range(data1$V1, data2$V1)
  yrng = range(data1$V2, data2$V2)
  d1 = MASS::kde2d(data1$V1, data1$V2, lims=c(xrng, yrng), n=200)
  d2 = MASS::kde2d(data2$V1, data2$V2, lims=c(xrng, yrng), n=200)

  # Calculate the difference between the 2d density estimates
  diff12 = d1
  diff12$z = d1$z - d2$z

  ## Melt data into long format
  # First, add row and column names (x and y grid values) to the z-value matrix
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y

  # Now melt it to long format
  diff12.m = melt(diff12$z, id.var=rownames(diff12))
  names(diff12.m) = c("V1","V2","z")

  diff12.m$z[diff12.m$z > 0.007] <- 0.007
  diff12.m$z[diff12.m$z < (-0.007)] <- (-0.007)
  limits = c(-0.007, 0.007)

  # Plot difference between densities
  ggplot(diff12.m, aes(V1, V2, z=z, fill=z)) +
    geom_tile() +
    stat_contour(aes(colour=..level..), binwidth=0.0005) +
    scale_fill_gradient2(name=expression(Delta*Density), low="#4575b4",mid="white", high="#d73027", midpoint=0,limits=limits, breaks=scales::pretty_breaks(5)) + 
    scale_colour_gradient2(low=scales::muted("#4575b4"), mid="white", high=scales::muted("#d73027"), midpoint=0) +
    coord_cartesian(xlim=xrng, ylim=yrng) +
    guides(colour=FALSE) + 
    ylab("") + xlab("") + 
    theme_adrc()+ theme(aspect.ratio = 1,legend.position="none")
}
