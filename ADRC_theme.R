
theme_adrc <- function(){ 
    font <- "sans"   #assign font family up front
    
    theme_cowplot(10) %+replace%    #replace elements we want to change
    
    theme(
      #text elements
      plot.title = element_text(             #title
                   family = font,            #set font family
                   size = 20,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0,                #left align
                   vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
                   family = font,            #font family
                   size = 14),               #font size
      
      plot.caption = element_text(           #caption
                   family = font,            #font family
                   size = 9,                 #font size
                   hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
                   family = font,            #font family
                   size = 8), #10              #font size
      
      axis.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = 6),                #font size
      
      strip.text.x = element_text( 
                   margin = margin(1,0,1,0, "mm"),             #axis text
                   family = font,            #axis famuly
                   size = 6),
      strip.text.y = element_text( 
                   margin = margin(0,1,0,1, "mm"),             #axis text
                   family = font,            #axis famuly
                   size = 6, angle = -90),
      
      legend.title=element_text(size=8),
      legend.text=element_text(size=6), 
      
      panel.spacing = unit(0.1,"line")
    )
}
