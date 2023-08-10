theme.clean <- function(){
  theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11, face = "plain"),             
          axis.title.y = element_text(size = 11, face = "plain"),             
          #panel.grid.major.x = element_blank(),                                          
          #panel.grid.minor.x = element_blank(),
          #panel.grid.minor.y = element_blank(),
          #panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 15, vjust = 1, hjust = 0.5, 
                                    face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5), 
          legend.text = element_text(size = 10, face = "plain"),          
          legend.position = "bottom")
}
