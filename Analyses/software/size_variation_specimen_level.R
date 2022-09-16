
# Load packages
library(tidyverse)

# Import and log-transform morphological data
morpho_all <- read.table('data/morphology/morpho_pristurus.csv', sep = ';', 
                         dec = '.', header = TRUE)

morpho_all_log <- morpho_all %>%
  mutate(across(where(is.numeric), log10))


# Set the theme
theme.clean <- function(){
  theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "plain"),             
          axis.title.y = element_text(size = 14, face = "plain"),             
          #          panel.grid.major.x = element_blank(),                                          
          #          panel.grid.minor.x = element_blank(),
          #          panel.grid.minor.y = element_blank(),
          #          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 12, face = "plain"),          
          legend.position = "right")
}


# Set the colors
habitat_broad_colors0 <- c(ground = 'brown', rock = 'gray', tree = 'darkgreen')
habitat_broad_colors <- c(ground = '#e0710b', rock = '#6C586E', tree = '#119616')
habitat_colors0 <- c('hard-ground' = 'brown', 'soft-ground' = 'orange', 
                     rock = 'gray', tree = 'darkgreen')
habitat_colors <- c('hard-ground' = '#D63916', 'soft-ground' = '#E9A800', 
                    rock = '#6C586E', tree = '#119616')
land_colors0 <- c(mainland = 'coral2', island = 'seagreen')
land_colors1 <- c(mainland = '#A54E29', island = '#42978A')
land_colors2 <- c(mainland = '#C56542', island = '#67D1B7')
land_colors <- c(mainland = '#C56542', island = '#42978A')



# ORDERED BOXPLOTS ----

# Boxplot habitat_broad
namevec_habitat_broad <- arrange(morpho_all_log, habitat_broad, SVL) %>%
  pull(species) 
morpho_all_log$species <- factor(morpho_all_log$species, level = unique(namevec_habitat_broad))

boxplot_habitat_broad <- ggplot(morpho_all_log) +
  geom_boxplot(aes(x = species, y = SVL, fill = habitat_broad)) +
  scale_fill_manual(values = habitat_broad_colors) +
  theme.clean() +
  #  theme(legend.position = 'none') +
  theme(legend.position = 'right')

# Boxplot habitat
namevec_habitat <- arrange(morpho_all_log, habitat, SVL) %>%
  pull(species) 
morpho_all_log$species <- factor(morpho_all_log$species, level = unique(namevec_habitat))

boxplot_habitat <- ggplot(morpho_all_log) +
  geom_boxplot(aes(x = species, y = SVL, fill = habitat)) +
  scale_fill_manual(values = habitat_colors) +
  theme.clean() +
  #  theme(legend.position = 'none') +
  theme(legend.position = 'right')

# Boxplot land
namevec_land <- arrange(morpho_all_log, land, SVL) %>%
  pull(species) 
morpho_all_log$species <- factor(morpho_all_log$species, level = unique(namevec_land))

boxplot_land <- ggplot(morpho_all_log) +
  geom_boxplot(aes(x = species, y = SVL, fill = land)) +
  scale_fill_manual(values = land_colors) +
  theme.clean() +
  #  theme(legend.position = 'none') +
  theme(legend.position = 'right')


# Save plots ----
ggsave("plots/boxplots_specimen/boxplot_svl_land_specimen.pdf", plot = boxplot_land)
ggsave("plots/boxplots_specimen/boxplot_svl_habitat_broad_specimen.pdf", plot = boxplot_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_svl_habitat_specimen.pdf", plot = boxplot_habitat)






