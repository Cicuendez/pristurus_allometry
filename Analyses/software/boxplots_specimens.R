# Boxplots size and shape specimens 

libs <- c('tidyverse', 'phytools', 'broom', 'patchwork')
lapply(libs, require, character.only = TRUE)


# Import morpho ----
pca_table_specimens <- readRDS("objects/pca_specimens/pca_table_specimens.rds")
sort(table(pca_table_specimens$species))
?desc
nrow(pca_table_specimens)

# Set the theme ----
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


# Set the colors ----
habitat_broad_colors0 <- c(ground = 'brown', rock = 'gray', tree = 'darkgreen')
habitat_broad_colors <- c(ground = '#e0710b', rock = '#6C586E', tree = '#119616')
habitat_broad_colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")
plot(1:3, 1:3, col = habitat_broad_colors, cex = 10, pch = 16)

habitat_colors0 <- c('hard-ground' = 'brown', 'soft-ground' = 'orange', 
                     rock = 'gray', tree = 'darkgreen')
habitat_colors <- c('hard-ground' = '#D63916', 'soft-ground' = '#E9A800', 
                    rock = '#6C586E', tree = '#119616')

land_colors0 <- c(mainland = 'coral2', island = 'seagreen')
land_colors1 <- c(mainland = '#A54E29', island = '#42978A')
land_colors2 <- c(mainland = '#C56542', island = '#67D1B7')
land_colors <- c(mainland = '#C56542', island = '#42978A')


# SVL ----
# Boxplot habitat_broad
namevec_habitat_broad <- arrange(pca_table_specimens, habitat_broad, SVL) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_habitat_broad))

boxplot_svl_habitat_broad <- ggplot(pca_table_specimens) +
  geom_boxplot(aes(x = species, y = SVL, fill = habitat_broad)) +
  scale_fill_manual(values = habitat_broad_colors) +
  theme.clean() +
  #  theme(legend.position = 'none') +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  coord_flip()

# Boxplot habitat
namevec_svl_habitat <- arrange(pca_table_specimens, habitat, SVL) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_svl_habitat))

boxplot_svl_habitat <- ggplot(pca_table_specimens) +
  geom_boxplot(aes(x = species, y = SVL, fill = habitat)) +
  scale_fill_manual(values = habitat_colors) +
  theme.clean() +
  #  theme(legend.position = 'none') +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  coord_flip()

# Boxplot land
namevec_land <- arrange(pca_table_specimens, land, SVL) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_land))

boxplot_svl_land <- ggplot(pca_table_specimens) +
  geom_boxplot(aes(x = species, y = SVL, fill = land)) +
  scale_fill_manual(values = land_colors) +
  theme.clean() +
  #  theme(legend.position = 'none') +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  coord_flip()

# PC1 ----

# Land
namevec_land <- arrange(pca_table_specimens, land, .fittedPC1) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_land))

boxplot_pc1_land <- pca_table_specimens %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = land)) +
  scale_fill_manual(values = land_colors) +
  theme.clean() +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  labs(y = "PC1 (42%)") +
  coord_flip()

# Habitat broad
namevec_habitat_broad <- arrange(pca_table_specimens, habitat_broad, .fittedPC1) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_habitat_broad))

boxplot_pc1_habitat_broad <- pca_table_specimens %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = habitat_broad)) +
  scale_fill_manual(values = habitat_broad_colors) +
  theme.clean() +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  labs(y = "PC1 (42%)") +
  coord_flip()

# Habitat
namevec_habitat <- arrange(pca_table_specimens, habitat, .fittedPC1) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_habitat))

boxplot_pc1_habitat <- pca_table_specimens %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = habitat)) +
  scale_fill_manual(values = habitat_colors) +
  theme.clean() +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  labs(y = "PC1 (42%)") +
  coord_flip()

# PC2 ----

# Land
namevec_land <- arrange(pca_table_specimens, land, .fittedPC2) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_land))

boxplot_pc2_land <- pca_table_specimens %>%
  ggplot(aes(x = species, y = .fittedPC2)) +
  geom_boxplot(aes(fill = land)) +
  scale_fill_manual(values = land_colors) +
  theme.clean() +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  labs(y = "PC2 (17%)") +
  coord_flip()

# Habitat broad
namevec_habitat_broad <- arrange(pca_table_specimens, habitat_broad, .fittedPC2) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_habitat_broad))

boxplot_pc2_habitat_broad <- pca_table_specimens %>%
  ggplot(aes(x = species, y = .fittedPC2)) +
  geom_boxplot(aes(fill = habitat_broad)) +
  scale_fill_manual(values = habitat_broad_colors) +
  theme.clean() +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  labs(y = "PC2 (17%)") +
  coord_flip()

# Habitat
namevec_habitat <- arrange(pca_table_specimens, habitat, .fittedPC2) %>%
  pull(species) 
pca_table_specimens$species <- factor(pca_table_specimens$species, level = unique(namevec_habitat))

boxplot_pc2_habitat <- pca_table_specimens %>%
  ggplot(aes(x = species, y = .fittedPC2)) +
  geom_boxplot(aes(fill = habitat)) +
  scale_fill_manual(values = habitat_colors) +
  theme.clean() +
  theme(
    axis.text.x = element_text(angle = 0)
  ) +
  labs(y = "PC2 (17%)") +
  coord_flip()


# Plot size and shape (pc1 and pc2) for both land and habitat ----
boxplots_specimens <- 
(boxplot_svl_land + boxplot_pc1_land + boxplot_pc2_land) / 
  (boxplot_svl_habitat_broad + boxplot_pc1_habitat_broad + boxplot_pc2_habitat_broad) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', 
        legend.title = element_blank())

ggsave('plots/boxplots_specimen/boxplot_all_specimens.pdf', 
       plot = boxplots_specimens,
       width = 15, height = 10)


# Save plots ----
ggsave("plots/boxplots_specimen/boxplot_svl_land_specimen.pdf", plot = boxplot_svl_land)
ggsave("plots/boxplots_specimen/boxplot_svl_habitat_broad_specimen.pdf", plot = boxplot_svl_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_svl_habitat_specimen.pdf", plot = boxplot_svl_habitat)
ggsave("plots/morphospace_specimen/morphospace_land_specimen.pdf", plot = morphospace_land)
ggsave("plots/morphospace_specimen/morphospace_habitat_specimen.pdf", plot = morphospace_habitat)
ggsave("plots/morphospace_specimen/morphospace_habitat_broad_specimen.pdf", plot = morphospace_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_pc1_land_specimen.pdf", plot = boxplot_pc1_land)
ggsave("plots/boxplots_specimen/boxplot_pc1_habitat_broad_specimen.pdf", plot = boxplot_pc1_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_pc1_habitat_specimen.pdf", plot = boxplot_pc1_habitat)
ggsave("plots/boxplots_specimen/boxplot_pc2_land_specimen.pdf", plot = boxplot_pc2_land)
ggsave("plots/boxplots_specimen/boxplot_pc2_habitat_broad_specimen.pdf", plot = boxplot_pc2_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_pc2_habitat_specimen.pdf", plot = boxplot_pc2_habitat)
