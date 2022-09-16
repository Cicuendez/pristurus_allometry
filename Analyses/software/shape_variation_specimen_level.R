
# Shape variation 

# Load packages
library(tidyverse)
library(broom)

# Import and log-transform morphological data
morpho_all <- read.table('data/morphology/morpho_pristurus.csv', sep = ';', 
                         dec = '.', header = TRUE)

morpho_all_log <- morpho_all %>%
  mutate(across(where(is.numeric), log10))

# Remove the effect of size ----
colnames(morpho_all_log)
morpho_all_resid <- morpho_all_log %>%
  mutate(TrL = resid(lm(TrL ~ SVL))) %>%
  mutate(HL = resid(lm(HL ~ SVL))) %>% 
  mutate(HW = resid(lm(HW ~ SVL))) %>% 
  mutate(HH = resid(lm(HH ~ SVL))) %>%
  mutate(Lhu = resid(lm(Lhu ~ SVL))) %>%
  mutate(Lun = resid(lm(Lun ~ SVL))) %>%
  mutate(Lfe = resid(lm(Lfe ~ SVL))) %>%
  mutate(Ltb = resid(lm(Ltb ~ SVL))) 

# Set the theme
theme.clean <- function(){
  theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11, face = "plain"),             
          axis.title.y = element_text(size = 11, face = "plain"),             
          #          panel.grid.major.x = element_blank(),                                          
          #          panel.grid.minor.x = element_blank(),
          #          panel.grid.minor.y = element_blank(),
          #          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 10, face = "plain"),          
          legend.position = "right")
}
theme_set(theme.clean())


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


# Perform PCA WITH size ----
pca_with_size <- morpho_all_resid %>%
  select(where(is.numeric)) %>% # retain only numeric columns
  scale() %>% # scale to zero mean and unit variance
  prcomp()

# Since we have removed the effect of size from the rest of variables, 
# when we include size in the PCA it appears as one single PC (in this case, PC3). 
# We'll keep going with the PCA without size. 

# Plot
pca_with_size %>%
  # add PCs to the original dataset
  augment(morpho_all_resid) %>%
  ggplot(aes(.fittedPC1, .fittedPC2)) +
  geom_point(aes(color = habitat_broad)) +
  scale_color_manual(values = habitat_broad_colors)

# Plot rotation matrix (no funciona)
arrow_style <- arrow(
  angle = 20, length = grid::unit(8, "pt"),
  ends = "first", type = "closed"
)

pca_with_size %>%
  # extract rotation matrix
  tidy(matrix = "rotation") %>%
  pivot_wider(
    names_from = "PC", values_from = "value",
    names_prefix = "PC"
  ) %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(
    xend = 0, yend = 0,
    arrow = arrow_style
  ) +
  geom_text(aes(label = column), hjust = 1) +
  xlim(-1.5, 0.5) + ylim(-1, 1) + 
  coord_fixed()

# Plot variance explained
pca_with_size %>%
  # extract eigenvalues
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) + 
  geom_col() + 
  scale_x_continuous(
    # create one axis tick per PC
    breaks = 1:6
  ) +
  scale_y_continuous(
    name = "variance explained",
    # format y axis ticks as percent values
    label = scales::label_percent(accuracy = 1)
  )


# Perform PCA WITHOUT size ----
pca_no_size <- morpho_all_resid %>%
  select(where(is.numeric)) %>% # retain only numeric columns
  select(-SVL) %>% # remove size
  scale() %>% # scale to zero mean and unit variance
  prcomp() # do PCA

# Loadings
pca_no_size$rotation

write.table(pca_no_size$rotation, 'objects/pca_specimens/pca_loadings.csv', 
            row.names = TRUE, 
            quote = FALSE, col.names = NA, dec = '.', sep = ';')
# PC1: limb dimensions (Lhu, Lun, Lfe, Ltb) (negative values). Longer limbs in negative values of PC1.
# PC2: head dimensions (HL, HW, HH) (positive values). Larger heads in positive values of PC2.
# We'll keep PC1 and PC2. 

# Variance explained
summary(pca_no_size)
pca_summary <- summary(pca_no_size)
pca_summary$importance
par(mar = c(3,3, 2,2))
pdf('plots/morphospace_specimen/variance_explained.pdf')
barplot(pca_summary$importance[2,])
title('Variance explained by specimen PCA')
dev.off()

# Plot
pca_no_size %>%
  # add PCs to the original dataset
  augment(morpho_all_resid) %>%
  ggplot(aes(.fittedPC1, .fittedPC2)) +
  geom_point(aes(color = habitat_broad)) +
  scale_color_manual(values = habitat_broad_colors)

# Plot rotation matrix (no funciona)
arrow_style <- arrow(
  angle = 20, length = grid::unit(8, "pt"),
  ends = "first", type = "closed"
)

pca_no_size %>%
  # extract rotation matrix
  tidy(matrix = "rotation") %>%
  pivot_wider(
    names_from = "PC", values_from = "value",
    names_prefix = "PC"
  ) %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(
    xend = 0, yend = 0,
    arrow = arrow_style
  ) +
  geom_text(aes(label = column), hjust = 1) +
  xlim(-1.5, 0.5) + ylim(-1, 1) + 
  coord_fixed()

# Plot variance explained
pca_no_size %>%
  # extract eigenvalues
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) + 
  geom_col(fill = 'gray70') + 
  scale_x_continuous(
    # create one axis tick per PC
    breaks = 1:8
  ) +
  scale_y_continuous(
    name = "variance explained",
    # format y axis ticks as percent values
    label = scales::label_percent(accuracy = 1)
  ) +
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 1), 
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.line.y = element_line(colour = "black", size = 0.3), 
        axis.ticks.y = element_line(),
        axis.ticks.length = unit(.2, "cm"))

# PLOT MORPHOSPACE (PC1 and PC2) ----
morphospace_land <- pca_no_size %>%
  # add PCs to the original dataset
  augment(morpho_all_resid) %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = land)) +
  scale_color_manual(values = land_colors) +
  labs(x = "PC1", y = "PC2")

morphospace_habitat <- pca_no_size %>%
  # add PCs to the original dataset
  augment(morpho_all_resid) %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = habitat)) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "PC1", y = "PC2")

morphospace_habitat_broad <- pca_no_size %>%
  # add PCs to the original dataset
  augment(morpho_all_resid) %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = habitat_broad)) +
  scale_color_manual(values = habitat_broad_colors) +
  labs(x = "PC1", y = "PC2")


# BOXPLOTS SHAPE (PC1 AND PC2) ----
pca_table <- pca_no_size %>%
  # add PCs to the original dataset
  augment(morpho_all_resid)


saveRDS(pca_table, "objects/pca_specimens/pca_table_specimens.rds")




# PC1 ----

# Land
namevec_land <- arrange(pca_table, land, .fittedPC1) %>%
  pull(species) 
pca_table$species <- factor(pca_table$species, level = unique(namevec_land))

boxplot_pc1_land <- pca_table %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = land), color = 'black') +
  scale_fill_manual(values = land_colors) +
  theme(
    axis.text.x = element_text(angle = 45)
  ) +
  labs(y = "PC1 (42%)")

# Habitat broad
namevec_habitat_broad <- arrange(pca_table, habitat_broad, .fittedPC1) %>%
  pull(species) 
pca_table$species <- factor(pca_table$species, level = unique(namevec_habitat_broad))

boxplot_pc1_habitat_broad <- pca_table %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = habitat_broad), color = 'black') +
  scale_fill_manual(values = habitat_broad_colors) +
  theme(
    axis.text.x = element_text(angle = 45)
    ) +
   labs(y = "PC1 (42%)")

# Habitat
namevec_habitat <- arrange(pca_table, habitat, .fittedPC1) %>%
  pull(species) 
pca_table$species <- factor(pca_table$species, level = unique(namevec_habitat))

boxplot_pc1_habitat <- pca_table %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = habitat), color = 'black') +
  scale_fill_manual(values = habitat_colors) +
  theme(
    axis.text.x = element_text(angle = 45)
  ) +
  labs(y = "PC1 (42%)")

# PC2 ----

# Land
namevec_land <- arrange(pca_table, land, .fittedPC2) %>%
  pull(species) 
pca_table$species <- factor(pca_table$species, level = unique(namevec_land))

boxplot_pc2_land <- pca_table %>%
  ggplot(aes(x = species, y = .fittedPC1)) +
  geom_boxplot(aes(fill = land), color = 'black') +
  scale_fill_manual(values = land_colors) +
  theme(
    axis.text.x = element_text(angle = 45)
  ) +
  labs(y = "PC2 (17%)") +
  coord_flip()

# Habitat broad
namevec_habitat_broad <- arrange(pca_table, habitat_broad, .fittedPC2) %>%
  pull(species) 
pca_table$species <- factor(pca_table$species, level = unique(namevec_habitat_broad))

boxplot_pc2_habitat_broad <- pca_table %>%
  ggplot(aes(x = species, y = .fittedPC2)) +
  geom_boxplot(aes(fill = habitat_broad), color = 'black') +
  scale_fill_manual(values = habitat_broad_colors) +
  theme(
    axis.text.x = element_text(angle = 45)
  ) +
  labs(y = "PC2 (17%)")

# Habitat
namevec_habitat <- arrange(pca_table, habitat, .fittedPC2) %>%
  pull(species) 
pca_table$species <- factor(pca_table$species, level = unique(namevec_habitat))

boxplot_pc2_habitat <- pca_table %>%
  ggplot(aes(x = species, y = .fittedPC2)) +
  geom_boxplot(aes(fill = habitat), color = 'black') +
  scale_fill_manual(values = habitat_colors) +
  theme(
    axis.text.x = element_text(angle = 45)
  ) +
  labs(y = "PC2 (17%)")


# Save plots ----
ggsave("plots/morphospace_specimen/morphospace_land_specimen.pdf", plot = morphospace_land)
ggsave("plots/morphospace_specimen/morphospace_habitat_specimen.pdf", plot = morphospace_habitat)
ggsave("plots/morphospace_specimen/morphospace_habitat_broad_specimen.pdf", plot = morphospace_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_pc1_land_specimen.pdf", plot = boxplot_pc1_land)
ggsave("plots/boxplots_specimen/boxplot_pc1_habitat_broad_specimen.pdf", plot = boxplot_pc1_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_pc1_habitat_specimen.pdf", plot = boxplot_pc1_habitat)
ggsave("plots/boxplots_specimen/boxplot_pc2_land_specimen.pdf", plot = boxplot_pc2_land)
ggsave("plots/boxplots_specimen/boxplot_pc2_habitat_broad_specimen.pdf", plot = boxplot_pc2_habitat_broad)
ggsave("plots/boxplots_specimen/boxplot_pc2_habitat_specimen.pdf", plot = boxplot_pc2_habitat)



################
################
################

# Plot without scaling
morpho_all_resid %>%
  ggplot() +
  aes(SVL, TrL) +
  geom_point(aes(color = habitat_broad), alpha = 0.9, shape = 16) +
  scale_color_manual(values = habitat_broad_colors) 


# Plot with scaling
morpho_all_resid %>%
  # scale all numeric columns
  mutate(across(where(is.numeric), scale)) %>%
  ggplot() +
  aes(SVL, TrL) +
  geom_point(aes(color = habitat), alpha = 0.9, shape = 16) +
  scale_color_manual(values = habitat_colors)







