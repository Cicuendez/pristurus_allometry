


libs <- c('tidyverse', 'phytools', 'broom')
lapply(libs, require, character.only = TRUE)

par(mfrow=c(1,1))

# Import tree(s)
tree <- read.nexus('data/phylogeny/pristurus_tree_final.nex')


# Import morpho ----
morpho <- read.table('data/morphology/morpho_sp_final.csv', sep=';', dec='.', 
                     header = TRUE)
rownames(morpho) <- morpho$species

# log10-transform morpho data
morpho_log <- morpho %>%
  mutate(across(where(is.numeric), log10))
rownames(morpho_log) <- morpho_log$species
as_tibble(morpho_log)


# Set the theme and colors ----
# Theme
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

# Colors
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
land_colors <- c(mainland = '#F37338', island = '#801638')
land_colors <- c(mainland = '#C56542', island = '#42978A')
plot(1:2, 1:2, col = land_colors, cex = 6, pch = 16)

# Phylogenetic PCA ----
# Phylogenetic regression and residuals
?phyl.resid
svl <- morpho_log$SVL
names(svl) <- rownames(morpho_log)
phyresid_bm <- phyl.resid(tree = tree, x = svl, Y = morpho_log[,7:14], method = 'BM')
phyresid_lambda <- phyl.resid(tree = tree, x = svl, Y = morpho_log[,7:14], method = 'lambda')
# We use lambda residuals, where the regressions are made after scaling the phylogeny 
# instead of assuming Brownian motion. 



# Phylogenetic PCA of the residuals
?phyl.pca
# We perform phylogenetic PCA on the variance-covariance matrix with the method 'lambda'.
phypca <- phyl.pca(tree = tree, Y = phyresid_lambda$resid, method = 'lambda', mode = 'cov')
phypca$Eval
phypca$Evec
phypca$S
phypca$L

# Variance explained
phypca_summary <- summary(phypca)
par(mar = c(3,3, 2,2))
pdf('plots/morphospace_species/variance_explained.pdf')
barplot(phypca_summary$importance[2,])
title('Variance explained by pPCA')
dev.off()
# PC1: 0.61
# PC2: 0.16
# PC3: 0.07
# We keep PC1 and PC2. 

# Loadings
phypca$L
write.table(phypca$L, 'objects/phypca/pca_loadings.csv', row.names = TRUE, 
            quote = FALSE, col.names = NA, dec = '.', sep = ';')
# PC1: limb dimensions (Lhu, Lun, Lfe, Ltb) (negative values). Longer limbs in negative values of PC1.
# PC2: head dimensions (HL, HW, HH) (negative values). Larger heads in negative values of PC2.

# PLOT MORPHOSPACE (PC1 and PC2) ----
phypca$S
rownames(phypca$S)
rownames(morpho_log)
identical(rownames(morpho_log), rownames(phypca$S))


phypca_scores <- cbind(morpho_log[, 1:6], phypca$S)
write.table(phypca_scores, 'objects/phypca/phypca_scores.csv', sep = ';', dec = '.', 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


morphospace_land_phy <- ggplot(phypca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = land, shape = land), size = 4) +
  scale_color_manual(values = land_colors) +
  labs(x = "PC1", y = "PC2")

morphospace_habitat_broad_phy <- ggplot(phypca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = habitat_broad, shape = habitat_broad), size = 4) +
  scale_color_manual(values = habitat_broad_colors) +
  labs(x = "PC1", y = "PC2")

morphospace_habitat_phy <- ggplot(phypca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = habitat, shape = habitat), size = 4) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "PC1", y = "PC2")

# Plot grid
three_morphospace_phy_plotlist <- list(morphospace_land_phy, 
                                       morphospace_habitat_phy,
                                       morphospace_habitat_broad_phy)

three_morphospace_phy <- cowplot::plot_grid(plotlist =three_morphospace_phy_plotlist, 
                                            ncol = 3, 
                                            labels = 'auto')
ggsave('plots/phylomorphospace/morphospace_phy.pdf', three_morphospace_phy, paper = 'a4r',
       height = 3, width = 15)


# Phylomorphospace (3 panels: insularity, habitat, habitat broad) ----
?phylomorphospace
phylomorphospace(tree = tree, X = phypca_scores[, 7:8], label = 'off')

# Import trees mapped with discrete characters
pd_habitat_cons <- readRDS("objects/anc_rec/pd_habitat_cons.rds")
pd_habitat_broad_cons <- readRDS("objects/anc_rec/pd_habitat_broad_cons.rds")
pd_land_cons <- readRDS("objects/anc_rec/pd_land_cons.rds")

pdf(file = 'plots/phylomorphospace/phylomorphospace_all.pdf', width = 15, height = 5)
par(mfrow=c(1, 3))

# Mainland/island
phymorphospace_land <- phylomorphospace(tree = pd_land_cons$tree[[1]], X = phypca_scores[, 7:8], label = 'off', 
                 colors = land_colors, bty = 'n', ftype = 'off', node.by.map = TRUE,
                 node.size = c(0, 1.7), 
                 xlab = 'PC1 (-limb dimensions)', 
                 ylab = 'PC2 (-head dimensions)')
points(x = c(-0.14, -0.14), y = c(0.04, 0.05), col = land_colors, cex = 2, pch = 16)
text(x = c(-0.135, -0.135), y = c(0.04, 0.05), labels = c('mainland', 'island'), adj = 0)
#mtext('Phylomorphospace of body shape in Pristurus', side = 3, line = 1, cex = 1, font = 2)
#title(main="Phylomorphospace of body shape in Pristurus", font.main=1)

# Habitat
phymorphospace_habitat <- phylomorphospace(tree = pd_habitat_cons$tree[[1]], X = phypca_scores[, 7:8], label = 'off', 
                                        colors = habitat_colors, bty = 'n', ftype = 'off', node.by.map = TRUE,
                                        node.size = c(0, 1.7), 
                                        xlab = 'PC1 (-limb dimensions)', 
                                        ylab = 'PC2 (-head dimensions)')
points(x = c(-0.14, -0.14, -0.14, -0.14), 
       y = c(0.03, 0.04, 0.05, 0.06), 
       col = habitat_colors, cex = 2, pch = 16)
text(x = c(-0.135, -0.135, -0.135, -0.135), 
     y = c(0.03, 0.04, 0.05, 0.06), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phylomorphospace of body shape in Pristurus', side = 3, line = 1, cex = 1, font = 2)

# Habitat broad
phymorphospace_habitat_broad <- phylomorphospace(tree = pd_habitat_broad_cons$tree[[1]], 
                                                 X = phypca_scores[, 7:8], label = 'off', 
                                                 colors = habitat_broad_colors, bty = 'n', ftype = 'off', node.by.map = TRUE,
                                                 node.size = c(0, 1.7), 
                                                 xlab = 'PC1 (-limb dimensions)', 
                                                 ylab = 'PC2 (-head dimensions)')
points(x = c(-0.14, -0.14, -0.14), 
       y = c(0.03, 0.04, 0.05), 
       col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(-0.135, -0.135, -0.135), 
     y = c(0.03, 0.04, 0.05), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
#mtext('Phylomorphospace of body shape in Pristurus', side = 3, line = 1, cex = 1, font = 2)

dev.off()

# Phylomorphospace (2 panels: insularity,  habitat broad) ----
?phylomorphospace
phylomorphospace(tree = tree, X = phypca_scores[, 7:8], label = 'off')

# Import trees mapped with discrete characters
pd_habitat_cons <- readRDS("objects/anc_rec/pd_habitat_cons.rds")
pd_habitat_broad_cons <- readRDS("objects/anc_rec/pd_habitat_broad_cons.rds")
pd_land_cons <- readRDS("objects/anc_rec/pd_land_cons.rds")

pdf(file = 'plots/phylomorphospace/phylomorphospace_2panels.pdf', width = 12, height = 5)
par(mfrow=c(1, 2))

# Mainland/island
phymorphospace_land <- phylomorphospace(tree = pd_land_cons$tree[[1]], X = phypca_scores[, 7:8], label = 'off', 
                                        colors = land_colors, bty = 'n', ftype = 'off', node.by.map = TRUE,
                                        node.size = c(0, 1.7), 
                                        xlab = 'PC1 (-limb dimensions)', 
                                        ylab = 'PC2 (-head dimensions)')
points(x = c(-0.14, -0.14), y = c(0.04, 0.05), col = land_colors, cex = 2, pch = 16)
text(x = c(-0.135, -0.135), y = c(0.04, 0.05), labels = c('mainland', 'island'), adj = 0)
#mtext('Phylomorphospace of body shape in Pristurus', side = 3, line = 1, cex = 1, font = 2)
#title(main="Phylomorphospace of body shape in Pristurus", font.main=1)

# Habitat broad
phymorphospace_habitat_broad <- phylomorphospace(tree = pd_habitat_broad_cons$tree[[1]], 
                                                 X = phypca_scores[, 7:8], label = 'off', 
                                                 colors = habitat_broad_colors, bty = 'n', ftype = 'off', node.by.map = TRUE,
                                                 node.size = c(0, 1.7), 
                                                 xlab = 'PC1 (-limb dimensions)', 
                                                 ylab = 'PC2 (-head dimensions)')
points(x = c(-0.14, -0.14, -0.14), 
       y = c(0.03, 0.04, 0.05), 
       col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(-0.135, -0.135, -0.135), 
     y = c(0.03, 0.04, 0.05), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
#mtext('Phylomorphospace of body shape in Pristurus', side = 3, line = 1, cex = 1, font = 2)

dev.off()



# Phenogram body shape ----
?phenogram
# PC1 ----
pc1 <- phypca_scores$PC1
names(pc1) <- rownames(phypca_scores)

# PC1 land
phenogram(tree = pd_land_cons$tree[[1]], x = pc1, 
          colors = land_colors, 
#          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0), y = c(-0.08, -0.10), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(-0.08, -0.10), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram PC1 Pristurus', side = 3, line = 1, cex = 1, font = 2)

# PC1 habitat
phenogram(tree = pd_habitat_cons$tree[[1]], x = pc1, 
          colors = habitat_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0, 0, 0), y = c(-0.11, -0.09, -0.07, -0.05), col = habitat_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2, 2), y = c(-0.11, -0.09, -0.07, -0.05), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC1 Pristurus', side = 3, line = 1, cex = 1, font = 2)

# PC1 habitat broad
phenogram(tree = pd_habitat_broad_cons$tree[[1]], x = pc1, 
          colors = habitat_broad_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0, 0), y = c(-0.11, -0.09, -0.07), col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2), y = c(-0.11, -0.09, -0.07), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC1 Pristurus', side = 3, line = 1, cex = 1, font = 2)



# PC2 ----
pc2 <- phypca_scores$PC2
names(pc2) <- rownames(phypca_scores)

# PC2 Land
phenogram(tree = pd_land_cons$tree[[1]], x = pc2, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0), y = c(-0.06, -0.07), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(-0.06, -0.07), labels = c('island', 'mainland'), adj = 0)
mtext('Phenogram PC2 Pristurus', side = 3, line = 1, cex = 1, font = 2)

# PC2 habitat
phenogram(tree = pd_habitat_cons$tree[[1]], x = pc2, 
          colors = habitat_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0, 0, 0), y = c(-0.06, -0.07, -0.05, -0.04), col = habitat_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2, 2), y = c(-0.06, -0.07, -0.05, -0.04), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC2 Pristurus', side = 3, line = 1, cex = 1, font = 2)

# PC2 habitat broad
phenogram(tree = pd_habitat_broad_cons$tree[[1]], x = pc2, 
          colors = habitat_broad_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0, 0), y = c(-0.06, -0.07, -0.05), col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2), y = c(-0.06, -0.07, -0.05), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC2 Pristurus', side = 3, line = 1, cex = 1, font = 2)




# Perform (non-phylogenetic) PCA WITHOUT size ----
# Remove the effect of size
colnames(morpho_log)
morpho_resid <- morpho_log %>% 
  mutate(TrL = resid(lm(TrL ~ SVL))) %>%
  mutate(HL = resid(lm(HL ~ SVL))) %>% 
  mutate(HW = resid(lm(HW ~ SVL))) %>% 
  mutate(HH = resid(lm(HH ~ SVL))) %>%
  mutate(Lhu = resid(lm(Lhu ~ SVL))) %>%
  mutate(Lun = resid(lm(Lun ~ SVL))) %>%
  mutate(Lfe = resid(lm(Lfe ~ SVL))) %>%
  mutate(Ltb = resid(lm(Ltb ~ SVL))) 

pca <- morpho_resid %>%
  select(where(is.numeric)) %>% # retain only numeric columns
  select(-SVL) %>% # remove size
  scale() %>% # scale to zero mean and unit variance
  prcomp() # do PCA

pca$rotation
# PC1: limb dimensions (Lhu, Lun, Lfe, Ltb).
# PC2: head dimensions (HL, HW, HH). 
# PC3: trunk length (TrL) and head length (HL). 
# We'll keep PC1 and PC2. 

summary(pca)

# Plot variance explained
pca %>%
  # extract eigenvalues
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) + 
  geom_col(fill = 'gray70', color = 'gray20') + 
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
morphospace_land <- pca %>%
  # add PCs to the original dataset
  augment(morpho_resid) %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = land, shape = land), size = 4) +
  scale_color_manual(values = land_colors) +
  labs(x = "PC1", y = "PC2")

morphospace_habitat <- pca %>%
  # add PCs to the original dataset
  augment(morpho_resid) %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = habitat, shape = habitat), size = 4) +
  scale_color_manual(values = habitat_colors) +
  labs(x = "PC1", y = "PC2")

morphospace_habitat_broad <- pca %>%
  # add PCs to the original dataset
  augment(morpho_resid) %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_point(aes(color = habitat_broad, shape = habitat_broad), size = 4) +
  scale_color_manual(values = habitat_broad_colors) +
  labs(x = "PC1", y = "PC2")


three_morphospace_plotlist <- list(morphospace_land, 
                                   morphospace_habitat,
                                   morphospace_habitat_broad)

three_morphospace <- cowplot::plot_grid(plotlist = three_morphospace_plotlist, 
                                        ncol = 3, 
                                        labels = 'auto')
ggsave('plots/morphospace_species/morphospace.pdf', three_morphospace, paper = 'a4r',
       height = 3, width = 15)











