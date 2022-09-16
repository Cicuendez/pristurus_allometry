# Plot phenograms of body size and shape, colored by different traits

libs <- c('tidyverse', 'treeio', 'phytools', 'geiger', 'ggtree', 'rBEAST')
lapply(libs, require, character.only = TRUE)


# Import morpho ----
morpho <- read.table('objects/phypca/phypca_scores.csv', sep = ";", dec = '.', 
                     header = TRUE, row.names = 1)
svl <- morpho$SVL
names(svl) <- rownames(morpho)
pc1 <- morpho$PC1
names(pc1) <- rownames(morpho)
pc2 <- morpho$PC2
names(pc2) <- rownames(morpho)

# Import mapped trees ----
pd_land <- readRDS('objects/anc_rec/pd_land_cons.rds')
pd_land_post <- readRDS('objects/anc_rec/pd_land_post.rds')
# pd_land_post is a list with ntrees simmap summaries of nsim simulations each.
pd_habitat_broad <- readRDS('objects/anc_rec/pd_habitat_broad_cons.rds')
pd_habitat_broad_post <- readRDS('objects/anc_rec/pd_habitat_broad_post.rds')
pd_habitat <- readRDS('objects/anc_rec/pd_habitat_cons.rds')
pd_habitat_post <- readRDS('objects/anc_rec/pd_habitat_post.rds')

# Colors
habitat_broad_colors0 <- c(ground = 'brown', rock = 'gray', tree = 'darkgreen')
habitat_broad_colors <- c(ground = '#e0710b', rock = '#6C586E', tree = '#119616')
habitat_broad_colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")
plot(1:3, 1:3, col = habitat_broad_colors, cex = 10, pch = 16)

habitat_colors0 <- c('hard-ground' = 'brown', 'soft-ground' = 'orange', 
                     rock = 'gray', tree = 'darkgreen')
habitat_colors <- c('hard-ground' = '#D63916', 'soft-ground' = '#E9A800', 
                    rock = '#6C586E', tree = '#119616')
plot(1:4, 1:4, col = habitat_colors, cex = 10, pch = 16)
land_colors0 <- c(mainland = 'coral2', island = 'seagreen')
land_colors1 <- c(mainland = '#A54E29', island = '#42978A')
land_colors2 <- c(mainland = '#C56542', island = '#67D1B7')
land_colors <- c(mainland = '#C56542', island = '#42978A')
plot(1:2, 1:2, col = land_colors, cex = 6, pch = 16)

# Phenogram body size SVL ----
# SVL land ----
pdf("plots/phenograms/phenogram_svl_land.pdf")
phenogram(tree = pd_land$tree[[1]], x = svl, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
points(x = c(0, 0), y = c(1.33, 1.3), col = NULL, bg = land_colors, cex = 2, pch = 21)
text(x = c(2, 2), y = c(1.33, 1.3), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram SVL Pristurus', side = 3, line = 1, cex = 1, font = 2)
dev.off()

# SVL habitat ----
pdf("plots/phenograms/phenogram_svl_habitat.pdf")
phenogram(tree = pd_habitat$tree[[1]], x = svl, 
          colors = habitat_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
points(x = c(0, 0, 0, 0), y = c(1.39, 1.36, 1.33, 1.3), col = habitat_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2, 2), y = c(1.39, 1.36, 1.33, 1.3), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram SVL Pristurus', side = 3, line = 1, cex = 1, font = 2)
dev.off()

# SVL habitat broad ----
par(mar = c(5,5,2,2))
pdf("plots/phenograms/phenogram_svl_habitat_broad.pdf")
phenogram(tree = pd_habitat_broad$tree[[14]], x = svl, 
          colors = habitat_broad_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
points(x = c(0, 0, 0), y = c(1.36, 1.33, 1.3), col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2), y = c(1.36, 1.33, 1.3), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram SVL Pristurus', side = 3, line = 1, cex = 1, font = 2)
dev.off()



# Phenogram body shape ----
?phenogram
# PC1 ----
# PC1 land ----
pdf("plots/phenograms/phenogram_pc1_land.pdf")
phenogram(tree = pd_land$tree[[1]], x = pc1, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0), y = c(-0.08, -0.10), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(-0.08, -0.10), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram PC1 Pristurus', side = 3, line = 1, cex = 1, font = 2)
dev.off()

# PC1 habitat ----
pdf("plots/phenograms/phenogram_pc1_habitat.pdf")
phenogram(tree = pd_habitat$tree[[1]], x = pc1, 
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
dev.off()

# PC1 habitat broad ----
pdf("plots/phenograms/phenogram_pc1_habitat_broad.pdf")
phenogram(tree = pd_habitat_broad$tree[[1]], x = pc1, 
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
dev.off()


# PC2 ----
# PC2 land ----
pdf("plots/phenograms/phenogram_pc2_land.pdf")
phenogram(tree = pd_land$tree[[1]], x = pc2, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0), y = c(-0.06, -0.07), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(-0.06, -0.07), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram PC2 Pristurus', side = 3, line = 1, cex = 1, font = 2)
dev.off()

# PC2 habitat ----
pdf("plots/phenograms/phenogram_pc2_habitat.pdf")
phenogram(tree = pd_habitat$tree[[1]], x = pc2, 
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
dev.off()

# PC2 habitat broad ----
pdf("plots/phenograms/phenogram_pc2_habitat_broad.pdf")
phenogram(tree = pd_habitat_broad$tree[[1]], x = pc2, 
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
dev.off()

# PLOT ALL PDF ----
paintBranches(tree, edge, state, anc.state)
paintBranches()
?paintBranches

pdf("plots/phenograms/all.pdf", width = 7, height = 10)
par(mfrow = c(3, 2))

# PC1 land
phenogram(tree = pd_land$tree[[1]], x = pc1, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0), y = c(-0.08, -0.10), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(-0.08, -0.10), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram PC1', side = 3, line = 1, cex = 1, font = 2)

# PC1 habitat
phenogram(tree = pd_habitat_broad$tree[[1]], x = pc1, 
          colors = habitat_broad_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0, 0), y = c(-0.11, -0.09, -0.07), col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2), y = c(-0.11, -0.09, -0.07), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC1', side = 3, line = 1, cex = 1, font = 2)


# PC2 land
phenogram(tree = pd_land$tree[[1]], x = pc2, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0), y = c(-0.06, -0.07), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(-0.06, -0.07), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram PC2', side = 3, line = 1, cex = 1, font = 2)


# PC2 habitat
phenogram(tree = pd_habitat_broad$tree[[1]], x = pc2, 
          colors = habitat_broad_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0, 0), y = c(-0.06, -0.07, -0.05), col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2), y = c(-0.06, -0.07, -0.05), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC2', side = 3, line = 1, cex = 1, font = 2)

# SVL land
phenogram(tree = pd_land$tree[[1]], x = svl, 
          colors = land_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
points(x = c(0, 0), y = c(1.33, 1.3), col = land_colors, cex = 2, pch = 16)
text(x = c(2, 2), y = c(1.33, 1.3), labels = c('mainland', 'island'), adj = 0)
mtext('Phenogram SVL', side = 3, line = 1, cex = 1, font = 2)

phenogram(tree = pd_habitat_broad$tree[[14]], x = svl, 
          colors = habitat_broad_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
points(x = c(0, 0, 0), y = c(1.36, 1.33, 1.3), col = habitat_broad_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2), y = c(1.36, 1.33, 1.3), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram SVL', side = 3, line = 1, cex = 1, font = 2)


dev.off()

# Plot 4-categories habitat svl, pc1 and pc2 ----
pdf("plots/phenograms/4habitats.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))

# svl
phenogram(tree = pd_habitat$tree[[1]], x = svl, 
          colors = habitat_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
points(x = c(0, 0, 0, 0), y = c(1.39, 1.36, 1.33, 1.3), col = habitat_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2, 2), y = c(1.39, 1.36, 1.33, 1.3), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram SVL', side = 3, line = 1, cex = 1, font = 2)

# pc1
phenogram(tree = pd_habitat$tree[[1]], x = pc1, 
          colors = habitat_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC1 (-limb dimensions)')
points(x = c(0, 0, 0, 0), y = c(-0.11, -0.09, -0.07, -0.05), col = habitat_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2, 2), y = c(-0.11, -0.09, -0.07, -0.05), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC1', side = 3, line = 1, cex = 1, font = 2)

# pc2
phenogram(tree = pd_habitat$tree[[1]], x = pc2, 
          colors = habitat_colors, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'PC2 (-head dimensions)')
points(x = c(0, 0, 0, 0), y = c(-0.06, -0.07, -0.05, -0.04), col = habitat_colors, cex = 2, pch = 16)
text(x = c(2, 2, 2, 2), y = c(-0.06, -0.07, -0.05, -0.04), 
     labels = c('hard-ground', 'soft-ground', 'rock', 'tree'), adj = 0)
mtext('Phenogram PC2', side = 3, line = 1, cex = 1, font = 2)

dev.off()








