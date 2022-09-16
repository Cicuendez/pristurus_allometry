

libs <- c('tidyverse', 'treeio', 'phytools', 'geiger', 'ggtree', 'OUwie', 
          'patchwork', 'mvMORPH')
lapply(libs, require, character.only = TRUE)


# Read morpho
morpho <- read.table('objects/phypca/phypca_scores.csv', sep = ";", dec = '.', 
                     header = TRUE)

rownames(morpho) <- morpho$species
svl <- morpho$SVL
names(svl) <- rownames(morpho)
pPC1 <- morpho$PC1
names(pPC1) <- rownames(morpho)
pPC2 <- morpho$PC2
names(pPC2) <- rownames(morpho)

svl_land_df <- morpho %>%
  dplyr::select(species, land, SVL)

# Read tree
pris <- read.nexus('data/phylogeny/pristurus_tree_final.nex')

# Establish node labels in the tree: mainland and island nodes 
# (based on the ancestral reconstruction)
nodes <- rep(1, 29)
names(nodes) <- 31:59
plot(pris)
nodelabels()
tiplabels()
pris$tip.label

island_nodes <- c('57', '58', '59', '55', '56')
nodes[island_nodes] <- 2
pris$node.label <- nodes


trait_pris <- data.frame(sp = pris$tip.label, 
                         reg = c(2, 1)[svl_land_df$land], 
                         X = svl_land_df$SVL)



?OUwie
fitted <- OUwie(pris, trait_pris, model=c("BMS"),
                algorithm = "three.point")
recon <- OUwie.anc(fitted, knowledge = TRUE)
plot(recon)

# use contMap function to plot it (to be able to change colors)
?contMap
cm_bms <- contMap(pris, svl, method = 'user', anc.states = recon$NodeRecon)

par(mar=c(5,5,2,2))
col_scale <- rev(brewer.pal(10, "RdYlBu"))
plot(c(1:10), c(1:10), col=col_scale, pch=16, cex=3)
palette <- colorRampPalette(col_scale)
cm_bms$cols[] <- palette(1001)

plot(cm_bms, outline = F, leg.txt = "log(SVL)")


# svl (reconstructed) for nodes
r <- recon$NodeRecon
names(r) <- c(31:59)

# svl (empirical) for tips
m <- svl_land_df$SVL[match(pris$tip.label, svl_land_df$species)]
names(m) <- c(1:30)

# svl for tips and nodes
valSVL <- c(m,r)

###
parent_node_bodySize <- valSVL[match(pris$edge[,1],names(valSVL))]
daughter_node_bodySize <- valSVL[match(pris$edge[,2],names(valSVL))]
branch_dist <- abs(daughter_node_bodySize - parent_node_bodySize)
branch_dist_noabs <- daughter_node_bodySize - parent_node_bodySize
names(branch_dist) <- c(1:58)

palespectral <- colorRampPalette(c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", 
                                   "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD"))
spectral_rev <- rev(c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", 
                          "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD"))
palspectral_rev <- colorRampPalette(spectral_rev)

branch_rate_tree <- plotBranchbyTrait(pris, branch_dist, cex=0.8, 
                                      palette = palspectral_rev)
edgelabels()
plot(pris)
edgelabels()
nodelabels()

# Branches with transition mainland to island:
# 51, 46, 7
branches_w_change <- c(7, 46, 51)

max(branch_dist)
branch_dist[7]

# paired t test all ----
paired_test <- t.test(parent_node_bodySize, daughter_node_bodySize, 
       paired = TRUE, alternative = 'two.sided')
summary(paired_test)
hist(parent_node_bodySize, col = 'gray90')
hist(daughter_node_bodySize, col = 'red', add = TRUE)



# PARENT AND DESCENDENT DATAFRAME ----
parent_descend <- data.frame(parent = parent_node_bodySize, descend = daughter_node_bodySize)
rownames(parent_descend) <- 1:58


par(mar=c(5,5,2,2))
boxplot(parent_descend)
?boxplot

# Plot paired data
install.packages('PairedData')
library(PairedData)
p <- paired(before, after)
plot(pd, type = "profile") + theme_bw()

island_change <- rep(1, 58)
names(island_change) <- 1:58
# Branches with transition mainland to island:
# 51, 46, 7
branches_w_change <- c(7, 46, 51)
island_change[branches_w_change] <- 2
parent_descend <- cbind(parent_descend, island_change)
parent_descend$island_change <- as.factor(parent_descend$island_change)
parent_descend <- cbind(parent_descend, branch_dist)


library(ggpubr)
ggpaired(parent_descend, cond1 = 'parent', cond2 = 'descend', 
         color = 'island_change', line.color = 'island_change', 
         palette = c('gray80', 'darkred'), line.size = 0.2, 
         xlab = 'node type', ylab = 'log body size')

boxplot_ggpaired <- parent_descend %>%
  filter(island_change == 2) %>%
  ggpaired(cond1 = 'parent', cond2 = 'descend', 
 #          color = 'island_change', 
 #          line.color = 'island_change', 
 #          palette = c('gray80', 'darkred'), line.size = 0.2, 
           xlab = 'node type', ylab = 'log body size')

# ANOVA ----
# es la distancia entre nodos distinta cuando se pasa de mainland a island?
# (comparado con mainland-mainland o island-island)
summary(aov(branch_dist ~ island_change, data = parent_descend))

# Histogram of branch dist (SVL difference between parent and descendent nodes)
hist(branch_dist)
abline(v = branch_dist[branches_w_change], col = 'red')

hist(branch_dist_noabs)
shapiro.test(branch_dist_noabs)
qqnorm(branch_dist_noabs)
?shapiro.test

# paired t test NOT MAINLAND ISLAND ----
parent_descend[-branches_w_change,]
t.test(parent_descend[-branches_w_change,]$parent, 
       parent_descend[-branches_w_change,]$descend, 
       paired = TRUE, alternative = 'two.sided')

# paired t test ONLY MAINLAND ISLAND ----
parent_descend[branches_w_change,]
t.test(parent_descend[branches_w_change,]$parent, 
       parent_descend[branches_w_change,]$descend, 
       paired = TRUE, alternative = 'two.sided')


#pgls
identical(as.character(svl_land_df$species), pris$tip.label)
phyl.anova_insularity <- phylANOVA(tree = pris, x = svl_land_df$land, y = svl_land_df$SVL)
phyl.anova_insularity$Pf

boxplot(svl_land_df$SVL ~ svl_land_df$land)

identical(as.character(morpho$species), pris$tip.label)
phyl.anova_habitat <- phylANOVA(tree = pris, x = morpho$habitat_broad, 
                        y = morpho$SVL)
boxplot(morpho$SVL ~ morpho$habitat)

# Plot contmap and paired boxplot ----
par(mfrow = c(1, 1))
pdf('plots/contmapSVL_BMS.pdf')
plot(cm_bms, outline = F, leg.txt = "log(SVL)")
dev.off()

pdf('plots/boxplot_ggpaired.pdf')
boxplot_ggpaired
dev.off()



