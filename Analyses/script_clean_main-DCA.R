# R-script for primary analyses

libs <- c('RRPP', 'phytools', 'geiger', 'tidyverse')
easypackages::libraries(libs)

data0 <- read.table('data/morpho/morpho_pristurus.csv', sep = ';', dec = '.', header = TRUE)
  sp.to.keep <- names(which(table(data0$species) >= 5))
data <- data0[data0$species %in% sp.to.keep, ]
  data$species <- as.factor(data$species)
  data$habitat_broad <- as.factor(data$habitat_broad)  
  data$SVL <- log(data$SVL)
shape <- as.matrix(log(data[, 8:ncol(data)]))
rdf <- rrpp.data.frame(svl = data$SVL, shape = shape, habitat = data$habitat_broad)
tree0 <- read.nexus('data/phylogeny/pristurus_tree_final.nex')

# Comparison of multivariate allometry among habitat types
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
  anova(fit.hab)
pw.hab <- pairwise(fit.hab, groups = rdf$habitat, covariate = rdf$svl)
  summary(pw.hab, type = 'VC', stat.table = FALSE)

#Slopes by habitat
fit.coef <- fit.hab$LM$coefficients
ind.coef <- rbind(fit.coef[2,], fit.coef[2,]+fit.coef[5,], fit.coef[2,]+fit.coef[6,])
rownames(ind.coef) <- c("Ground","Rock", "Tree")

  # Analysis using species means
LS.mns <- pairwise(lm.rrpp(shape~species, data = rdf, iter=0), groups = rdf$species)$LS.means[[1]]
sz.mn <- tapply(rdf$svl,rdf$species,mean)
hab.mn <- as.factor(by(rdf$habitat,rdf$species,unique))
 levels(hab.mn) <- levels(rdf$habitat)
tree <- treedata(phy = tree0, data = LS.mns)$phy
C <- vcv.phylo(tree)

#Similar anova
fit.hab.sp <- lm.rrpp(LS.mns~sz.mn*hab.mn)
anova(fit.hab.sp)

#pw.hab.sp <- pairwise(fit.hab.sp,groups = hab.mn, covariate = sz.mn)
#summary(pw.hab.sp, type = 'VC', stat.table = FALSE) 
  #similar pattern: ground vs. rock ***, rock vs. tree NS (also ground vs. tree NS)

### NOTE: when accounting for phylogeny, evol. allometry still present, but not the habitat effect
fit.hab.sp.phy <- lm.rrpp(LS.mns~sz.mn*hab.mn, Cov = C)
anova(fit.hab.sp.phy)


fit.coef.sp <- fit.hab.sp$LM$coefficients
sp.coef <- rbind(fit.coef.sp[2,], fit.coef.sp[2,]+fit.coef.sp[5,], fit.coef.sp[2,]+fit.coef.sp[6,])
rownames(sp.coef) <- c("Ground","Rock", "Tree")

ind.coef
sp.coef
plot(fit.hab.sp, type = "regression",predictor = sz.mn)
##########################################

# get species habitats
# Drop species with less than 5 individuals from the tree
data.sp <- data.sp %>% 
  filter(species %in% sp.to.keep)
rownames(data.sp) <- data.sp$species

hab.sp <- data.sp$habitat_broad
names(hab.sp) <- data.sp$species

# size model (species data) ----
# Prepare data 
# Drop species with less than 5 individuals from the tree
data.sp <- data.sp %>% 
  filter(species %in% sp.to.keep)


tree0 <- read.nexus('data/phylogeny/pristurus_tree_final.nex')


dat.tree.morpho <- treedata(phy = tree0, data = data.sp)
tree <- dat.tree.morpho$phy
dat.morpho <- as.data.frame(dat.tree.morpho$data)
class(dat.morpho$SVL)
# class is character, let's do it numeric
dat.num <- dat.morpho[, c(6:ncol(dat.morpho))]
for (i in 1:ncol(dat.num)){
  dat.num[,i] <- as.numeric(as.character(dat.num[,i]))
}
dat.morpho.num <- cbind(dat.morpho[, 1:5], dat.num)
svl.sp <- log(dat.morpho.num$SVL)
names(svl.sp) <- dat.morpho.num$species

# shape~size phylogenetic regression ----
shape.sp <- as.matrix(log(dat.morpho.num[, 7:ncol(dat.morpho.num)]))
svl.sp
identical(names(svl.sp), rownames(shape.sp))

# Prepare data
rdf.size <- rrpp.data.frame(svl = svl.sp, shape = shape.sp)

# Multivariate linear model
fit.svl <- lm.rrpp(shape~svl, data = rdf.size, Cov = vcv(tree))
anova(fit.svl)

shape.res <- residuals(fit.svl)


# PHYLOMORPHOSPACE ggplot2 ----
# PLOT shape residuals PCA.OLS with ggplot2
# PCA shape residuals
# ordinary PCA
pca.shape.ols <- gm.prcomp(shape.res, phy = tree)
x <- rownames_to_column(as.data.frame(pca.shape.ols$x), 
                        var = 'species')
svl.sp
identical(x$species, names(svl.sp))

x$size <- svl.sp

dat.ggphylo <- x %>% left_join(data.sp) %>% 
  mutate(PC1 = Comp1) %>% 
  mutate(PC2 = Comp2) %>% 
  mutate(taxon = species) %>%
  mutate(group = habitat_broad) %>%
  select(taxon, PC1, PC2, group, size)

rownames(dat.ggphylo) <- dat.ggphylo$taxon

# We slightly modified the function ggphylomorpho from the package with 
# the same name, to allow for annotating only certain species and to 
# additionally for changing the point size based on a certain variable (this 
# latter option in the function 'ggphylomorpho_size_HTC')
ggphylomorpho_HTC <- function (tree, tipinfo, 
                               taxa = taxon, 
                               xvar = PC1, 
                               yvar = PC2, 
                               factorvar = group, 
                               labelvar = splabel, 
                               title.plot = "Phylomorphospace", 
                               label.x.axis = "PC1", 
                               label.y.axis = "PC2", 
                               repel = TRUE, 
                               edge.width = 1, 
                               fontface = "italic", 
                               tree.alpha = 0.7)
{
  require(ggplot2)
  require(phytools)
  require(ggrepel)
  mat <- cbind(eval(substitute(xvar), tipinfo), 
               eval(substitute(yvar), tipinfo))
  rownames(mat) <- eval(substitute(taxa), tipinfo)
  stopifnot(length(setdiff(tree$tip.label, rownames(mat))) == 0)
  
  xAnc <- fastAnc(tree, mat[, 1])
  yAnc <- fastAnc(tree, mat[, 2])
  all_node_coords <- data.frame(x = c(mat[tree$tip.label, 1], xAnc), 
                                y = c(mat[tree$tip.label, 2], yAnc), 
                                nodeid = 1:(tree$Nnode + length(tree$tip.label)))
  edges <- data.frame(tree$edge)
  names(edges) <- c("node1", "node2")
  
  edgecoords <- merge(merge(edges, all_node_coords, by.x = "node1", 
                            by.y = "nodeid"), 
                      all_node_coords, by.x = "node2", by.y = "nodeid")
  
  pointsForPlot <- data.frame(x = eval(substitute(xvar), tipinfo), 
                              y = eval(substitute(yvar), tipinfo), 
                              color = eval(substitute(factorvar), tipinfo), 
                              label = eval(substitute(labelvar), tipinfo))
  
  
  theplot <- ggplot() + 
    geom_segment(data = edgecoords, 
                 aes(x = x.x, xend = x.y, y = y.x, yend = y.y), 
                 size = edge.width, 
                 alpha = tree.alpha) + 
    geom_point(data = pointsForPlot, 
               aes(x = x, y = y, color = color), 
               size = 5) + 
    labs(title = title.plot, x = label.x.axis, y = label.y.axis) + 
    theme_bw(10) + 
    theme(legend.position = "bottom", 
          plot.title = element_text(size = 15))
  
  if (repel) {
    theplot <- theplot + 
      geom_text_repel(data = pointsForPlot, 
                      aes(x = x, y = y, label = label), 
                      segment.alpha = 0.5, 
                      fontface = fontface)
  }
  else {
    theplot <- theplot + 
      geom_text(data = pointsForPlot, 
                aes(x = x, y = y, label = label), 
                fontface = fontface)
  }
  return(theplot)
}

ggphylomorpho_size_HTC <- function (tree, tipinfo, 
                                    taxa = taxon, 
                                    xvar = PC1, 
                                    yvar = PC2, 
                                    factorvar = group, 
                                    labelvar = splabel, 
                                    point.size = size, 
                                    title.plot = "Phylomorphospace", 
                                    label.x.axis = "PC1", 
                                    label.y.axis = "PC2", 
                                    repel = TRUE, 
                                    edge.width = 1, 
                                    fontface = "italic", 
                                    tree.alpha = 0.7) 
{
  require(ggplot2)
  require(phytools)
  require(ggrepel)
  mat <- cbind(eval(substitute(xvar), tipinfo), 
               eval(substitute(yvar), tipinfo))
  rownames(mat) <- eval(substitute(taxa), tipinfo)
  stopifnot(length(setdiff(tree$tip.label, rownames(mat))) == 0)
  
  xAnc <- fastAnc(tree, mat[, 1])
  yAnc <- fastAnc(tree, mat[, 2])
  all_node_coords <- data.frame(x = c(mat[tree$tip.label, 1], xAnc), 
                                y = c(mat[tree$tip.label, 2], yAnc), 
                                nodeid = 1:(tree$Nnode + length(tree$tip.label)))
  edges <- data.frame(tree$edge)
  names(edges) <- c("node1", "node2")
  
  edgecoords <- merge(merge(edges, all_node_coords, by.x = "node1", 
                            by.y = "nodeid"), 
                      all_node_coords, by.x = "node2", by.y = "nodeid")
  
  pointsForPlot <- data.frame(x = eval(substitute(xvar), tipinfo), 
                              y = eval(substitute(yvar), tipinfo), 
                              color = eval(substitute(factorvar), tipinfo), 
                              label = eval(substitute(labelvar), tipinfo), 
                              size = eval(substitute(size), tipinfo))
  
  
  theplot <- ggplot() + 
    geom_segment(data = edgecoords, 
                 aes(x = x.x, xend = x.y, y = y.x, yend = y.y), 
                 size = edge.width, 
                 alpha = tree.alpha) + 
    geom_point(data = pointsForPlot, 
               aes(x = x, y = y, color = color, size = size)) + 
    labs(title = title.plot, x = label.x.axis, y = label.y.axis) + 
    theme_bw(5) + 
    theme(legend.position = "bottom")
  
  if (repel){
    theplot <- theplot + 
      geom_text_repel(data = pointsForPlot, 
                      aes(x = x, y = y, label = label), 
                      segment.alpha = 0.5, 
                      fontface = fontface)
  }
  else {
    theplot <- theplot + 
      geom_text(data = pointsForPlot, 
                aes(x = x, y = y, label = label), 
                fontface = fontface)
  }
  return(theplot)
}

hab.colors
phymor.plot <- ggphylomorpho_HTC(tree = tree, tipinfo = dat.ggphylo, 
                             labelvar = taxon, edge.width = 0.3)
phymor.plot + 
  scale_color_manual(values = hab.colors) +
  coord_fixed() +
  theme.clean() + 
  theme(legend.position = 'bottom')

# PLOT PCA.OLS large species ----
large.sp.rock <- c('Pristurus_insignis', 'Pristurus_insignoides')
large.sp.ground <- c('Pristurus_carteri', 'Pristurus_ornithocephalus', 
                     'Pristurus_collaris')
large.sp <- c(large.sp.rock, large.sp.ground)

'%nin%' <- Negate('%in%')
dat.ggphylo.large <- dat.ggphylo %>% 
  mutate(group_large = case_when(group == 'rock' & taxon %in% large.sp.rock ~ 'rock',
                                 group == 'rock' & taxon %nin% large.sp.rock ~ 'rock_faded',
                                 group == 'tree' ~ 'tree_faded', 
                                 group == 'ground' & taxon %in% large.sp.ground ~ 'ground', 
                                 group == 'ground' & taxon %nin% large.sp.ground ~ 'ground_faded')) %>% 
  mutate(splabel = case_when(taxon %in% large.sp ~ taxon, TRUE ~ '')) %>% 
  mutate(splabel = gsub('Pristurus_', 'P. ', .$splabel))

# define new colors by adding transparency (to highlight only certain species)
hab.colors.faded <- c(ground = "#F1B670", rock = "#683B5E", 
                      tree = "#E93F7B", 
                      tree_faded = '#E93F7B40', 
                      ground_faded = '#F1B67040', 
                      rock_faded = '#683B5E40')
dat.ggphylo.large$splabel
phymor.plot.large <- ggphylomorpho_HTC(tree = tree, 
                                       tipinfo = dat.ggphylo.large, 
                                       taxa = taxon, 
                                       xvar = PC1, 
                                       yvar = PC2, 
                                       factorvar = group_large, 
                                       labelvar = splabel, 
                                       title.plot = "Phylomorphospace", 
                                       label.x.axis = "PC1", 
                                       label.y.axis = "PC2", 
                                       repel = TRUE, 
                                       edge.width = 0.3, 
                                       fontface = "italic", 
                                       tree.alpha = 0.7)

large.phylomorpho.plot <- phymor.plot.large + 
  scale_color_manual(values = hab.colors.faded) +
  coord_fixed() +
  theme.clean() + 
  theme(legend.position = 'bottom')



# PLOT PCA.OLS small species ----
small.sp.rock <- c('Pristurus_sp5', 'Pristurus_sp1', 'Pristurus_rupestris', 
                   'Pristurus_sp3', 'Pristurus_sp2')
small.sp.ground <- c('Pristurus_masirahensis', 'Pristurus_minimus')
small.sp <- c(small.sp.rock, small.sp.ground)

dat.ggphylo.small <- dat.ggphylo %>% 
  mutate(group_small = case_when(group == 'rock' & taxon %in% small.sp.rock ~ 'rock',
                                 group == 'rock' & taxon %nin% small.sp.rock ~ 'rock_faded',
                                 group == 'tree' ~ 'tree_faded', 
                                 group == 'ground' & taxon %in% small.sp.ground ~ 'ground', 
                                 group == 'ground' & taxon %nin% small.sp.ground ~ 'ground_faded')) %>% 
  mutate(splabel = case_when(taxon %in% small.sp ~ taxon, TRUE ~ '')) %>% 
  mutate(splabel = gsub('Pristurus_', 'P. ', .$splabel))

phymor.plot.small <- ggphylomorpho_HTC(tree = tree, 
                                       tipinfo = dat.ggphylo.small, 
                                       taxa = taxon, 
                                       xvar = PC1, 
                                       yvar = PC2, 
                                       factorvar = group_small, 
                                       labelvar = splabel, 
                                       title.plot = "Phylomorphospace", 
                                       label.x.axis = "PC1", 
                                       label.y.axis = "PC2", 
                                       repel = TRUE, 
                                       edge.width = 0.3, 
                                       fontface = "italic", 
                                       tree.alpha = 0.7)

small.phylomorpho.plot <- phymor.plot.small + 
  scale_color_manual(values = hab.colors.faded) +
  coord_fixed() +
  theme.clean() + 
  theme(legend.position = 'bottom')

library(patchwork)
large.phylomorpho.plot | small.phylomorpho.plot

# PLOT PCA.OLS large & small species ----
# With body size mapped to point size
large.sp.rock 
large.sp.ground 
large.sp

small.sp.rock
small.sp.ground 
small.sp

large_small.rock <- c(large.sp.rock, small.sp.rock)
large_small.ground <- c(large.sp.ground, small.sp.ground)
large_small.sp <- c(large.sp, small.sp)

dat.ggphylo.large_small <- dat.ggphylo %>% 
  mutate(group_ls = case_when(group == 'rock' & taxon %in% large_small.sp ~ 'rock',
                                 group == 'rock' & taxon %nin% large_small.sp ~ 'rock_faded',
                                 group == 'tree' ~ 'tree_faded', 
                                 group == 'ground' & taxon %in% large_small.sp ~ 'ground', 
                                 group == 'ground' & taxon %nin% large_small.sp ~ 'ground_faded')) %>% 
  mutate(splabel = case_when(taxon %in% large_small.sp ~ taxon, TRUE ~ '')) %>% 
  mutate(splabel = gsub('Pristurus_', 'P. ', .$splabel))

phymor.plot.ls <- ggphylomorpho_size_HTC(tree = tree, 
                                       tipinfo = dat.ggphylo.large_small, 
                                       taxa = taxon, 
                                       xvar = PC1, 
                                       yvar = PC2, 
                                       point.size = size,
                                       factorvar = group_ls, 
                                       labelvar = splabel, 
                                       title.plot = "Phylomorphospace", 
                                       label.x.axis = "PC1", 
                                       label.y.axis = "PC2", 
                                       repel = TRUE, 
                                       edge.width = 0.1, 
                                       fontface = "italic", 
                                       tree.alpha = 0.7)

large_small.phylomorpho.plot <- phymor.plot.ls + 
  scale_color_manual(values = hab.colors.faded, 
                     name = 'habitat', 
                     breaks = c('ground', 'rock', 'tree')) +
  guides(size = 'none') +
  coord_fixed() +
  labs(title = 'Phylomorphospace', 
       subtitle = '(largest and smallest species highlighted)') +
  theme.clean() + 
  theme(legend.position = 'bottom')

ggsave('plots/phylomorphospace_large_small.pdf', 
       plot = large_small.phylomorpho.plot)
ggsave('plots/phylomorphospace_large_small.png', 
       plot = large_small.phylomorpho.plot)

?scale_color_manual



# Plot large and small with base R ----
dat.ggphylo.large$group_large
hab.colors.large <- c(ground = "#F1B670", rock = "#683B5E", 
                      tree = "#E93F7B30", ground_faded = '#F1B67030', 
                      rock_faded = '#683B5E30')
hab.outline.large <- setNames(c('black', 'black', 'gray80', 
                                'gray80', 'gray80'), 
                              names(hab.colors.large))
plot.options <- list(node.labels = FALSE, anc.states = FALSE, 
                     edge.color = 'gray40', edge.width = 0.5, 
                     tip.labels = FALSE)
plot(pca.shape.ols, phylo = TRUE, pch = 21, phylo.par = plot.options, 
     bg = hab.colors.faded[dat.ggphylo.large$group_large], 
     col = hab.colors.faded[dat.ggphylo.large$group_large], 
     include.axes = FALSE, 
     cex = svl.sp-2)




# shape ~ habitat model (species data) ----
shape.res
hab.sp <- dat.morpho.num$habitat_broad
names(hab.sp) <- rownames(dat.morpho.num)
identical(names(hab.sp), rownames(shape.res))
hab.sp.fctr <- as.factor(hab.sp)

rdf.shape.hab <- rrpp.data.frame(shape = shape.res, hab = hab.sp.fctr)

# Phylogenetic regression
fit.shape.hab <- lm.rrpp(shape~hab, data = rdf.shape.hab, Cov = vcv(tree))
anova(fit.shape.hab)

# Pairwise differences in shape variance ----
pw.shape.hab <- pairwise(fit.shape.hab, groups = rdf.shape.hab$hab)
pw.shape.hab_df <- summary(pw.shape.hab, type = 'var', stat.table = FALSE)
?lm.rrpp

# head ~ habitat model ----
hab.sp.fctr
head.res <- shape.res[, 2:4]
rdf.head.hab <- rrpp.data.frame(head = head.res, hab = hab.sp.fctr)
fit.head.hab <- lm.rrpp(head~hab, data = rdf.head.hab, Cov = vcv(tree))
anova(fit.head.hab)
pw.head.hab <- pairwise(fit.head.hab, groups = rdf.head.hab$hab)
pw.head.hab_df <- summary(pw.head.hab, type = 'var', stat.table = FALSE)

# limb ~ habitat model ----
limb.res <- shape.res[, 5:8]
rdf.limb.hab <- rrpp.data.frame(limb = limb.res, hab = hab.sp.fctr)
fit.limb.hab <- lm.rrpp(limb~hab, data = rdf.limb.hab, Cov = vcv(tree))
anova(fit.limb.hab)
pw.limb.hab <- pairwise(fit.limb.hab, groups = rdf.limb.hab$hab)
pw.limb.hab_df <- summary(pw.limb.hab, type = 'var', stat.table = FALSE)

# PLS ----
# PLS head vs size 
svl <- log(data$SVL)
shape.head <- shape[, c(2:4)]

pls.head <- two.b.pls(shape.head, svl)

# PLS limbs vs size
shape.limb <- shape[, 5:8]
pls.limb <- two.b.pls(shape.limb, svl)

# Get scores
limb.scores <- pls.limb$XScores[, 1]
head.scores <- pls.head$XScores[, 1]

par(mfrow=c(1,2))
plot(svl, pls.limb$XScores[, 1], pch = 22, bg = hab.colors[habitat.fctr])
plot(svl, pls.head$XScores[, 1], pch = 25, 
       bg = hab.colors[habitat.fctr])

# LM scores ~ svl*sp 
# Get a slope per species for head and limbs 
# (instead of for each variable as we did before)

fit.head.sp <- lm.rrpp(head.scores ~ svl*species.fctr)
fit.limb.sp <- lm.rrpp(limb.scores ~ svl*species.fctr)
anova(fit.head.sp)
anova(fit.limb.sp)

plot(fit.head.sp, predictor = svl, type = 'regression', pch = 16, 
     col = species.fctr)
plot(fit.limb.sp, predictor = svl, type = 'regression', pch = 16, 
     col = species.fctr)

# regression coefficients
coef.head <- fit.head.sp$LM$coefficients
coef.limb <- fit.limb.sp$LM$coefficients

# get per-species slopes
head.slp <- coef.head[grep('svl', rownames(coef.head)), ]
head.slp[-1] <- head.slp[-1] + head.slp[1]
names(head.slp) <- gsub('svl:species.fctr', '', names(head.slp))
names(head.slp) <- gsub('svl', 'Pristurus_abdelkuri', names(head.slp))

limb.slp <- coef.limb[grep('svl', rownames(coef.limb)), ]
limb.slp[-1] <- limb.slp[-1] + limb.slp[1]
names(limb.slp) <- gsub('svl:species.fctr', '', names(limb.slp))
names(limb.slp) <- gsub('svl', 'Pristurus_abdelkuri', names(limb.slp))

# We already have the tree and the data matched 
identical(tree$tip.label, names(head.slp))

# contMap head and limb slopes ----
ramp <- colorRampPalette(c("#00929c","gray80",  "#d62e31"))
cm.head <- contMap(tree = tree, x = head.slp, outline = FALSE)
cm.head$cols[] <- ramp(1001)
pdf('plots/cm_head.pdf')
plot(cm.head, outline = FALSE)
dev.off()

cm.limb <- contMap(tree = tree, x = limb.slp, outline = FALSE)
cm.limb$cols[] <- ramp(1001)
pdf('plots/cm_limb.pdf')
plot(cm.limb, outline = FALSE)
dev.off()

# phenograms head and limb slopes ----
pdf('plots/phenogram_head.pdf')
phenogram(tree = cm.head$tree, x = svl.sp, 
          colors = cm.head$cols, 
          #          ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
add.color.bar(prompt = FALSE, cols = cm.head$cols, leg = 20, 
              title = 'head slope', 
              subtitle = '',
              lwd = 8, lims = cm.head$lims, 
              outline = FALSE, 
              x = 5, y = 3, fsize = 0.8)
dev.off()

pdf('plots/phenogram_limb_reverse.pdf')
phenogram(tree = cm.limb$tree, x = svl.sp, 
          colors = cm.limb$cols, 
          #ftype = 'off', 
          ftype = 'i',
          fsize = 0.7, 
          spread.cost=c(1,0),
          xlab = 'time (ma)', ylab = 'logSVL')
add.color.bar(prompt = FALSE, cols = cm.limb$cols, leg = 20, 
              title = 'limb slope', 
              subtitle = '',
              lwd = 8, lims = cm.limb$lims, 
              outline = FALSE, 
              x = 5, y = 3, fsize = 0.8)
dev.off()

# These two plots, the phenograms of the head and the limbs slopes, 
# have been combined in Illustrator for the final version of the figure.




######### THE END ############################




