

# Packages ----
libs <- c('geomorph', 'RRPP', 'phytools', 'geiger', 'tidyverse', 
          'ggphylomorpho')
easypackages::libraries(libs)

# Morpho data ----
# Both species data and specimen data
data.sp <- read.table('data/morpho/morpho_sp_final.csv', sep = ';', 
                      dec = '.', header = TRUE)
data0 <- read.table('data/morpho/morpho_pristurus.csv', sep = ';', 
                    dec = '.', header = TRUE)

# Drop species with less than 5 individuals
sp.to.keep <- names(which(table(data0$species) > 5) == TRUE)
data <- data0[data0$species %in% sp.to.keep, ]

# number of species
n.sp <- length(table(data$species))

# Phylogeny ----
tree0 <- read.nexus('data/phylogeny/pristurus_tree_final.nex')

# Prepare specimen data for analysis ----
svl <- log(data$SVL)
shape <- as.matrix(log(data[, 8:ncol(data)]))
species.fctr <- as.factor(data$species)
habitat.fctr <- as.factor(data$habitat_broad)

rdf <- rrpp.data.frame(svl = svl, shape = shape, habitat = habitat.fctr, 
                       species = species.fctr)

# Habitat model ----
# Get slopes per habitat
# Multivariate linear model
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
anova(fit.hab)
##### regression coefficients
fit.coef <- fit.hab$LM$coefficients

rbind(fit.coef[1,], fit.coef[2,]) #ground
rbind(fit.coef[1,]+fit.coef[3,], fit.coef[2,]+fit.coef[5,]) #rock
rbind(fit.coef[1,]+fit.coef[4,], fit.coef[2,]+fit.coef[6,]) #tree 

# Pairwise differences in the angle ----
pw.hab <- pairwise(fit.hab, groups = rdf$habitat, covariate = rdf$svl)
pw.hab_df <- summary(pw.hab, type = 'VC', stat.table = FALSE)

pw.hab_df

# Set habitat colors ----
hab.colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")

# Set customized ggplot2 theme ----
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
          plot.title = element_text(size = 15, vjust = 1, hjust = 0.5, 
                                    face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5), 
          legend.text = element_text(size = 10, face = "plain"),          
          legend.position = "right")
}


# plot base ----
# plot predicted values
plot.base <- plot(fit.hab, predictor = rdf$svl, type = 'regression', pch = 16, 
     col = hab.colors[as.numeric(rdf$habitat)])
legend('topleft', levels(habitat.fctr), pt.bg = hab.colors, pch = 22)

plot(rdf$svl, plot.base$PredLine, col = hab.colors[as.numeric(rdf$habitat)], 
     pch = 16) # same thing
# This is like taking the fitted values of the model and doing a PCA; 
# the predicted line we see is the PC1.

# plot RegScores
plot(fit.hab, predictor = rdf$svl, type = 'regression', reg.type = 'RegScore', 
     pch = 16, 
     col = hab.colors[as.numeric(rdf$habitat)])
legend('topleft', levels(habitat.fctr), pt.bg = hab.colors, pch = 22)

plot(rdf$svl, plot.base$RegScore, col = hab.colors[as.numeric(rdf$habitat)], 
     pch = 16) # same thing

# plot ggplot2 ----
# Figure 2
# Plot habitat slopes with ggplot2 
fit.hab.ggplot.data <- data.frame(fitted_PC1 = plot.base$PredLine,
                                  RegScore = plot.base$RegScore[,1],
                                  svl = rdf$svl, 
                                  habitat = rdf$habitat)

habitat.slope.plot <- ggplot(data = fit.hab.ggplot.data, aes(x = svl)) +
  geom_point(aes(y = RegScore, color = habitat), alpha = 0.2, size = 2) +
  geom_line(aes(y = fitted_PC1, color = habitat), size = 1) +
  scale_color_manual(values = hab.colors) +
  labs(x = 'logSVL', y = 'Regression Scores') +
  theme.clean() +
  theme(legend.position = 'bottom')

ggsave('plots/figure_2_ggplot.png', habitat.slope.plot)



# Species model ----
# Get slopes per species
# Multivariate linear model
fit.sp <- lm.rrpp(shape~svl*species, data = rdf)
anova(fit.sp)

plot(fit.sp, predictor = rdf$svl, type = 'regression', pch = 16, 
     col = as.numeric(rdf$species))

# regression coefficients
fit.coef.sp <- fit.sp$LM$coefficients

# get species slopes ----
sp.slp <- fit.coef.sp[grep('svl', rownames(fit.coef.sp)), ]
sp.slp[-1,] <- sp.slp[-1,] + sp.slp[1,]

rownames(sp.slp) <- gsub('svl:species', '', rownames(sp.slp))
rownames(sp.slp) <- gsub('svl', 'Pristurus_abdelkuri', rownames(sp.slp))

# Plot species slopes on the phylogeny ----
# match the species in the tree and in the data
dat.tree.slp <- treedata(phy = tree0, data = sp.slp)
tree <- dat.tree.slp$phy
dat.slp <- dat.tree.slp$data

# contMap slopes per variable ----
ramp <- colorRampPalette(c("#00929c","gray80",  "#d62e31"))
cm.plots <- vector('list', length = ncol(dat.slp))
names(cm.plots) <- colnames(dat.slp)

for (i in 1:ncol(dat.slp)){
  print(names(cm.plots)[i])
  cm <- contMap(tree = tree, x = dat.slp[,i], outline = FALSE)
  cm$cols[] <- ramp(1001)
  cm.plots[[i]] <- cm
  plot.filename <- paste0('plots/cm_', names(cm.plots[i]), '.png')
  png(plot.filename)
  plot(cm, outline = FALSE)
  dev.off()
}

# contMap SVL ----
rownames(data.sp) <- data.sp$species
dat.tree.morpho <- treedata(phy = tree, data = data.sp)
dat.morpho <- as.data.frame(dat.tree.morpho$data)
class(dat.morpho$SVL)
dat.num <- dat.morpho[, c(6:ncol(dat.morpho))]
for (i in 1:ncol(dat.num)){
  dat.num[,i] <- as.numeric(as.character(dat.num[,i]))
}
dat.morpho.num <- cbind(dat.morpho[, 1:5], dat.num)
svl.sp <- log(dat.morpho.num$SVL)
names(svl.sp) <- dat.morpho.num$species

cm.svl <- contMap(tree = tree, x = svl.sp, outline = FALSE, plot = FALSE)
cm.svl$cols[] <- ramp(1001)

# Phenograms ----
# phenogram slopes and SVL
# __ Y = SVL, color = slope ----

# separate files
for (i in 1:length(cm.plots)){
  phname <- paste0('plots/phenogram_svl_', names(cm.plots)[i], '.png')
  png(phname)
  phenogram(tree = cm.plots[[i]]$tree, x = svl.sp, 
            colors = cm.plots[[i]]$cols, 
            #          ftype = 'off', 
            ftype = 'i',
            fsize = 0.7, 
            spread.cost=c(1,0),
            xlab = 'time (ma)', ylab = 'logSVL')
  points(x = 30.5, y = 3.65, bg = 'orchid3', col = 'black', 
         cex = 1.5, pch = 21)
  text(x = 31, y = 3.6, cex = 0.8, col = 'orchid3',
       labels = substitute(paste(bold('hard ground'))),
       adj = 0)
  add.color.bar(prompt = FALSE, cols = cm.plots[[i]]$cols, leg = 20, 
                title = paste0(names(cm.plots[i]), ' slope'), 
                subtitle = '',
                lwd = 8, lims = cm.plots[[i]]$lims, 
                outline = FALSE, 
                x = 5, y = 3, fsize = 0.8)
  dev.off()
}

# all in one file
png('plots/phenograms_all.png', width = 1200, height = 600)
par(mfrow = c(2, 4))
par(mar = c(5,5,2,2))
for (i in 1:length(cm.plots)){
  phenogram(tree = cm.plots[[i]]$tree, x = svl.sp, 
            colors = cm.plots[[i]]$cols, 
            ftype = 'off', 
            #ftype = 'i',
            #fsize = 0.7, 
            #spread.cost=c(1,0),
            xlab = 'time (ma)', ylab = 'logSVL')
  # points(x = 30.5, y = 3.65, bg = 'orchid3', col = 'black', cex = 1.5, pch = 21)
  add.color.bar(prompt = FALSE, cols = cm.plots[[i]]$cols, leg = 20, 
                title = paste0(names(cm.plots[i]), ' slope'), 
                subtitle = '',
                lwd = 8, lims = cm.plots[[i]]$lims, 
                outline = FALSE, 
                x = 5, y = 3, fsize = 1)
}
dev.off()


# PCA slopes ----
# Using gm.prcomp
# PCA (on correlation matrix)
pca.ols <- gm.prcomp(dat.slp, phy = tree)

# phyPCA (sensu Revell, 2009)
pca.gls <- gm.prcomp(dat.slp, phy = tree, GLS = TRUE)


# PaCA: Align to phylogenetic signal (in each case)
paca <- gm.prcomp(dat.slp, phy = tree, GLS = TRUE, align.to.phy = TRUE)

# get species habitats
hab.sp <- data.sp$habitat_broad[match(rownames(pca.ols$x), data.sp$species)]
names(hab.sp) <- data.sp$species[match(rownames(pca.ols$x), data.sp$species)]

# PCA slopes plot ----
#?plot.gm.prcomp
plot.options <- list(node.labels = FALSE, anc.states = FALSE, 
                     edge.color = 'gray40', edge.width = 0.5)

png('plots/PCA_shape_slopes.png', w = 1200, h = 400)
par(mfrow = c(1,3))
plot(pca.ols, phylo = TRUE, pch = 16, cex = 2, phylo.par = plot.options, 
     col = hab.colors[hab.sp], include.axes = FALSE)
points(x = c(-1.6, -1.6, -1.6), y = c(-1.4, -1.5, -1.6), 
       bg = hab.colors, col = 'transparent', cex = 2, pch = 21)
text(x = c(-1.5, -1.5, -1.5), y = c(-1.4, -1.5, -1.6), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
title(main = 'PCA.OLS')

plot(pca.gls, phylo = TRUE, pch = 16, cex = 2, phylo.par = plot.options,
     col = hab.colors[hab.sp], include.axes = FALSE)
points(x = c(-1.6, -1.6, -1.6), y = c(-1.4, -1.5, -1.6), 
       bg = hab.colors, col = 'transparent', cex = 2, pch = 21)
text(x = c(-1.5, -1.5, -1.5), y = c(-1.4, -1.5, -1.6), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
title(main = 'PCA.GLS')

plot(paca, phylo = TRUE, pch = 16, cex = 2, phylo.par = plot.options,
     col = hab.colors[hab.sp], include.axes = FALSE)
points(x = c(1, 1, 1), y = c(-1.4, -1.5, -1.6), 
       bg = hab.colors, col = 'transparent', cex = 2, pch = 21)
text(x = c(1.1, 1.1, 1.1), y = c(-1.4, -1.5, -1.6), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
title(main = 'PaCA')
dev.off()


# size model (species data) ----
# shape~size phylogenetic regression
shape.sp <- as.matrix(log(dat.morpho.num[, 7:ncol(dat.morpho.num)]))
svl.sp
identical(names(svl.sp), rownames(shape.sp))

# Prepare data
rdf.size <- rrpp.data.frame(svl = svl.sp, shape = shape.sp)

# Multivariate linear model
fit.svl <- lm.rrpp(shape~svl, data = rdf.size, Cov = vcv(tree))
anova(fit.svl)

shape.res <- residuals(fit.svl)

# PCA shape residuals ----
# ordinary PCA
pca.shape.ols <- gm.prcomp(shape.res, phy = tree)

# phylo PCA
pca.shape.gls <- gm.prcomp(shape.res, phy = tree, GLS = TRUE)

# PaCA
paca.shape <- gm.prcomp(shape.res, phy = tree, GLS = TRUE, align.to.phy = TRUE)

# plot PHYLOMORPHOSPACES ----
#hab.sp # habitat states for colors
#hab.colors

pdf('plots/pca_shape.pdf', w = 15, h = 5)
par(mfrow = c(1,3))

# OLS
plot(pca.shape.ols, phylo = TRUE, pch = 16, cex = 2, phylo.par = plot.options, 
     col = hab.colors[hab.sp], include.axes = FALSE)
points(x = c(-0.3, -0.3, -0.3), y = c(-0.18, -0.2, -0.22), 
       bg = hab.colors, col = 'transparent', cex = 2, pch = 21)
text(x = c(-0.28, -0.28, -0.28), y = c(-0.18, -0.2, -0.22), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
title(main = 'pca.shape.OLS')


# GLS
plot(pca.shape.gls, phylo = TRUE, pch = 16, cex = 2, phylo.par = plot.options,
     col = hab.colors[hab.sp], include.axes = FALSE)
points(x = c(-0.22, -0.22, -0.22), y = c(-0.24, -0.26, -0.28), 
       bg = hab.colors, col = 'transparent', cex = 2, pch = 21)
text(x = c(-0.2, -0.2, -0.2), y = c(-0.24, -0.26, -0.28), 
     labels = c('ground', 'rock', 'tree'), adj = 0)
title(main = 'pca.shape.GLS')

# PaCA
plot(paca.shape, phylo = TRUE, pch = 16, cex = 2, phylo.par = plot.options,
     col = hab.colors[hab.sp], include.axes = FALSE)
points(x = c(0.8, 0.8, 0.8), y = c(-1.4, -1.6, -1.8), 
       bg = hab.colors, col = 'transparent', cex = 2, pch = 21)
text(x = c(0.9, 0.9, 0.9), y = c(-1.4, -1.6, -1.8), 
     labels = c('ground', 'rock', 'tree'), adj = 0)

title(main = 'PaCA')

dev.off()

# PHYLOMORPHOSPACE ggplot2 ----
# PLOT shape residuals PCA.OLS with ggplot2
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
source('scripts/ggphylomorpho_HTC.R')
source('scripts/ggphylomorpho_size_HTC.R')
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




######### THE END ############################





