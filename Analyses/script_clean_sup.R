
# Supplementary analyses and figures

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
# Set habitat colors
hab.colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")

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






