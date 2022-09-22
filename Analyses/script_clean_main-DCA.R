# R-script for primary analyses

libs <- c('geomorph', 'RRPP', 'phytools', 'geiger', 'tidyverse')
easypackages::libraries(libs)

# 0: Data Prep
data0 <- read.table('data/morpho/morpho_pristurus.csv', sep = ';', dec = '.', header = TRUE, stringsAsFactors = TRUE)
  sp.to.keep <- names(which(table(data0$species) >= 5))
data <- data0[data0$species %in% sp.to.keep, ]
  data$species <- droplevels(data$species)
  data$SVL <- log(data$SVL)
shape <- as.matrix(log(data[, 8:ncol(data)]))
rdf <- rrpp.data.frame(svl = data$SVL, shape = shape, habitat = data$habitat_broad, species = data$species)
tree0 <- read.nexus('data/phylogeny/pristurus_tree_final.nex')
LS.mns <- pairwise(lm.rrpp(shape~species, data = rdf, iter=0), groups = rdf$species)$LS.means[[1]]
sz.mn <- tapply(rdf$svl,rdf$species,mean)
hab.mn <- as.factor(by(rdf$habitat,rdf$species,unique))
levels(hab.mn) <- levels(rdf$habitat)
tree <- treedata(phy = tree0, data = LS.mns)$phy
C <- vcv.phylo(tree)

# 1: Evolutionary Allometry
allom.sp <- lm.rrpp(LS.mns~sz.mn, Cov = C)
allom.ind <- lm.rrpp(shape~svl, data = rdf)
anova(allom.sp)
anova(allom.ind)

M <-rbind(coef.sp <- allom.sp$LM$gls.coefficients[2,],
        coef.ind <- allom.ind$LM$coefficients[2,])

acos(RRPP:::vec.cor.matrix(M))*180/pi  #virtually parallel (angle of 1.49 degrees)

# 2: Comparison of multivariate allometry among habitat types
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
  anova(fit.hab)
pw.hab <- pairwise(fit.hab, groups = rdf$habitat, covariate = rdf$svl)
  summary(pw.hab, type = 'VC', stat.table = FALSE)

#Slopes by habitat
fit.coef <- fit.hab$LM$coefficients
ind.coef <- rbind(fit.coef[2,], fit.coef[2,]+fit.coef[5,], fit.coef[2,]+fit.coef[6,])
rownames(ind.coef) <- c("Ground","Rock", "Tree")
ind.coef

# 3: Map allometry slopes on phylogeny
head.scores <- two.b.pls(shape[, c(2:4)], rdf$svl)$XScores[, 1]
limb.scores <- two.b.pls(shape[, 5:8], rdf$svl)$XScores[, 1]

coef.head <- lm.rrpp(head.scores ~ rdf$svl*rdf$species)$LM$coefficients
coef.limb <- lm.rrpp(limb.scores ~ rdf$svl*rdf$species)$LM$coefficients

head.slp <- coef.head[grep('svl', rownames(coef.head)), ]
  head.slp[-1] <- head.slp[-1] + head.slp[1]
limb.slp <- coef.limb[grep('svl', rownames(coef.limb)), ]
  limb.slp[-1] <- limb.slp[-1] + limb.slp[1]
names(limb.slp) <- names(head.slp) <- levels(rdf$species)


contMap(tree = tree, x = head.slp, outline = FALSE)
cm.limb <- contMap(tree = tree, x = limb.slp, outline = FALSE)

# 4: phylomorphospace of size-standardized data (residuals)
shape.res <- residuals(allom.sp)
pca.w.phylo <- gm.prcomp(shape.res, phy = tree)
plot(pca.w.phylo, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))

