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

M <-rbind(coef.evol <- allom.sp$LM$gls.coefficients[2,],
        coef.ind <- allom.ind$LM$coefficients[2,])

acos(RRPP:::vec.cor.matrix(M))*180/pi  #virtually parallel (angle of 1.49 degrees)

# 2: MANCOVA & comparison of allometry among habitats 
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
  anova(fit.hab)
pw.hab1 <- pairwise(fit.hab, groups = rdf$habitat, covariate = rdf$svl)
  summary(pw.hab1, type = 'VC', stat.table = FALSE)
  
# 2A: Compare evolutionary and static (habitat) allometry 
fit.0 <- lm.rrpp(shape~habitat, data = rdf) #H_0: group differences only (no allometry) 
pw.hab <- pairwise(fit.hab, fit.null = fit.0, groups = rdf$habitat, covariate = rdf$svl)
summary(pw.hab, type = 'VC', stat.table = FALSE)
  slp.hab <- pw.hab$slopes[[1]]

#Slopes by habitat
slopes <- lapply(1:1000, function(j) rbind(coef.evol,pw.hab$slopes[[j]]))
slp.ang <- lapply(1:1000, function(j) acos(RRPP:::vec.cor.matrix(slopes[[j]]))*180/pi)

slp.hab.obs <- slp.ang[[1]]
slp.Z <- RRPP:::effect.list(slp.ang)
slp.P <- RRPP:::Pval.list(slp.ang)

slp.hab.obs
slp.Z
slp.P  #angles significantly smaller than expected

res <- cbind(slp.hab.obs[-1,1],slp.Z[-1,1],slp.P[-1,1])
colnames(res) <- c("Angle","Effect Size", "P-value")
rownames(res) <- c("Ground", "Rock", "Tree")
res

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

cor(head.slp,limb.slp)
plot(head.slp,limb.slp)

contMap(tree = tree, x = head.slp, outline = FALSE)
cm.limb <- contMap(tree = tree, x = limb.slp, outline = FALSE)

# 4: phylomorphospace of size-standardized data (residuals)
shape.res <- residuals(allom.sp)
pca.w.phylo <- gm.prcomp(shape.res, phy = tree)
plot(pca.w.phylo, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))

################################ NEW HERE

# 5: Compare Integration
lindims.gp <- lapply( split( shape[,1:ncol(shape)], rdf$habitat), matrix, ncol=ncol(shape))
Vrel.gp <- Map(function(x) integration.Vrel(x), lindims.gp) 
c(Vrel.gp$ground$ZR,Vrel.gp$rock$ZR,Vrel.gp$tree$ZR)
out <- compare.ZVrel(Vrel.gp$ground, Vrel.gp$rock, Vrel.gp$tree)
summary(out)


## plot
library(gplots)
Z.gp <- c(Vrel.gp$ground$ZR,Vrel.gp$rock$ZR,Vrel.gp$tree$ZR)
Z.var <- c(Vrel.gp$ground$ZR.var,Vrel.gp$rock$ZR.var,Vrel.gp$tree$ZR.var)
CI.gp<-qnorm(1- 0.05/2) * sqrt(Z.var)

plotCI(Z.gp,ui =(Z.gp+CI.gp), li = (Z.gp-CI.gp), 
       ylab="Integration Level (Z-Vrel)", xlab="",pch=21,cex=2, 
       xaxt='n',pt.bg="black",lwd=3,lty=1)
axis(1, at = 1:3,
     labels = c("Ground",
                "Rock",
                "Tree"))

##################################### FOR SI
# 4b: phylomorphospace of linear measures
pca.w.phylo2 <- gm.prcomp(LS.mns, phy = tree)
plot(pca.w.phylo2, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))

#PC1 is size*-1
plot(sz.mn,pca.w.phylo2$x[,1])
cor(sz.mn,pca.w.phylo2$x[,1]) #-0.987
#size drives the show

# Disparity among habitat groups
fit.lm <- lm.rrpp(LS.mns~hab.mn, Cov =C)
PW.lm <- pairwise(fit.lm,groups = hab.mn)
summary(PW.lm, test.type = "var")

# Integration of size-standardized variables (similar patterns; lower integration)
shape.2 <- shape - rdf$svl
shape.gp <- lapply( split( shape.2[,1:ncol(shape.2)], rdf$habitat), matrix, ncol=ncol(shape.2))
Vrel.gp.shp <- Map(function(x) integration.Vrel(x), shape.gp) 
out.shp <- compare.ZVrel(Vrel.gp.shp$ground, Vrel.gp.shp$rock, Vrel.gp.shp$tree)
summary(out.shp)

PCA.gp <- Map(function(x) prcomp(x), lindims.gp)
PCA.gp.shp <- Map(function(x) prcomp(x), shape.gp)
PCA.gp$ground$sdev
PCA.gp.shp$ground$sdev

##### Histograms of size: largest size range is in ground dwelling 
svl.gp <- split(rdf$svl, rdf$habitat)
range(svl.gp$ground)

hist(svl.gp$ground,xlim = c(2.5,4.5), ylim = c(0,150),xlab = "log(SVL)", col = "black")
hist(svl.gp$rock,add = TRUE, col="blue")
hist(svl.gp$tree,add = TRUE, col="red")
