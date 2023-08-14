# R-script for primary analyses

libs <- c('geomorph', 'RRPP', 'phytools', 'geiger', 'tidyverse')
easypackages::libraries(libs)

# 0: Data Prep ----
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

### NEW for Revision
SVL.resid <- resid(lm(data$SVL~data$species))
shape.resid <- resid(lm(shape~data$species))

rdf2 <- rrpp.data.frame(svl = SVL.resid, shape = shape.resid, 
                        habitat = data$habitat_broad, 
                        species = data$species)


# 1: Evolutionary Allometry ----
allom.sp <- lm.rrpp(LS.mns~sz.mn, Cov = C)
allom.ind <- lm.rrpp(shape~svl, data = rdf2)  #CHANGED 3/22/2023
anova(allom.sp)
anova(allom.ind)

M <-rbind(coef.evol <- allom.sp$LM$gls.coefficients[2,],
        coef.ind <- allom.ind$LM$coefficients[2,])

acos(RRPP:::vec.cor.matrix(M))*180/pi  #virtually parallel (angle of 1.49 degrees)
  #ANGLE 5.6 degrees

# 2: MANCOVA & comparison of allometry among habitats ----
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf2)  #CHANGED 3/22/2023
  anova(fit.hab)

# 2A: Compare habitat vectors versus isometry and to each other
  #H_0: isometry as common slope model
mn.sz <- tapply(rdf$svl,rdf$habitat,mean)
mn.shape <- rowsum(rdf2$shape, rdf$habitat)/as.vector(table(rdf$habitat)) #CHANGED 3/22/2023
coef.iso <- c(1,1,1,1,1,1,1,1)
intercepts <- mn.shape - t(tcrossprod(coef.iso,mn.sz))
X <- model.matrix(~rdf$svl+rdf$habitat)
b <- rbind(intercepts[1,],coef.iso,intercepts[2,]-intercepts[1,],
           intercepts[3,]-intercepts[1,])
preds <- X%*%b
E.iso <- rdf2$shape - preds 
perms <- RRPP:::perm.index(n = fit.hab$LM$n, iter = 999)
slopes <- list()
for(j in 1:1000){ #CHANGED 3/22/2023
  slopes[[j]] <- pairwise(lm.rrpp((preds+E.iso[perms[[j]],]) ~ rdf2$svl*rdf2$habitat, iter=0), 
           groups = rdf2$habitat,covariate = rdf2$svl)$slopes[[1]]
  
}  
slp.ang <- lapply(1:1000, function(j) acos(RRPP:::vec.cor.matrix(rbind(slopes[[j]],coef.iso)))*180/pi)

slp.hab.obs <- slp.ang[[1]]
slp.Z <- RRPP:::effect.list(slp.ang)
slp.P <- RRPP:::Pval.list(slp.ang)

slp.hab.obs
slp.Z
slp.P #all different from isometry, and ground different from rock and tree

# 2B: Compare evolutionary and static (habitat) allometry 
  #H_0: common slope isometry
slp.ang.ev <- lapply(1:1000, function(j) acos(RRPP:::vec.cor.matrix(rbind(coef.evol,slopes[[j]])))*180/pi)

slp.hab.ev.obs <- slp.ang.ev[[1]]
slp.Z.ev <- RRPP:::effect.list(slp.ang.ev)
slp.P.ev <- RRPP:::Pval.list(slp.ang.ev)

slp.hab.ev.obs
slp.Z.ev
slp.P.ev  #not different from evol. allometry

res <- cbind(slp.hab.ev.obs[-1,1],slp.Z.ev[-1,1],slp.P.ev[-1,1])
colnames(res) <- c("Angle","Effect Size", "P-value")
rownames(res) <- c("Ev vs. Ground", "Ev vs. Rock", "Ev vs. Tree")
res

########################################  DID NOT IMPLEMENT YET

### NOTE: SEE REV 2. May just simplify this one and description

# 3: Map allometry slopes on phylogeny ----
#head.scores <- two.b.pls(shape[, c(2:4)], rdf2$svl)$XScores[, 1]
#limb.scores <- two.b.pls(shape[, 5:8], rdf2$svl)$XScores[, 1]

head.scores <- plot(lm.rrpp(shape[, c(2:4)]~ rdf$svl),
                    type = "regression", predictor = rdf$svl, reg.type = "RegScore")$RegScore
limb.scores <- plot(lm.rrpp(shape[, c(5:8)]~ rdf$svl),
                    type = "regression", predictor = rdf$svl, reg.type = "RegScore")$RegScore

coef.head <- lm.rrpp(head.scores ~ rdf$svl*rdf$species)$LM$coefficients
coef.limb <- lm.rrpp(limb.scores ~ rdf$svl*rdf$species)$LM$coefficients

head.slp <- coef.head[grep('svl', rownames(coef.head)), ]
  head.slp[-1] <- head.slp[-1] + head.slp[1]
limb.slp <- coef.limb[grep('svl', rownames(coef.limb)), ]
  limb.slp[-1] <- limb.slp[-1] + limb.slp[1]
names(limb.slp) <- names(head.slp) <- levels(rdf$species)

cor(head.slp,limb.slp)
plot(head.slp,limb.slp)

cm.head <- contMap(tree = tree, x = head.slp, outline = FALSE)
cm.limb <- contMap(tree = tree, x = limb.slp, outline = FALSE)

slope.dat$x <- cbind(head.slp,limb.slp)    
  slope.dat$alignment = "principal"
  slope.dat$transform = "FALSE"
  slope.dat$GLS = "FALSE"
  slope.dat$phy <- tree
class(slope.dat) <- "ordinate"

P2 <- plot(slope.dat, phylo = TRUE, pch = 21, bg = 'black', 
           phylo.par = list(node.labels = FALSE))
add.tree(P2, tree, edge.col = 4)


################################


# 4: Compare Integration ----  #CHANGED 3/22/2023
lindims.gp <- lapply( split( shape.resid[,1:ncol(shape.resid)], rdf$habitat), matrix, ncol=ncol(shape.resid))
Vrel.gp <- Map(function(x) integration.Vrel(x), lindims.gp) 
c(Vrel.gp$ground$ZR,Vrel.gp$rock$ZR,Vrel.gp$tree$ZR)
out <- compare.ZVrel(Vrel.gp$ground, Vrel.gp$rock, Vrel.gp$tree)
summary(out)

########################################### DID NOT IMPLEMENT YET (IS IT NEEDED?)

# 5: phylomorphospace of size-standardized data (residuals) ----
shape.res <- residuals(allom.sp)
pca.w.phylo <- gm.prcomp(shape.res, phy = tree)
plot(pca.w.phylo, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))

###### NEW PHYLOMORPHOSPACE: no phylogeny
allom2.sp <- lm.rrpp(LS.mns~sz.mn)
shape2.res <- residuals(allom2.sp)
pca.w.phylo2 <- gm.prcomp(shape2.res, phy = tree)
plot(pca.w.phylo2, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))



##################################### FOR SI ----
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
