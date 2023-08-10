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

acos(RRPP:::vec.cor.matrix(M))*180/pi  #nearly parallel (angle of 5.6 degrees)

# 2: MANCOVA & comparison of allometry among habitats ----
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf2)  #CHANGED 3/22/2023
  anova(fit.hab)

# 2A: Compare habitat vectors versus isometry and to each other
  #H_0: isometry as common slope model
mn.sz <- tapply(rdf2$svl,rdf2$habitat,mean)
mn.shape <- rowsum(rdf2$shape, rdf2$habitat)/as.vector(table(rdf2$habitat)) #CHANGED 3/22/2023
coef.iso <- c(1,1,1,1,1,1,1,1)
intercepts <- mn.shape - t(tcrossprod(coef.iso,mn.sz))
X <- model.matrix(~rdf2$svl+rdf2$habitat)
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
slp.P.ev  #rock different from evol. allometry

res <- cbind(slp.hab.ev.obs[-1,1],slp.Z.ev[-1,1],slp.P.ev[-1,1])
colnames(res) <- c("Angle","Effect Size", "P-value")
rownames(res) <- c("Ev vs. Ground", "Ev vs. Rock", "Ev vs. Tree")
res

# 3: Map allometry slopes on phylogeny ----
 # use multivariate regression
head.scores.old <- two.b.pls(shape[, c(2:4)], rdf$svl)$XScores[, 1]
limb.scores.old <- two.b.pls(shape[, 5:8], rdf$svl)$XScores[, 1]

#CHANGED 3/31/2023 DCA 
#Instead of the PLS, get the regression scores as suggested by reviewer
head.scores <- plot(lm.rrpp(shape[, c(2:4)]~ rdf$svl), 
                 type = "regression", predictor = rdf$svl, reg.type = "RegScore")$RegScore

limb.scores <- plot(lm.rrpp(shape[, c(5:8)]~ rdf$svl), 
                    type = "regression", predictor = rdf$svl, reg.type = "RegScore")$RegScore

#Reg Scores Identical to our PLS score approach
plot(head.scores.old, head.scores)  
plot(limb.scores.old, limb.scores)

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

# ALTERNATIVE TO CONTMAPS: SCATTERPLOT WITH PHYLOGENY ----
slope.dat <- list()
slope.dat$x <- cbind(head.slp,limb.slp)
slope.dat$alignment = "principal"
slope.dat$transform = FALSE
slope.dat$GLS = FALSE
slope.dat$phy <- tree
class(slope.dat) <- "ordinate"

par(mar = c(5,5,2,2))
pdf('../Figs/scatterplot_slopes.pdf')
P2 <- plot(slope.dat, phylo = TRUE, pch = 21, 
           phylo.par = list(node.labels = FALSE), 
           bg = hab.colors[hab.mn[rownames(slope.dat$x)]],
           col = 'black', 
           cex = 1.7, 
           #xlim = c(0, 3), ylim = c(0, 3),
           xlab = 'head slope', ylab = 'limb slope')
add.tree(P2, tree, edge.col = 'black')
text(slope.dat$x[,1], slope.dat$x[,2], labels = rownames(slope.dat$x), 
     #col = hab.colors[hab.mn[rownames(slope.dat$x)]], 
     col = 'black',
     cex = 0.7)
dev.off()


slope.dat$x





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

###### NEW PHYLOMORPHOSPACE: no phylogeny ----
allom2.sp <- lm.rrpp(LS.mns~sz.mn)
shape2.res <- residuals(allom2.sp)
pca.w.phylo2 <- gm.prcomp(shape2.res, phy = tree)
plot(pca.w.phylo2, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))

###### NEW PLOT PHYLOMORPHOSPACE ----
# a) Source the function and the customized theme for the plot
source('function_ggphylomorpho_size.R')
source('function_theme_clean.R')

# b) Prepare data
pca.w.phylo2
x <- rownames_to_column(as.data.frame(pca.w.phylo2$x), 
                        var = 'species')

identical(x$species, names(sz.mn))

x$size <- sz.mn
dat.ggphylo <- cbind(x, hab.mn)

dat.ggphylo <- dat.ggphylo %>% 
  mutate(PC1 = Comp1) %>% 
  mutate(PC2 = Comp2) %>% 
  mutate(taxon = species) %>%
  mutate(group = hab.mn) %>%
  select(taxon, PC1, PC2, group, size)

rownames(dat.ggphylo) <- dat.ggphylo$taxon

# c) Set habitat colors
hab.colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")

# define new colors by adding transparency (to highlight only certain species)
hab.colors.faded <- c(ground = "#F1B670", rock = "#683B5E", 
                      tree = "#E93F7B", 
                      tree_faded = '#E93F7B40', 
                      ground_faded = '#F1B67040', 
                      rock_faded = '#683B5E40')

# d) Set large and small species from rock and ground habitats 
'%nin%' <- Negate('%in%')
large.sp.rock <- c('Pristurus_insignis', 'Pristurus_insignoides')
large.sp.ground <- c('Pristurus_carteri', 'Pristurus_ornithocephalus', 
                     'Pristurus_collaris')
large.sp <- c(large.sp.rock, large.sp.ground)

small.sp.rock <- c('Pristurus_sp5', 'Pristurus_sp1', 'Pristurus_rupestris', 
                   'Pristurus_sp3', 'Pristurus_sp2')
small.sp.ground <- c('Pristurus_masirahensis', 'Pristurus_minimus')
small.sp <- c(small.sp.rock, small.sp.ground)

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


# e) Plot
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
                     breaks = c('ground', 'rock', 'tree_faded'), 
                     labels = c('ground', 'rock', 'tree')) +
  guides(size = 'none') +
  coord_fixed() +
  labs(title = 'Phylomorphospace', 
       subtitle = '(largest and smallest species highlighted)') +
  theme.clean() + 
  theme(legend.title = element_blank())

ggsave('../Figs/new_phylomorphospace.pdf', 
       plot = large_small.phylomorpho.plot)
ggsave('../Figs/new_phylomorphospace.png', 
       plot = large_small.phylomorpho.plot)



##################################### FOR SI ----
# 4b: phylomorphospace of linear measures ----
pca.w.phylo3 <- gm.prcomp(LS.mns, phy = tree)
plot(pca.w.phylo3, phylo = TRUE, pch = 21, bg = 'black', phylo.par = list(node.labels = FALSE))

#PC1 is size*-1
plot(sz.mn,pca.w.phylo3$x[,1])
cor(sz.mn,pca.w.phylo3$x[,1]) #-0.987
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


# Phylogenetic signal of habitat ----
# Borges et al. 2019 Bioinformatics
# https://github.com/mrborges23/delta_statistic
#remotes::install_github("stoufferlab/phyloint")
library(ape)
source('delta_statistic_code.R')
tree # phylogenetic tree
hab.mn # trait vector
identical(names(hab.mn), tree$tip.label) # same order

# estimate delta parameter 
deltaA <- delta(trait = hab.mn, tree = tree, 
                lambda0 = 0.1, se = 0.0589, 
                sim = 10000, thin = 10, burn = 100)

# estimate p-value
nsim <- 1000
random_delta <- rep(NA,nsim)
for (i in 1:nsim){
  print(paste0('simulation ', i))
  rtrait <- sample(hab.mn) # random order of trait values
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")

# histogram p-value
hist(random_delta, breaks = 25, col = 'gray80', 
     main = '')
arrows(x0 = deltaA, x1 = deltaA, y0 = 20, y1 = 0, col = 'red', 
       lwd = 4, length = 0.1)

# No phylogenetic signal of habitat




