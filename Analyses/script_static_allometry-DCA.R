
libs <- c('geomorph', 'RRPP', 'phytools', 'geiger', 'tidyverse')
easypackages::libraries(libs)

# Evolutionary allometry ---- 
# (from individuals dataset)
data0 <- read.table('data/morpho/morpho_pristurus.csv', sep = ';', dec = '.', header = TRUE, stringsAsFactors = TRUE)
sp.to.keep <- names(which(table(data0$species) >= 5))
data <- data0[data0$species %in% sp.to.keep, ]
data$species <- droplevels(data$species)
data$SVL <- log(data$SVL)
species.fctr <- as.factor(data$species)
habitat.fctr <- as.factor(data$habitat_broad)
shape <- as.matrix(log(data[, 8:ncol(data)]))

rdf <- rrpp.data.frame(svl = data$SVL, shape = shape, 
                       habitat = data$habitat_broad, 
                       species = data$species)

####################### DCA added from other script  #####################
tree0 <- read.nexus('data/phylogeny/pristurus_tree_final.nex')
LS.mns <- pairwise(lm.rrpp(shape~species, data = rdf, iter=0), groups = rdf$species)$LS.means[[1]]
sz.mn <- tapply(rdf$svl,rdf$species,mean)
hab.mn <- as.factor(by(rdf$habitat,rdf$species,unique))
levels(hab.mn) <- levels(rdf$habitat)
tree <- treedata(phy = tree0, data = LS.mns)$phy
C <- vcv.phylo(tree)

######   EVOLUTIONARY ALLOMETRY (is based on species means)
allom.evol <- lm.rrpp(LS.mns~sz.mn, Cov = C)
allom.ind <- lm.rrpp(shape~svl, data = rdf)  #overall individual allometry

#comparison of angles
M <-rbind(coef.evol <- allom.evol$LM$gls.coefficients[2,],
          coef.ind <- allom.ind$LM$coefficients[2,])
acos(RRPP:::vec.cor.matrix(M))*180/pi  #virtually parallel (angle of 1.49 degrees)

#Static allometry in each species 
fit.sp <- lm.rrpp(shape~svl*species, data = rdf)
fit.coef.sp <- fit.sp$LM$coefficients
sp.slp <- fit.coef.sp[grep('svl', rownames(fit.coef.sp)), ]
sp.slp[-1,] <- sp.slp[-1,] + sp.slp[1,]
rownames(sp.slp) <- gsub('svl:species', '', rownames(sp.slp))
rownames(sp.slp) <- gsub('svl', 'Pristurus_abdelkuri', rownames(sp.slp))

# 1: Compare evol and static allometry ----
slp.comp <-rbind(coef.evol, sp.slp)
slp.angles <- acos(RRPP:::vec.cor.matrix(slp.comp)[,1])*180/pi
mean(slp.angles)
hist(slp.angles)

sort(slp.angles, decreasing = FALSE)

# 2: Compare evol and static (habitat) allometry ----
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
coef.hab <- fit.hab$LM$coefficients
slp.hab <- rbind(coef.hab[2,], coef.hab[2,]+coef.hab[5,], 
                 coef.hab[2,]+coef.hab[6,])
rownames(slp.hab) <- c("ground","rock", "tree")
slp.hab <-rbind(coef.evol, slp.hab)
slp.hab.angles <- acos(RRPP:::vec.cor.matrix(slp.hab)[,1])*180/pi
mean(slp.hab.angles)
hist(slp.hab.angles)
slp.hab.angles

#####################################################################
#################################   From the script you sent #########
allom.ind <- lm.rrpp(shape~svl, data = rdf)
coef.ind <- allom.ind$LM$coefficients[2,]
coef.ind

# Static allometry ----
fit.sp <- lm.rrpp(shape~svl*species, data = rdf)
fit.coef.sp <- fit.sp$LM$coefficients
sp.slp <- fit.coef.sp[grep('svl', rownames(fit.coef.sp)), ]
sp.slp[-1,] <- sp.slp[-1,] + sp.slp[1,]
rownames(sp.slp) <- gsub('svl:species', '', rownames(sp.slp))
rownames(sp.slp) <- gsub('svl', 'Pristurus_abdelkuri', rownames(sp.slp))

# Compare evol and static allometry ----
slp.comp <-rbind(evol_alom = coef.ind, sp.slp)
slp.correlations <- vec.cor.matrix(slp.comp)[1,]
slp.angles <- (acos(RRPP:::vec.cor.matrix(slp.comp))*180/pi)[1,]

sort(slp.correlations, decreasing = TRUE)





# Evolutionary allometry PER HABITAT ---- 
# (from individuals dataset)
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
coef.hab <- fit.hab$LM$coefficients
slp.hab <- rbind(coef.hab[2,], coef.hab[2,]+coef.hab[5,], 
                 coef.hab[2,]+coef.hab[6,])
rownames(slp.hab) <- c("ground","rock", "tree")
slp.hab

# Static allometry PER HABITAT ----
sp.slp

hab.mn <- as.factor(by(rdf$habitat,rdf$species,unique))
levels(hab.mn) <- levels(rdf$habitat)

identical(rownames(sp.slp), names(hab.mn))

sp.slp[hab.mn == 'ground', ]
sp.slp[hab.mn == 'rock', ]
sp.slp[hab.mn == 'tree', ]

# Compare evol and static allometry PER HABITAT ----
slp.comp.ground <-rbind(evol_alom = slp.hab['ground',], sp.slp[hab.mn == 'ground', ])
slp.correlations.ground <- vec.cor.matrix(slp.comp.ground)[1,]
slp.angles.ground <- (acos(RRPP:::vec.cor.matrix(slp.comp.ground))*180/pi)[1,]

slp.comp.rock <-rbind(evol_alom = slp.hab['rock',], sp.slp[hab.mn == 'rock', ])
slp.correlations.rock <- vec.cor.matrix(slp.comp.rock)[1,]
slp.angles.rock <- (acos(RRPP:::vec.cor.matrix(slp.comp.rock))*180/pi)[1,]

slp.comp.tree <-rbind(evol_alom = slp.hab['tree',], sp.slp[hab.mn == 'tree', ])
slp.correlations.tree <- vec.cor.matrix(slp.comp.tree)[1,]
slp.angles.tree <- (acos(RRPP:::vec.cor.matrix(slp.comp.tree))*180/pi)[1,]


# SEPARATE HEAD AND LIMBS ----
svl <- log(data$SVL)
shape.head <- shape[, c(2:4)]
shape.limb <- shape[, 5:8]

# EVOLUTIONARY ALLOMETRY head and limbs ----
# Linear Models and slopes ----
# lm(PLS.scores ~ svl*habitat)
fit.headshape.hab <- lm.rrpp(shape.head ~ svl*habitat.fctr)
coef.headshape.hab <- fit.headshape.hab$LM$coefficients
slp.headshape.hab <- rbind(coef.headshape.hab[2,], 
                           coef.headshape.hab[2,]+coef.headshape.hab[5,], 
                           coef.headshape.hab[2,]+coef.headshape.hab[6,])
rownames(slp.headshape.hab) <- c("ground","rock", "tree")
slp.headshape.hab

fit.limbshape.hab <- lm.rrpp(shape.limb ~ svl*habitat.fctr)
coef.limbshape.hab <- fit.limbshape.hab$LM$coefficients
slp.limbshape.hab <- rbind(coef.limbshape.hab[2,], 
                           coef.limbshape.hab[2,]+coef.limbshape.hab[5,], 
                           coef.limbshape.hab[2,]+coef.limbshape.hab[6,])
rownames(slp.limbshape.hab) <- c("ground","rock", "tree")
slp.limbshape.hab



# STATIC ALLOMETRY head and limbs ----
# Linear Models ----
# lm(shape.vars ~ svl*sp)
fit.headshape.sp <- lm.rrpp(shape.head ~ svl*species.fctr)
plot.headshape.fit <- plot(fit.headshape.sp, predictor = svl, type = 'regression', pch = 16, 
                           col = species.fctr)
plot.headshape.fit$PredLine
plot.headshape.fit$RegScore
plot(fit.headshape.sp, predictor = svl, type = 'regression', pch = 16, 
     reg.type = 'RegScore', 
     col = species.fctr)

fit.limbshape.sp <- lm.rrpp(shape.limb ~ svl*species.fctr)
plot.limbshape.fit <- plot(fit.limbshape.sp, predictor = svl, type = 'regression', pch = 16, 
                           col = species.fctr)
plot.limbshape.fit$PredLine
plot.limbshape.fit$RegScore
plot(fit.limbshape.sp, predictor = svl, type = 'regression', pch = 16, 
     reg.type = 'RegScore', 
     col = species.fctr)

# Regression coefficients ----
coef.headshape <- fit.headshape.sp$LM$coefficients
coef.limbshape <- fit.limbshape.sp$LM$coefficients

# Per-species slopes ----
headshape.slp <- coef.headshape[grep('svl', rownames(coef.headshape)), ]
headshape.slp[-1,] <- headshape.slp[-1,] + headshape.slp[1,]
rownames(headshape.slp) <- gsub('svl:species.fctr', '', rownames(headshape.slp))
rownames(headshape.slp) <- gsub('svl', 'Pristurus_abdelkuri', rownames(headshape.slp))

limbshape.slp <- coef.limbshape[grep('svl', rownames(coef.limbshape)), ]
limbshape.slp[-1,] <- limbshape.slp[-1,] + limbshape.slp[1,]
rownames(limbshape.slp) <- gsub('svl:species.fctr', '', rownames(limbshape.slp))
rownames(limbshape.slp) <- gsub('svl', 'Pristurus_abdelkuri', rownames(limbshape.slp))

# Compare evol and static allometry head and limbs per habitat ----
headshape.slp
limbshape.slp
slp.headshape.hab
slp.limbshape.hab

# head ground
slp.headshape.ground <-rbind(evol_alom = slp.headshape.hab['ground',], 
                             headshape.slp[hab.mn == 'ground', ])
slp.headshape.correlations.ground <- vec.cor.matrix(slp.headshape.ground)[1,]
slp.headshape.angles.ground <- (acos(RRPP:::vec.cor.matrix(slp.headshape.ground))*180/pi)[1,]

# head rock
slp.headshape.rock <-rbind(evol_alom = slp.headshape.hab['rock',], 
                           headshape.slp[hab.mn == 'rock', ])
slp.headshape.correlations.rock <- vec.cor.matrix(slp.headshape.rock)[1,]
slp.headshape.angles.rock <- (acos(RRPP:::vec.cor.matrix(slp.headshape.rock))*180/pi)[1,]

# head tree
slp.headshape.tree <-rbind(evol_alom = slp.headshape.hab['tree',], 
                           headshape.slp[hab.mn == 'tree', ])
slp.headshape.correlations.tree <- vec.cor.matrix(slp.headshape.tree)[1,]
slp.headshape.angles.tree <- (acos(RRPP:::vec.cor.matrix(slp.headshape.tree))*180/pi)[1,]


# limb ground
slp.limbshape.ground <-rbind(evol_alom = slp.limbshape.hab['ground',], 
                             limbshape.slp[hab.mn == 'ground', ])
slp.limbshape.correlations.ground <- vec.cor.matrix(slp.limbshape.ground)[1,]
slp.limbshape.angles.ground <- (acos(RRPP:::vec.cor.matrix(slp.limbshape.ground))*180/pi)[1,]

# limb rock
slp.limbshape.rock <-rbind(evol_alom = slp.limbshape.hab['rock',], 
                           limbshape.slp[hab.mn == 'rock', ])
slp.limbshape.correlations.rock <- vec.cor.matrix(slp.limbshape.rock)[1,]
slp.limbshape.angles.rock <- (acos(RRPP:::vec.cor.matrix(slp.limbshape.rock))*180/pi)[1,]

# limb tree
slp.limbshape.tree <-rbind(evol_alom = slp.limbshape.hab['tree',], 
                           limbshape.slp[hab.mn == 'tree', ])
slp.limbshape.correlations.tree <- vec.cor.matrix(slp.limbshape.tree)[1,]
slp.limbshape.angles.tree <- (acos(RRPP:::vec.cor.matrix(slp.limbshape.tree))*180/pi)[1,]

# TABLE CORRELATION AND ANGLES ----
# Get a table with the correlation and angle of each species to the 
# evolutionary allometry of their respective habitat
cor.angle.table <- data.frame(species = c(names(slp.headshape.angles.ground), 
                                          names(slp.headshape.angles.rock), 
                                          names(slp.headshape.angles.tree)), 
                              hab = c(rep('ground', length(slp.headshape.angles.ground)), 
                                      rep('rock', length(slp.headshape.angles.rock)), 
                                      rep('tree', length(slp.headshape.angles.tree))),
                              cor_head = c(slp.headshape.correlations.ground, 
                                           slp.headshape.correlations.rock,
                                           slp.headshape.correlations.tree),
                              angle_head = c(slp.headshape.angles.ground, 
                                             slp.headshape.angles.rock, 
                                             slp.headshape.angles.tree),
                              cor_limb = c(slp.limbshape.correlations.ground, 
                                           slp.limbshape.correlations.rock,
                                           slp.limbshape.correlations.tree),
                              angle_limb = c(slp.limbshape.angles.ground, 
                                             slp.limbshape.angles.rock, 
                                             slp.limbshape.angles.tree))
cor.angle.table <- cor.angle.table[cor.angle.table$species != 'evol_alom', ]

write.table(cor.angle.table, '../Tables/cor_angles.csv', sep = ';', 
            dec = '.', quote = FALSE, row.names = FALSE)

# COMPARE BETWEEN SPECIES ----
# Sort by habitat
slp.sort <- sp.slp[as.character(cor.angle.table$species),]

# Get correlation matrix between species 
slp.sp.cor <- round(vec.cor.matrix(slp.sort), 2)
slp.sp.angle <- round(acos(RRPP:::vec.cor.matrix(slp.sort))*180/pi, 2)

write.table(slp.sp.cor, '../Tables/cor_between_species_slopes.csv', sep = ';', 
            dec = '.', quote = FALSE, row.names = TRUE, col.names = NA)

write.table(slp.sp.angle, '../Tables/angles_between_species_slopes.csv', sep = ';', 
            dec = '.', quote = FALSE, row.names = TRUE, col.names = NA)




# PLS ----
# ... PLS head vs size ----
svl <- log(data$SVL)
shape.head <- shape[, c(2:4)]

pls.head <- two.b.pls(shape.head, svl)

# ... PLS limbs vs size ----
shape.limb <- shape[, 5:8]
pls.limb <- two.b.pls(shape.limb, svl)

# ... PLS whole shape vs size ----
pls.shape <- two.b.pls(shape, svl)

# Get PLS scores ----
limb.scores <- pls.limb$XScores[, 1]
head.scores <- pls.head$XScores[, 1]
shape.scores <- pls.shape$XScores[, 1]

hab.colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")

par(mfrow=c(1,3))
plot(svl, pls.limb$XScores[, 1], pch = 22, bg = hab.colors[habitat.fctr])
plot(svl, pls.head$XScores[, 1], pch = 25, 
     bg = hab.colors[habitat.fctr])
plot(svl, pls.shape$XScores[, 1], pch = 21, 
     bg = hab.colors[habitat.fctr])

# SLOPES ----
# LM scores ~ svl*sp 
# Get a slope per species for head and limbs 
# (instead of for each variable as we did before)

fit.head.sp <- lm.rrpp(head.scores ~ svl*species.fctr)
fit.limb.sp <- lm.rrpp(limb.scores ~ svl*species.fctr)
fit.shape.sp <- lm.rrpp(shape.scores ~ svl*species.fctr)
anova(fit.head.sp)
anova(fit.limb.sp)
anova(fit.shape.sp)

plot.pls.head <- plot(fit.head.sp, predictor = svl, type = 'regression', pch = 16, 
                      col = species.fctr)
plot.pls.head$PredLine
plot.pls.head$RegScore
plot(fit.head.sp, predictor = svl, type = 'regression', pch = 16, 
     reg.type = 'RegScore', 
     col = species.fctr)

plot.pls.limb <- plot(fit.limb.sp, predictor = svl, type = 'regression', pch = 16, 
                      col = species.fctr)
plot.pls.limb$PredLine
plot.pls.limb$RegScore
plot(fit.limb.sp, predictor = svl, type = 'regression', pch = 16, 
     reg.type = 'RegScore', 
     col = species.fctr)


plot(fit.shape.sp, predictor = svl, type = 'regression', pch = 16, 
     col = species.fctr)

# regression coefficients
coef.head <- fit.head.sp$LM$coefficients
coef.limb <- fit.limb.sp$LM$coefficients
coef.shape <- fit.shape.sp$LM$coefficients

# get per-species slopes
head.slp <- coef.head[grep('svl', rownames(coef.head)), ]
head.slp[-1] <- head.slp[-1] + head.slp[1]
names(head.slp) <- gsub('svl:species.fctr', '', names(head.slp))
names(head.slp) <- gsub('svl', 'Pristurus_abdelkuri', names(head.slp))

limb.slp <- coef.limb[grep('svl', rownames(coef.limb)), ]
limb.slp[-1] <- limb.slp[-1] + limb.slp[1]
names(limb.slp) <- gsub('svl:species.fctr', '', names(limb.slp))
names(limb.slp) <- gsub('svl', 'Pristurus_abdelkuri', names(limb.slp))

shape.slp <- coef.shape[grep('svl', rownames(coef.shape)), ]
shape.slp[-1] <- shape.slp[-1] + shape.slp[1]
names(shape.slp) <- gsub('svl:species.fctr', '', names(shape.slp))
names(shape.slp) <- gsub('svl', 'Pristurus_abdelkuri', names(shape.slp))



# PLOT STATIC ALLOMETRY PER HABITAT ----

fit.sp.ggplot.data <- data.frame(pred_head = plot.pls.head$PredLine, 
                                 regscore_head = plot.pls.head$RegScore, 
                                 slope_head = head.slp[rdf$species],
                                 pred_limb = plot.pls.limb$PredLine, 
                                 regscore_limb = plot.pls.limb$RegScore, 
                                 slope_limb = limb.slp[rdf$species],
                                 svl = rdf$svl,
                                 habitat = rdf$habitat, 
                                 sp = rdf$species)

range(fit.sp.ggplot.data$slope_head)

# Set the color palette
ramp <- colorRampPalette(c("#00929c","gray80",  "#d62e31"))


head_static_plot <- ggplot(data = fit.sp.ggplot.data, aes(x = svl)) +
  facet_wrap(.~habitat, ncol = 3) +
  geom_point(aes(y = regscore_head, color = as.factor(slope_head)), 
             alpha = 0.1, size = 2) +
  geom_line(aes(y = pred_head, color = as.factor(slope_head)), size = 1) +
  #  geom_line(aes(y = iso.line), lty = 'dashed', size = 0.5) +
  #scale_color_manual(values = hab.colors) +
  scale_color_manual(values = ramp(25)) +
  labs(x = 'logSVL', y = 'Regression Scores',
       title = 'HEAD STATIC ALLOMETRY') +
  theme_bw() +
  theme(legend.position = '', 
        plot.title = element_text(size = 15, vjust = 1, hjust = 0.5, 
                                  face = 'bold')
  )

limb_static_plot <- ggplot(data = fit.sp.ggplot.data, aes(x = svl)) +
  facet_wrap(.~habitat, ncol = 3) +
  geom_point(aes(y = regscore_limb, color = as.factor(slope_limb)), 
             alpha = 0.1, size = 2) +
  geom_line(aes(y = pred_limb, color = as.factor(slope_limb)), size = 1) +
  #  geom_line(aes(y = iso.line), lty = 'dashed', size = 0.5) +
  #scale_color_manual(values = hab.colors) +
  scale_color_manual(values = ramp(25)) +
  labs(x = 'logSVL', y = 'Regression Scores',
       title = 'LIMB STATIC ALLOMETRY') +
  theme_bw() +
  theme(legend.position = '', 
        plot.title = element_text(size = 15, vjust = 1, hjust = 0.5, 
                                  face = 'bold')
  )

library(patchwork)

static_allometry <- head_static_plot / limb_static_plot

ggsave('../Figs/static_allometry.png', plot = static_allometry)


