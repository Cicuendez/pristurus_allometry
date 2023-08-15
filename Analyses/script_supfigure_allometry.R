
# purpose ----
# Visualize the different types of allometry used in the paper 
# for explanatory purposes.

# packages ----
libs <- c('geomorph', 'RRPP', 'phytools', 'geiger', 'tidyverse', 
          'ggphylomorpho')
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

# Set title size for plots ----
title_sz <- 12

# Evolutionary allometry ----
allom.sp <- lm.rrpp(LS.mns~sz.mn, Cov = C)
allom.sp.plot <- plot(allom.sp, predictor = as.numeric(sz.mn), 
                      type = 'regression', reg.type = 'RegScore', pch = 16)

gg_df_evol <- data.frame(pred = -allom.sp.plot$PredLine, 
           RegScore = allom.sp.plot$RegScore[,1], 
           SVL = sz.mn)

gg_evol <- ggplot(data = gg_df_evol, aes(x = SVL)) +
  geom_point(aes(y = RegScore), color = 'transparent',
             fill = 'gray60', size = 5, pch = 21, alpha = 0.7) +
  geom_line(aes(y = pred), color = 'gray30', size = 0.8) +
  labs(title = 'Evolutionary allometry', x = 'size', y = 'shape') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = title_sz))

# Individual allometry ----
allom.ind <- lm.rrpp(shape~svl, data = rdf)
allom.ind.plot <- plot(allom.ind, predictor = rdf$svl, 
                       type = 'regression', reg.type = 'RegScore', pch = 16)
gg_df_ind <- data.frame(pred = allom.ind.plot$PredLine, 
                        RegScore = allom.ind.plot$RegScore[,1], 
                        SVL = rdf$svl)

gg_ind <- ggplot(data = gg_df_ind, aes(x = SVL)) +
  geom_point(aes(y = RegScore),
             color = 'gray60', size = 1, pch = 16, alpha = 0.7) +
  geom_line(aes(y = pred), color = 'black', size = 0.8) +
  labs(title = 'Individual allometry', x = 'size', y = 'shape') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = title_sz))


# Evolutionary allometry -- Habitat ----
allom.sp.hab <- lm.rrpp(LS.mns~sz.mn*hab.mn, Cov = C)
allom.sp.hab.plot <- plot(allom.sp.hab, predictor = as.numeric(sz.mn), 
                      type = 'regression', reg.type = 'RegScore', pch = 16)
gg_df_evol_hab <- data.frame(pred = -allom.sp.hab.plot$PredLine, 
                             RegScore = allom.sp.hab.plot$RegScore[,1], 
                             SVL = sz.mn, 
                             habitat = hab.mn)

hab.colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")

gg_evol_hab <- ggplot(data = gg_df_evol_hab, aes(x = SVL)) +
  geom_point(aes(y = RegScore, fill = habitat), color = 'transparent',
             size = 5, pch = 21, alpha = 0.7) +
  geom_line(aes(y = pred, color = habitat), size = 1) +
  scale_fill_manual(values = hab.colors) +
  scale_color_manual(values = hab.colors) +
  labs(title = 'Evolutionary allometry per habitat', x = 'size', y = 'shape') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = title_sz), 
        legend.position = 'none')


# Intraspecific allometry -- Habitat ----
allom.ind.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
allom.ind.hab.plot <- plot(allom.ind.hab, predictor = rdf$svl, 
                           type = 'regression', pch = 16, 
                           col = hab.colors[as.numeric(rdf$habitat)])

gg_df_ind_hab <- data.frame(pred = allom.ind.hab.plot$PredLine, 
                            RegScore = allom.ind.hab.plot$RegScore[,1], 
                            SVL = rdf$svl, 
                            habitat = rdf$habitat)

gg_ind_hab <- ggplot(data = gg_df_ind_hab, aes(x = SVL)) +
  geom_point(aes(y = RegScore, color = habitat), 
             size = 1, pch = 16, alpha = 0.7) +
  geom_line(aes(y = pred, color = habitat), size = 1) +
  scale_color_manual(values = hab.colors) +
  labs(title = 'Habitat-based individual allometry', x = 'size', y = 'shape') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = title_sz), 
        legend.position = 'none')

# Intraspecific allometry -- Species ----
allom.ind.sp <- lm.rrpp(shape~svl*species, data = rdf)
allom.ind.sp.plot <- plot(allom.ind.sp, predictor = rdf$svl, 
                           type = 'regression', pch = 16, 
                           col = hab.colors[as.numeric(rdf$habitat)])
gg_df_ind_sp <- data.frame(pred = allom.ind.sp.plot$PredLine, 
                           RegScore = allom.ind.sp.plot$RegScore[,1], 
                           SVL = rdf$svl, 
                           species = rdf$species,
                           habitat = rdf$habitat)

gg_ind_sp <- ggplot(data = gg_df_ind_sp, aes(x = SVL)) +
  geom_point(aes(y = RegScore, color = habitat), 
             size = 1, pch = 16, alpha = 0.7) +
  geom_line(aes(y = pred, color = habitat, group = species), 
            size = 1) +
  scale_color_manual(values = hab.colors) +
  labs(title = 'Intraspecific allometry', 
       subtitle = '(colored by habitat)', x = 'size', y = 'shape') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = title_sz), 
        plot.subtitle = element_text(hjust = 0.5, face = 'plain'), 
        legend.position = 'none')

# Plot everything
library(patchwork)
layout <- '
AB
CD
#E
'
wrap_plots(A = gg_evol, B = gg_ind, 
           C = gg_evol_hab, D = gg_ind_hab, 
           E = gg_ind_sp, design = layout)

# Save figure ----
ggsave('../Figs/figure_sup_allometries.pdf', plot = last_plot(), 
       height = 9, width = 7)
ggsave('../Figs/figure_sup_allometries.png', plot = last_plot(), 
       height = 9, width = 7)



