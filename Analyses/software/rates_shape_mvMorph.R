

# Explanation ----
# Fit models of multivariate trait evolution. We use the first two PCs 
# from the shape phyPCA, representing limb (PC1) and head (PC2) proportions.
# We fit single-regime models (BM1, OU1) on the consensus tree and on
# 100 posterior trees. 
# We fit multi-rate regime models (BMM, OUM) on 1,000 stochastic character 
# maps across the consensus tree, and on 100 stochastic character maps across 
# each of the 100 posterior trees. 
# We do this with island-mainland categories, and with the three-habitat categories 
# (habitat_broad): ground, rock and tree. 
# We do model selection based on AICc. 
# If the best-fit model is a multi-rate model (BMM or OUM), then we need to see 
# in which trait state the rates are higher (island vs mainland; ground vs rock vs tree).

# Packages ----
libs <- c('tidyverse', 'treeio', 'phytools', 'geiger', 'ggtree', 'OUwie', 
          'patchwork', 'mvMORPH', 'ggridges', 'ggpubr')
lapply(libs, require, character.only = TRUE)

parallel::detectCores()
doParallel::registerDoParallel(cores = 4)

# Import morpho ----
morpho <- read.table('objects/phypca/phypca_scores.csv', sep = ";", dec = '.', 
                     header = TRUE)
rownames(morpho) <- morpho$species
svl <- morpho$SVL
names(svl) <- rownames(morpho)
pPC1 <- morpho$PC1
names(pPC1) <- rownames(morpho)
pPC2 <- morpho$PC2
names(pPC2) <- rownames(morpho)

# Import trees ----
tree <- read.nexus('data/phylogeny/pristurus_tree_final.nex')
postrees <- readRDS('data/phylogeny/posterior_pristurus_final.rds')
name.check(tree, morpho)

# Import mapped trees ----
pd_land <- readRDS('objects/anc_rec/pd_land_cons.rds')
pd_land_post <- readRDS('objects/anc_rec/pd_land_post.rds')
# pd_land_post is a list with ntrees simmap summaries of nsim simulations each.
pd_habitat_broad <- readRDS('objects/anc_rec/pd_habitat_broad_cons.rds')
pd_habitat_broad_post <- readRDS('objects/anc_rec/pd_habitat_broad_post.rds')
pd_habitat <- readRDS('objects/anc_rec/pd_habitat_cons.rds')
pd_habitat_post <- readRDS('objects/anc_rec/pd_habitat_post.rds')

# Colors
habitat_broad_colors0 <- c(ground = 'brown', rock = 'gray', tree = 'darkgreen')
habitat_broad_colors <- c(ground = '#e0710b', rock = '#6C586E', tree = '#119616')
habitat_colors0 <- c('hard-ground' = 'brown', 'soft-ground' = 'orange', 
                     rock = 'gray', tree = 'darkgreen')
habitat_colors <- c('hard-ground' = '#D63916', 'soft-ground' = '#E9A800', 
                    rock = '#6C586E', tree = '#119616')
land_colors0 <- c(mainland = 'coral2', island = 'seagreen')
land_colors1 <- c(mainland = '#A54E29', island = '#42978A')
land_colors2 <- c(mainland = '#C56542', island = '#67D1B7')
land_colors <- c(mainland = '#C56542', island = '#42978A')

# Set ggplot2 theme
theme.clean <- function(){
  theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 0.5),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11, face = "plain"),             
          axis.title.y = element_text(size = 11, face = "plain"),  
          panel.grid = element_blank(), 
          #          panel.grid.major.x = element_blank(),                                          
          #          panel.grid.minor.x = element_blank(),
          #          panel.grid.minor.y = element_blank(),
          #          panel.grid.major.y = element_blank(),  
          axis.line.y = element_line(colour = "black", size = 0.3), 
          axis.line.x = element_line(colour = "black", size = 0.3), 
          axis.ticks.length = unit(.15, "cm"),
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 10, face = "plain"),          
          legend.position = "right")
}
theme_set(theme.clean())



# mvMorph fit models ----
# ______________________________ ----
# Land ----
# LAND CONSENSUS ----
# single-regime models: one fit per topology
# Land consensus BM1 ----
fitBM1_land_cons <- mvBM(pd_land$tree[[1]],  morpho[,7:8], 
                         model = "BM1", param=list(decomp="diagonal"))
fitBM1_land_cons$AICc

# Land consensus OU1 ----
fitOU1_land_cons <- mvOU(pd_land$tree[[1]],  morpho[,7:8], 
                         model = "OU1", param=list(decomp="diagonal"))
fitOU1_land_cons$AICc

# Land consensus EB ----
#fitEB_land_cons <- mvEB(pd_land$tree[[1]],  
#                        morpho[,7:8], param=list(decomp="diagonal"))
#fitEB_land_cons$AICc


# multirate models: one fit per character map
# # Land consensus BMM ----
ntrees <- length(pd_land$tree)
fitBMM_land_cons <- list()

for (i in 1:ntrees){
  print(paste0('Land consensus BMM -- tree ', i))
  fitBMM_land_cons[[i]] <- mvBM(pd_land$tree[[i]], morpho[,7:8], 
                                model = 'BMM', param=list(decomp='diagonal'))
}
saveRDS(fitBMM_land_cons, 'objects/mvMorph/shape_land_cons_BMM.rds')
fitBMM_land_cons <- readRDS('objects/mvMorph/shape_land_cons_BMM.rds')

# Get AICc values
fitBMM_land_cons_AICc <- c()
for (i in 1:ntrees){
  fitBMM_land_cons_AICc[i] <- fitBMM_land_cons[[i]]$AICc
}
saveRDS(fitBMM_land_cons_AICc, 'objects/mvMorph/shape_land_cons_BMM_AICc.rds')
fitBMM_land_cons_AICc <- readRDS('objects/mvMorph/shape_land_cons_BMM_AICc.rds')


hist(fitBMM_land_cons_AICc)
fitBMM_land_cons_AICc_mean <- mean(fitBMM_land_cons_AICc)

# ...Get rates BMM ----
sigmasq_island_cons_BMM_PC1 <- c()
sigmasq_island_cons_BMM_PC2 <- c()
sigmasq_mainland_cons_BMM_PC1 <- c()
sigmasq_mainland_cons_BMM_PC2 <- c()

for (i in 1:ntrees){
  print(paste0('Land consensus BMM RATES -- tree ', i))
  sigma <- fitBMM_land_cons[[i]]$sigma
  s <- diag(sigma[,,'island'])
  m <- diag(sigma[,,'mainland'])
  sigmasq_island_cons_BMM_PC1 <- c(sigmasq_island_cons_BMM_PC1, s['PC1'])
  sigmasq_island_cons_BMM_PC2 <- c(sigmasq_island_cons_BMM_PC2, s['PC2'])
  sigmasq_mainland_cons_BMM_PC1 <- c(sigmasq_mainland_cons_BMM_PC1, m['PC1'])
  sigmasq_mainland_cons_BMM_PC2 <- c(sigmasq_mainland_cons_BMM_PC2, m['PC2'])
}

# ...Plot rates BMM ----
sigmasq_island_cons_BMM_PC1_df <- data.frame(sigmasq = sigmasq_island_cons_BMM_PC1, 
                                             trait = 'land', PC ='PC1', land = 'island')
sigmasq_island_cons_BMM_PC2_df <- data.frame(sigmasq = sigmasq_island_cons_BMM_PC2, 
                                             trait = 'land', PC ='PC2', land = 'island')
sigmasq_mainland_cons_BMM_PC1_df <- data.frame(sigmasq = sigmasq_mainland_cons_BMM_PC1, 
                                               trait = 'land', PC ='PC1', land = 'mainland')
sigmasq_mainland_cons_BMM_PC2_df <- data.frame(sigmasq = sigmasq_mainland_cons_BMM_PC2, 
                                               trait = 'land', PC ='PC2', land = 'mainland')
rate_land_cons_BMM_df <- rbind(sigmasq_island_cons_BMM_PC1_df,
                               sigmasq_island_cons_BMM_PC2_df,
                               sigmasq_mainland_cons_BMM_PC1_df,
                               sigmasq_mainland_cons_BMM_PC2_df)
saveRDS(rate_land_cons_BMM_df, 'objects/mvMorph/shape_land_cons_BMM_rate_df.rds')
rate_land_cons_BMM_df <- readRDS('objects/mvMorph/shape_land_cons_BMM_rate_df.rds')
head(rate_land_cons_BMM_df)

rate_land_cons_BMM_df %>%
  filter(trait == 'PC1') %>%
  ggplot() +
  geom_density(aes(x = sigmasq, fill = land), color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = land_colors) +
  labs(title = 'PC1 rates consensus tree')

ggplot(rate_land_cons_BMM_df) +
  facet_grid(cols = vars(trait), scales = 'free_x') +
  geom_density(aes(x = sigmasq, fill = land), color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = land_colors) +
  labs(title = 'Rates of shape evolution', 
       subtitle = 'Consensus tree (1,000 character maps)') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))


# Land consensus OUM ----
ntrees <- length(pd_land$tree)
fitOUM_land_cons <- list()

for (i in 1:ntrees){
  print(paste0('Land consensus OUM -- tree ', i))
  fitOUM_land_cons[[i]] <- mvOU(pd_land$tree[[i]], morpho[,7:8], 
                                model = 'OUM', param=list(decomp='diagonal'))
}
saveRDS(fitOUM_land_cons, 'objects/mvMorph/shape_land_cons_OUM.rds')
fitOUM_land_cons <- readRDS('objects/mvMorph/shape_land_cons_OUM.rds')

# Get AICc values
fitOUM_land_cons_AICc <- c()
for (i in 1:ntrees){
  fitOUM_land_cons_AICc[i] <- fitOUM_land_cons[[i]]$AICc
}
saveRDS(fitOUM_land_cons_AICc, 'objects/mvMorph/shape_land_cons_OUM_AICc.rds')
fitOUM_land_cons_AICc <- readRDS('objects/mvMorph/shape_land_cons_OUM_AICc.rds')

hist(fitOUM_land_cons_AICc)
fitOUM_land_cons_AICc_mean <- mean(fitOUM_land_cons_AICc)



# Shape land consensus model selection ----
fitBM1_land_cons$AICc
fitOU1_land_cons$AICc
fitBMM_land_cons_AICc_mean
fitOUM_land_cons_AICc_mean
meanAICc_shape_land_cons <- data.frame(model = c('BM1', 'OU1', 'BMM', 'OUM'), 
                                       AICc = c(fitBM1_land_cons$AICc, 
                                                fitOU1_land_cons$AICc,
                                                fitBMM_land_cons_AICc_mean, 
                                                fitOUM_land_cons_AICc_mean))
saveRDS(meanAICc_shape_land_cons, 'objects/mvMorph/meanAICc_shape_land_cons.rds')
meanAICc_shape_land_cons <- readRDS('objects/mvMorph/meanAICc_shape_land_cons.rds')

meanAICc_shape_land_cons$model[meanAICc_shape_land_cons$AICc == min(meanAICc_shape_land_cons$AICc)]
# The minimum AICc value corresponds to the best fitting model: BM1

# Land consensus density plot AICc ----
fitBM1_land_cons_AICc_df <- data.frame(AICc = rep(fitBM1_land_cons$AICc, 100), 
                                       model = 'BM1', trait = 'land', tree = 'cons')
fitOU1_land_cons_AICc_df <- data.frame(AICc = rep(fitOU1_land_cons$AICc, 100), 
                                       model = 'OU1', trait = 'land', tree = 'cons')
fitBMM_land_cons_AICc_df <- data.frame(AICc = fitBMM_land_cons_AICc, 
                                       model = 'BMM', trait = 'land', tree = 'cons')
fitOUM_land_cons_AICc_df <- data.frame(AICc = fitOUM_land_cons_AICc, 
                                       model = 'OUM', trait = 'land', tree = 'cons')

aic_df_land_cons <- rbind(fitBM1_land_cons_AICc_df, fitOU1_land_cons_AICc_df,
                          fitBMM_land_cons_AICc_df, fitOUM_land_cons_AICc_df)
saveRDS(aic_df_land_cons, 'objects/mvMorph/aic_df_shape_land_cons.rds')
aic_df_land_cons <- readRDS('objects/mvMorph/aic_df_shape_land_cons.rds')


aic_plot_shape_land_cons <- ggplot(aic_df_land_cons) +
  geom_histogram(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.8) +
  scale_fill_viridis_d() +
  ggsave('plots/rates/shape/aic_plot_shape_land_cons.pdf')

# ______________________________ ----
# LAND POSTERIOR ----
# single-regime models: one fit per topology
# Land posterior BM1 ----
ntrees <- length(pd_land_post)
pd_land_post[[2]]$tree[[1]]

fitBM1_land_post <- list()

for (i in 1:ntrees){
  print(paste0('Land posterior BM1 -- tree ', i))
  fitBM1_land_post[[i]] <- mvBM(pd_land_post[[i]]$tree[[1]],  morpho[,7:8], 
                                model = "BM1", param=list(decomp="diagonal"))
}
saveRDS(fitBM1_land_post, 'objects/mvMorph/shape_land_post_BM1.rds', version = 2)
fitBM1_land_post <- readRDS('objects/mvMorph/shape_land_post_BM1.rds')

# Get AICc values
fitBM1_land_post_AICc <- c()
for (i in 1:ntrees){
  fitBM1_land_post_AICc[i] <- fitBM1_land_post[[i]]$AICc
}
saveRDS(fitBM1_land_post_AICc, 'objects/mvMorph/shape_land_post_BM1_AICc.rds', version = 2)
fitBM1_land_post_AICc <- readRDS('objects/mvMorph/shape_land_post_BM1_AICc.rds')
hist(fitBM1_land_post_AICc)
fitBM1_land_post_AICc_mean <- mean(fitBM1_land_post_AICc)

# Land posterior OU1 ----
ntrees <- length(pd_land_post)

fitOU1_land_post <- list()

for (i in 1:ntrees){
  print(paste0('Land posterior OU1 -- tree ', i))
  fitOU1_land_post[[i]] <- mvOU(pd_land_post[[i]]$tree[[1]],  morpho[,7:8], 
                                model = "OU1", param=list(decomp="diagonal"))
}
saveRDS(fitOU1_land_post, 'objects/mvMorph/shape_land_post_OU1.rds', version = 2)
fitOU1_land_post <- readRDS('objects/mvMorph/shape_land_post_OU1.rds')

# Get AICc values
fitOU1_land_post_AICc <- c()
for (i in 1:ntrees){
  fitOU1_land_post_AICc[i] <- fitOU1_land_post[[i]]$AICc
}
saveRDS(fitOU1_land_post_AICc, 'objects/mvMorph/shape_land_post_OU1_AICc.rds', version = 2)
fitOU1_land_post_AICc <- readRDS('objects/mvMorph/shape_land_post_OU1_AICc.rds')
hist(fitOU1_land_post_AICc)
fitOU1_land_post_AICc_mean <- mean(fitOU1_land_post_AICc)


# multirate models: one fit per character map
# Land posterior BMM ----
ntrees <- length(pd_land_post)
nsim <- length(pd_land_post[[1]]$tree)

fitBMM_land_post <- list()

for (i in 1:ntrees){
  fitBMM_land_post[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('Land posterior BMM -- tree ', i, '; sim ', j))
    fitBMM_land_post[[i]][[j]] <- mvBM(pd_land_post[[i]]$tree[[j]], morpho[,7:8], 
                                       model = 'BMM', param=list(decomp='diagonal'))
    
  }
}
saveRDS(fitBMM_land_post, 'objects/mvMorph/shape_land_post_BMM.rds', version = 2)
fitBMM_land_post <- readRDS('objects/mvMorph/shape_land_post_BMM.rds')

# Get AICc values
fitBMM_land_post_AICc <- c()

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('Land posterior BMM -- tree ', i, '; sim ', j))
    fitBMM_land_post_AICc <- c(fitBMM_land_post_AICc, fitBMM_land_post[[i]][[j]]$AICc)
  }
}
saveRDS(fitBMM_land_post_AICc, 'objects/mvMorph/shape_land_post_BMM_AICc.rds', version = 2)
fitBMM_land_post_AICc <- readRDS('objects/mvMorph/shape_land_post_BMM_AICc.rds')

hist(fitBMM_land_post_AICc)
fitBMM_land_post_AICc_mean <- mean(fitBMM_land_post_AICc)

# ...Get rates BMM ----
sigmasq_island_post_BMM_PC1 <- c()
sigmasq_island_post_BMM_PC2 <- c()
sigmasq_mainland_post_BMM_PC1 <- c()
sigmasq_mainland_post_BMM_PC2 <- c()

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('Land posterior BMM RATES -- tree ', i, '; nsim ', j))
    sigma <- fitBMM_land_post[[i]][[j]]$sigma
    s <- diag(sigma[,,'island'])
    m <- diag(sigma[,,'mainland'])
    sigmasq_island_post_BMM_PC1 <- c(sigmasq_island_post_BMM_PC1, s['PC1'])
    sigmasq_island_post_BMM_PC2 <- c(sigmasq_island_post_BMM_PC2, s['PC2'])
    sigmasq_mainland_post_BMM_PC1 <- c(sigmasq_mainland_post_BMM_PC1, m['PC1'])
    sigmasq_mainland_post_BMM_PC2 <- c(sigmasq_mainland_post_BMM_PC2, m['PC2'])
    
  }
}

# ...Plot rates BMM ----
sigmasq_island_post_BMM_PC1_df <- data.frame(sigmasq = sigmasq_island_post_BMM_PC1, 
                                             trait = 'land', PC ='PC1', land = 'island')
sigmasq_island_post_BMM_PC2_df <- data.frame(sigmasq = sigmasq_island_post_BMM_PC2, 
                                             trait = 'land', PC ='PC2', land = 'island')
sigmasq_mainland_post_BMM_PC1_df <- data.frame(sigmasq = sigmasq_mainland_post_BMM_PC1, 
                                               trait = 'land', PC ='PC1', land = 'mainland')
sigmasq_mainland_post_BMM_PC2_df <- data.frame(sigmasq = sigmasq_mainland_post_BMM_PC2, 
                                               trait = 'land', PC ='PC2', land = 'mainland')
rate_land_post_BMM_df <- rbind(sigmasq_island_post_BMM_PC1_df,
                               sigmasq_island_post_BMM_PC2_df,
                               sigmasq_mainland_post_BMM_PC1_df,
                               sigmasq_mainland_post_BMM_PC2_df)
saveRDS(rate_land_post_BMM_df, 'objects/mvMorph/shape_land_post_BMM_rate_df.rds', version = 2)
rate_land_post_BMM_df <- readRDS('objects/mvMorph/shape_land_post_BMM_rate_df.rds')

head(rate_land_post_BMM_df)

ggplot(rate_land_post_BMM_df) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(aes(x = sigmasq, fill = land), color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = land_colors) +
  #  xlim(0, 3e-04) +
  labs(title = 'Rates of shape evolution', 
       subtitle = 'Posterior trees (100 trees, 100 character maps)') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))
?facet_grid

# Land posterior OUM ----
ntrees <- length(pd_land_post)
nsim <- length(pd_land_post[[1]]$tree)

fitOUM_land_post <- list()

for (i in 1:ntrees){
  fitOUM_land_post[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('Land posterior OUM -- tree ', i, '; sim ', j))
    fitOUM_land_post[[i]][[j]] <- mvOU(pd_land_post[[i]]$tree[[j]], morpho[,7:8], 
                                       model = 'OUM', param=list(decomp='diagonal'))
    
  }
}
saveRDS(fitOUM_land_post, 'objects/mvMorph/shape_land_post_OUM.rds', version = 2)
fitOUM_land_post <- readRDS('objects/mvMorph/shape_land_post_OUM.rds')

# Get AICc values
fitOUM_land_post_AICc <- c()

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('Land posterior OUM -- tree ', i, '; sim ', j))
    fitOUM_land_post_AICc <- c(fitOUM_land_post_AICc, fitOUM_land_post[[i]][[j]]$AICc)
  }
}
saveRDS(fitOUM_land_post_AICc, 'objects/mvMorph/shape_land_post_OUM_AICc.rds', version = 2)
fitOUM_land_post_AICc <- readRDS('objects/mvMorph/shape_land_post_OUM_AICc.rds')


hist(fitOUM_land_post_AICc)
fitOUM_land_post_AICc_mean <- mean(fitOUM_land_post_AICc)

# Shape land posterior model selection ----
fitBM1_land_post_AICc_mean
fitOU1_land_post_AICc_mean
fitBMM_land_post_AICc_mean
fitOUM_land_post_AICc_mean
meanAICc_shape_land_post <- data.frame(model = c('BM1', 'OU1', 'BMM', 'OUM'), 
                                       AICc = c(fitBM1_land_post_AICc_mean, 
                                                fitOU1_land_post_AICc_mean,
                                                fitBMM_land_post_AICc_mean, 
                                                fitOUM_land_post_AICc_mean))
saveRDS(meanAICc_shape_land_post, 'objects/mvMorph/meanAICc_shape_land_post.rds')
meanAICc_shape_land_post <- readRDS('objects/mvMorph/meanAICc_shape_land_post.rds')

meanAICc_shape_land_post$model[meanAICc_shape_land_post$AICc == min(meanAICc_shape_land_post$AICc)]
# The minimum AICc value corresponds to the best fitting model: OU1 (and BM1)

# Land posterior density plot AICc ----
fitBM1_land_post_AICc_df <- data.frame(AICc = fitBM1_land_post_AICc, 
                                       model = 'BM1', trait = 'land', tree = 'post')
fitOU1_land_post_AICc_df <- data.frame(AICc = fitOU1_land_post_AICc, 
                                       model = 'OU1', trait = 'land', tree = 'post')
fitBMM_land_post_AICc_df <- data.frame(AICc = fitBMM_land_post_AICc, 
                                       model = 'BMM', trait = 'land', tree = 'post')
fitOUM_land_post_AICc_df <- data.frame(AICc = fitOUM_land_post_AICc, 
                                       model = 'OUM', trait = 'land', tree = 'post')

aic_df_land_post <- rbind(fitBM1_land_post_AICc_df, fitOU1_land_post_AICc_df,
                          fitBMM_land_post_AICc_df, fitOUM_land_post_AICc_df)
saveRDS(aic_df_land_post, 'objects/mvMorph/aic_df_shape_land_post.rds')
aic_df_land_post <- readRDS('objects/mvMorph/aic_df_shape_land_post.rds')

aic_plot_shape_land_post <- ggplot(aic_df_land_post) +
  geom_density(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.8) +
  scale_fill_viridis_d()

# ______________________________ ----
# _===================== ----
# ______________________________ ----
# Habitat_broad ----
# HABITAT_BROAD CONSENSUS ----
# single-regime models: one fit per topology
# habitat_broad consensus BM1 ----
fitBM1_habitat_broad_cons <- mvBM(pd_habitat_broad$tree[[1]], morpho[,7:8], 
                                  model = "BM1", param=list(decomp="diagonal"))
fitBM1_habitat_broad_cons$AICc

# habitat_broad consensus OU1 ----
fitOU1_habitat_broad_cons <- mvOU(pd_habitat_broad$tree[[1]],  morpho[,7:8], 
                                  model = "OU1", param=list(decomp="diagonal"))
fitOU1_habitat_broad_cons$AICc

# habitat_broad consensus EB ----
#fitEB_habitat_broad_cons <- mvEB(pd_habitat_broad$tree[[1]],  
#                        morpho[,7:8], param=list(decomp="diagonal"))
#fitEB_habitat_broad_cons$AICc


# multirate models: one fit per character map
# # habitat_broad consensus BMM ----
ntrees <- length(pd_habitat_broad$tree)
fitBMM_habitat_broad_cons <- list()

for (i in 1:ntrees){
  print(paste0('habitat_broad consensus BMM -- tree ', i))
  fitBMM_habitat_broad_cons[[i]] <- mvBM(pd_habitat_broad$tree[[i]], morpho[,7:8], 
                                         model = 'BMM', param=list(decomp='diagonal'))
}
saveRDS(fitBMM_habitat_broad_cons, 'objects/mvMorph/shape_habitat_broad_cons_BMM.rds')
fitBMM_habitat_broad_cons <- readRDS('objects/mvMorph/shape_habitat_broad_cons_BMM.rds')

# Get AICc values
fitBMM_habitat_broad_cons_AICc <- c()
for (i in 1:ntrees){
  fitBMM_habitat_broad_cons_AICc[i] <- fitBMM_habitat_broad_cons[[i]]$AICc
}
saveRDS(fitBMM_habitat_broad_cons_AICc, 'objects/mvMorph/shape_habitat_broad_cons_BMM_AICc.rds')
fitBMM_habitat_broad_cons_AICc <- readRDS('objects/mvMorph/shape_habitat_broad_cons_BMM_AICc.rds')


hist(fitBMM_habitat_broad_cons_AICc)
fitBMM_habitat_broad_cons_AICc_mean <- mean(fitBMM_habitat_broad_cons_AICc)

# ...Get rates BMM ----
sigmasq_ground_cons_BMM_PC1 <- c()
sigmasq_ground_cons_BMM_PC2 <- c()
sigmasq_rock_cons_BMM_PC1 <- c()
sigmasq_rock_cons_BMM_PC2 <- c()
sigmasq_tree_cons_BMM_PC1 <- c()
sigmasq_tree_cons_BMM_PC2 <- c()

for (i in 1:ntrees){
  print(paste0('habitat_broad consensus BMM RATES -- tree ', i))
  sigma <- fitBMM_habitat_broad_cons[[i]]$sigma
  g <- diag(sigma[,,'ground'])
  r <- diag(sigma[,,'rock'])
  t <- diag(sigma[,,'tree'])
  sigmasq_ground_cons_BMM_PC1 <- c(sigmasq_ground_cons_BMM_PC1, g['PC1'])
  sigmasq_ground_cons_BMM_PC2 <- c(sigmasq_ground_cons_BMM_PC2, g['PC2'])
  sigmasq_rock_cons_BMM_PC1 <- c(sigmasq_rock_cons_BMM_PC1, r['PC1'])
  sigmasq_rock_cons_BMM_PC2 <- c(sigmasq_rock_cons_BMM_PC2, r['PC2'])
  sigmasq_tree_cons_BMM_PC1 <- c(sigmasq_tree_cons_BMM_PC1, t['PC1'])
  sigmasq_tree_cons_BMM_PC2 <- c(sigmasq_tree_cons_BMM_PC2, t['PC2'])
}

# ...Plot rates BMM ----
sigmasq_ground_cons_BMM_PC1_df <- data.frame(sigmasq = sigmasq_ground_cons_BMM_PC1, 
                                             trait = 'habitat', PC = 'PC1', habitat_broad = 'ground')
sigmasq_ground_cons_BMM_PC2_df <- data.frame(sigmasq = sigmasq_ground_cons_BMM_PC2, 
                                             trait = 'habitat', PC = 'PC2', habitat_broad = 'ground')
sigmasq_rock_cons_BMM_PC1_df <- data.frame(sigmasq = sigmasq_rock_cons_BMM_PC1, 
                                           trait = 'habitat', PC = 'PC1', habitat_broad = 'rock')
sigmasq_rock_cons_BMM_PC2_df <- data.frame(sigmasq = sigmasq_rock_cons_BMM_PC2, 
                                           trait = 'habitat', PC = 'PC2', habitat_broad = 'rock')
sigmasq_tree_cons_BMM_PC1_df <- data.frame(sigmasq = sigmasq_tree_cons_BMM_PC1, 
                                           trait = 'habitat', PC = 'PC1', habitat_broad = 'tree')
sigmasq_tree_cons_BMM_PC2_df <- data.frame(sigmasq = sigmasq_tree_cons_BMM_PC2, 
                                           trait = 'habitat', PC = 'PC2', habitat_broad = 'tree')
rate_habitat_broad_cons_BMM_df <- rbind(sigmasq_ground_cons_BMM_PC1_df,
                                        sigmasq_ground_cons_BMM_PC2_df,
                                        sigmasq_rock_cons_BMM_PC1_df,
                                        sigmasq_rock_cons_BMM_PC2_df, 
                                        sigmasq_tree_cons_BMM_PC1_df, 
                                        sigmasq_tree_cons_BMM_PC2_df)
saveRDS(rate_habitat_broad_cons_BMM_df, 'objects/mvMorph/shape_habitat_broad_cons_BMM_rate_df.rds')
rate_habitat_broad_cons_BMM_df <- readRDS('objects/mvMorph/shape_habitat_broad_cons_BMM_rate_df.rds')


rate_habitat_broad_cons_BMM_df %>%
  filter(trait == 'PC1') %>%
  ggplot() +
  geom_density(aes(x = sigmasq, fill = habitat_broad), color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = habitat_broad_colors) +
  labs(title = 'PC1 rates consensus tree')

rate_plot_shape_habitat_cons <- ggplot(rate_habitat_broad_cons_BMM_df) +
  facet_grid(rows = vars(trait), scales = 'free', space = 'fixed') +
  geom_density(aes(x = sigmasq, fill = habitat_broad), color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = habitat_broad_colors) +
  #  xlim(0, 3e-04) +
  labs(title = 'Rates of shape evolution', 
       subtitle = 'Consensus tree (1,000 character maps)') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5)) +
  ggsave('plots/rates/shape/rate_plot_shape_habitat_cons.pdf')

# habitat_broad consensus OUM ----
ntrees <- length(pd_habitat_broad$tree)
fitOUM_habitat_broad_cons <- list()

for (i in 1:ntrees){
  print(paste0('habitat_broad consensus OUM -- tree ', i))
  fitOUM_habitat_broad_cons[[i]] <- mvOU(pd_habitat_broad$tree[[i]], morpho[,7:8], 
                                         model = 'OUM', param=list(decomp='diagonal'))
}
saveRDS(fitOUM_habitat_broad_cons, 'objects/mvMorph/shape_habitat_broad_cons_OUM.rds')
fitOUM_habitat_broad_cons <- readRDS('objects/mvMorph/shape_habitat_broad_cons_OUM.rds')

# Get AICc values
fitOUM_habitat_broad_cons_AICc <- c()
for (i in 1:ntrees){
  fitOUM_habitat_broad_cons_AICc[i] <- fitOUM_habitat_broad_cons[[i]]$AICc
}
saveRDS(fitOUM_habitat_broad_cons_AICc, 'objects/mvMorph/shape_habitat_broad_cons_OUM_AICc.rds')
fitOUM_habitat_broad_cons_AICc <- readRDS('objects/mvMorph/shape_habitat_broad_cons_OUM_AICc.rds')
hist(fitOUM_habitat_broad_cons_AICc)
fitOUM_habitat_broad_cons_AICc_mean <- mean(fitOUM_habitat_broad_cons_AICc)

# Shape habitat_broad consensus model selection ----
fitBM1_habitat_broad_cons$AICc
fitOU1_habitat_broad_cons$AICc
fitBMM_habitat_broad_cons_AICc_mean
fitOUM_habitat_broad_cons_AICc_mean
meanAICc_shape_habitat_broad_cons <- data.frame(model = c('BM1', 'OU1', 'BMM', 'OUM'), 
                                                AICc = c(fitBM1_habitat_broad_cons$AICc, 
                                                         fitOU1_habitat_broad_cons$AICc,
                                                         fitBMM_habitat_broad_cons_AICc_mean, 
                                                         fitOUM_habitat_broad_cons_AICc_mean))
saveRDS(meanAICc_shape_habitat_broad_cons, 'objects/mvMorph/meanAICc_shape_habitat_broad_cons.rds')
meanAICc_shape_habitat_broad_cons <- readRDS('objects/mvMorph/meanAICc_shape_habitat_broad_cons.rds')

meanAICc_shape_habitat_broad_cons$model[meanAICc_shape_habitat_broad_cons$AICc == min(meanAICc_shape_habitat_broad_cons$AICc)]
# The minimum AICc value corresponds to the best fitting model: BM1

# habitat_broad consensus density plot AICc ----
fitBM1_habitat_broad_cons_AICc_df <- data.frame(AICc = rep(fitBM1_habitat_broad_cons$AICc, 100), 
                                                model = 'BM1', trait = 'habitat', tree = 'cons')
fitOU1_habitat_broad_cons_AICc_df <- data.frame(AICc = rep(fitOU1_habitat_broad_cons$AICc, 100), 
                                                model = 'OU1', trait = 'habitat', tree = 'cons')
fitBMM_habitat_broad_cons_AICc_df <- data.frame(AICc = fitBMM_habitat_broad_cons_AICc, 
                                                model = 'BMM', trait = 'habitat', tree = 'cons')
fitOUM_habitat_broad_cons_AICc_df <- data.frame(AICc = fitOUM_habitat_broad_cons_AICc, 
                                                model = 'OUM', trait = 'habitat', tree = 'cons')

aic_df_habitat_cons <- rbind(fitBM1_habitat_broad_cons_AICc_df, fitOU1_habitat_broad_cons_AICc_df,
                             fitBMM_habitat_broad_cons_AICc_df, fitOUM_habitat_broad_cons_AICc_df)
saveRDS(aic_df_habitat_cons, 'objects/mvMorph/aic_df_habitat_cons.rds')
aic_df_habitat_cons <- readRDS('objects/mvMorph/aic_df_habitat_cons.rds')


aic_plot_shape_habitat_broad_cons <- ggplot(aic_df_habitat_cons) +
  geom_histogram(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.8) +
  scale_fill_viridis_d()


# ____________________________ ----
# HABITAT_BROAD POSTERIOR ----
# single-regime models: one fit per topology
# habitat_broad posterior BM1 ----
ntrees <- length(pd_habitat_broad_post)

fitBM1_habitat_broad_post <- list()

for (i in 1:ntrees){
  print(paste0('habitat_broad posterior BM1 -- tree ', i))
  fitBM1_habitat_broad_post[[i]] <- mvBM(pd_habitat_broad_post[[i]]$tree[[1]],  morpho[,7:8], 
                                         model = "BM1", param=list(decomp="diagonal"))
}
saveRDS(fitBM1_habitat_broad_post, 'objects/mvMorph/shape_habitat_broad_post_BM1.rds', version = 2)
fitBM1_habitat_broad_post <- readRDS('objects/mvMorph/shape_habitat_broad_post_BM1.rds')

# Get AICc values
fitBM1_habitat_broad_post_AICc <- c()
for (i in 1:ntrees){
  fitBM1_habitat_broad_post_AICc[i] <- fitBM1_habitat_broad_post[[i]]$AICc
}
saveRDS(fitBM1_habitat_broad_post_AICc, 'objects/mvMorph/shape_habitat_broad_post_BM1_AICc.rds', version = 2)
fitBM1_habitat_broad_post_AICc <- readRDS('objects/mvMorph/shape_habitat_broad_post_BM1_AICc.rds')
hist(fitBM1_habitat_broad_post_AICc)
fitBM1_habitat_broad_post_AICc_mean <- mean(fitBM1_habitat_broad_post_AICc)

# habitat_broad posterior OU1 ----
ntrees <- length(pd_habitat_broad_post)

fitOU1_habitat_broad_post <- list()

for (i in 1:ntrees){
  print(paste0('habitat_broad posterior OU1 -- tree ', i))
  fitOU1_habitat_broad_post[[i]] <- mvOU(pd_habitat_broad_post[[i]]$tree[[1]],  morpho[,7:8], 
                                         model = "OU1", param=list(decomp="diagonal"))
}
saveRDS(fitOU1_habitat_broad_post, 'objects/mvMorph/shape_habitat_broad_post_OU1.rds', version = 2)
fitOU1_habitat_broad_post <- readRDS('objects/mvMorph/shape_habitat_broad_post_OU1.rds')

# Get AICc values
fitOU1_habitat_broad_post_AICc <- c()
for (i in 1:ntrees){
  fitOU1_habitat_broad_post_AICc[i] <- fitOU1_habitat_broad_post[[i]]$AICc
}
saveRDS(fitOU1_habitat_broad_post_AICc, 'objects/mvMorph/shape_habitat_broad_post_OU1_AICc.rds', version = 2)
fitOU1_habitat_broad_post_AICc <- readRDS('objects/mvMorph/shape_habitat_broad_post_OU1_AICc.rds')
hist(fitOU1_habitat_broad_post_AICc)
fitOU1_habitat_broad_post_AICc_mean <- mean(fitOU1_habitat_broad_post_AICc)


# multirate models: one fit per character map
# habitat_broad posterior BMM ----
ntrees <- length(pd_habitat_broad_post)
nsim <- length(pd_habitat_broad_post[[1]]$tree)

fitBMM_habitat_broad_post <- list()

for (i in 1:ntrees){
  fitBMM_habitat_broad_post[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('habitat_broad posterior BMM -- tree ', i, '; sim ', j))
    fitBMM_habitat_broad_post[[i]][[j]] <- mvBM(pd_habitat_broad_post[[i]]$tree[[j]], morpho[,7:8], 
                                                model = 'BMM', param=list(decomp='diagonal'))
    
  }
}
saveRDS(fitBMM_habitat_broad_post, 'objects/mvMorph/shape_habitat_broad_post_BMM.rds', version = 2)
fitBMM_habitat_broad_post <- readRDS('objects/mvMorph/shape_habitat_broad_post_BMM.rds')

# Get AICc values
fitBMM_habitat_broad_post_AICc <- c()

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('habitat_broad posterior BMM -- tree ', i, '; sim ', j))
    fitBMM_habitat_broad_post_AICc <- c(fitBMM_habitat_broad_post_AICc, fitBMM_habitat_broad_post[[i]][[j]]$AICc)
  }
}
saveRDS(fitBMM_habitat_broad_post_AICc, 'objects/mvMorph/shape_habitat_broad_post_BMM_AICc.rds', version = 2)
fitBMM_habitat_broad_post_AICc <- readRDS('objects/mvMorph/shape_habitat_broad_post_BMM_AICc.rds')


hist(fitBMM_habitat_broad_post_AICc)
fitBMM_habitat_broad_post_AICc_mean <- mean(fitBMM_habitat_broad_post_AICc)

# ...Get rates BMM ----
sigmasq_ground_post_BMM_PC1 <- c()
sigmasq_ground_post_BMM_PC2 <- c()
sigmasq_rock_post_BMM_PC1 <- c()
sigmasq_rock_post_BMM_PC2 <- c()
sigmasq_tree_post_BMM_PC1 <- c()
sigmasq_tree_post_BMM_PC2 <- c()


for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('habitat_broad posterior BMM RATES -- tree ', i, '; nsim ', j))
    sigma <- fitBMM_habitat_broad_post[[i]][[j]]$sigma
    g <- diag(sigma[,,'ground'])
    r <- diag(sigma[,,'rock'])
    t <- diag(sigma[,,'tree'])
    sigmasq_ground_post_BMM_PC1 <- c(sigmasq_ground_post_BMM_PC1, g['PC1'])
    sigmasq_ground_post_BMM_PC2 <- c(sigmasq_ground_post_BMM_PC2, g['PC2'])
    sigmasq_rock_post_BMM_PC1 <- c(sigmasq_rock_post_BMM_PC1, r['PC1'])
    sigmasq_rock_post_BMM_PC2 <- c(sigmasq_rock_post_BMM_PC2, r['PC2'])
    sigmasq_tree_post_BMM_PC1 <- c(sigmasq_tree_post_BMM_PC1, t['PC1'])
    sigmasq_tree_post_BMM_PC2 <- c(sigmasq_tree_post_BMM_PC2, t['PC2'])
  }
}

# ...Plot rates BMM ----
sigmasq_ground_post_BMM_PC1_df <- data.frame(sigmasq = sigmasq_ground_post_BMM_PC1, 
                                             trait = 'habitat', PC = 'PC1', habitat_broad = 'ground')
sigmasq_ground_post_BMM_PC2_df <- data.frame(sigmasq = sigmasq_ground_post_BMM_PC2, 
                                             trait = 'habitat', PC = 'PC2', habitat_broad = 'ground')
sigmasq_rock_post_BMM_PC1_df <- data.frame(sigmasq = sigmasq_rock_post_BMM_PC1, 
                                           trait = 'habitat', PC = 'PC1', habitat_broad = 'rock')
sigmasq_rock_post_BMM_PC2_df <- data.frame(sigmasq = sigmasq_rock_post_BMM_PC2, 
                                           trait = 'habitat', PC = 'PC2', habitat_broad = 'rock')
sigmasq_tree_post_BMM_PC1_df <- data.frame(sigmasq = sigmasq_tree_post_BMM_PC1, 
                                           trait = 'habitat', PC = 'PC1', habitat_broad = 'tree')
sigmasq_tree_post_BMM_PC2_df <- data.frame(sigmasq = sigmasq_tree_post_BMM_PC2, 
                                           trait = 'habitat', PC = 'PC2', habitat_broad = 'tree')

rate_habitat_broad_post_BMM_df <- rbind(sigmasq_ground_post_BMM_PC1_df,
                                        sigmasq_ground_post_BMM_PC2_df,
                                        sigmasq_rock_post_BMM_PC1_df,
                                        sigmasq_rock_post_BMM_PC2_df, 
                                        sigmasq_tree_post_BMM_PC1_df,
                                        sigmasq_tree_post_BMM_PC2_df)
saveRDS(rate_habitat_broad_post_BMM_df, 'objects/mvMorph/shape_habitat_broad_post_BMM_rate_df.rds', version = 2)
rate_habitat_broad_post_BMM_df <- readRDS('objects/mvMorph/shape_habitat_broad_post_BMM_rate_df.rds')


head(rate_habitat_broad_post_BMM_df)

ggplot(rate_habitat_broad_post_BMM_df) +
  facet_wrap(facets = vars(PC), scales = 'free') +
  geom_density(aes(x = sigmasq, fill = habitat_broad), 
               color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = habitat_broad_colors) +
  #  xlim(0, 3e-04) +
  labs(title = 'Rates of shape evolution', 
       subtitle = 'Posterior trees (100 trees, 100 character maps)') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))
?facet_grid
?facet_wrap

# habitat_broad posterior OUM ----
ntrees <- length(pd_habitat_broad_post)
nsim <- length(pd_habitat_broad_post[[1]]$tree)

fitOUM_habitat_broad_post <- list()

for (i in 1:ntrees){
  fitOUM_habitat_broad_post[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('habitat_broad posterior OUM -- tree ', i, '; sim ', j))
    fitOUM_habitat_broad_post[[i]][[j]] <- mvOU(pd_habitat_broad_post[[i]]$tree[[j]], morpho[,7:8], 
                                                model = 'OUM', param=list(decomp='diagonal'))
    
  }
}
saveRDS(fitOUM_habitat_broad_post, 'objects/mvMorph/shape_habitat_broad_post_OUM.rds', version = 2)
fitOUM_habitat_broad_post <- readRDS('objects/mvMorph/shape_habitat_broad_post_OUM.rds')

# Get AICc values
fitOUM_habitat_broad_post_AICc <- c()

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('habitat_broad posterior OUM -- tree ', i, '; sim ', j))
    fitOUM_habitat_broad_post_AICc <- c(fitOUM_habitat_broad_post_AICc, fitOUM_habitat_broad_post[[i]][[j]]$AICc)
  }
}
saveRDS(fitOUM_habitat_broad_post_AICc, 'objects/mvMorph/shape_habitat_broad_post_OUM_AICc.rds', version = 2)
fitOUM_habitat_broad_post_AICc <- readRDS('objects/mvMorph/shape_habitat_broad_post_OUM_AICc.rds')
hist(fitOUM_habitat_broad_post_AICc)
fitOUM_habitat_broad_post_AICc_mean <- mean(fitOUM_habitat_broad_post_AICc)


# Shape habitat_broad posterior model selection ----
fitBM1_habitat_broad_post_AICc_mean
fitOU1_habitat_broad_post_AICc_mean
fitBMM_habitat_broad_post_AICc_mean
fitOUM_habitat_broad_post_AICc_mean
meanAICc_shape_habitat_broad_post <- data.frame(model = c('BM1', 'OU1', 'BMM', 'OUM'), 
                                                AICc = c(fitBM1_habitat_broad_post_AICc_mean, 
                                                         fitOU1_habitat_broad_post_AICc_mean,
                                                         fitBMM_habitat_broad_post_AICc_mean, 
                                                         fitOUM_habitat_broad_post_AICc_mean))
saveRDS(meanAICc_shape_habitat_broad_post, 'objects/mvMorph/meanAICc_shape_habitat_broad_post.rds')
meanAICc_shape_habitat_broad_post <- readRDS('objects/mvMorph/meanAICc_shape_habitat_broad_post.rds')
meanAICc_shape_habitat_broad_post$model[meanAICc_shape_habitat_broad_post$AICc == min(meanAICc_shape_habitat_broad_post$AICc)]
# The minimum AICc value corresponds to the best fitting model: OU1 (and BM1)

# habitat_broad posterior density plot AICc ----
fitBM1_habitat_broad_post_AICc_df <- data.frame(AICc = fitBM1_habitat_broad_post_AICc, 
                                                model = 'BM1', trait = 'habitat', tree = 'post')
fitOU1_habitat_broad_post_AICc_df <- data.frame(AICc = fitOU1_habitat_broad_post_AICc, 
                                                model = 'OU1', trait = 'habitat', tree = 'post')
fitBMM_habitat_broad_post_AICc_df <- data.frame(AICc = fitBMM_habitat_broad_post_AICc, 
                                                model = 'BMM', trait = 'habitat', tree = 'post')
fitOUM_habitat_broad_post_AICc_df <- data.frame(AICc = fitOUM_habitat_broad_post_AICc, 
                                                model = 'OUM', trait = 'habitat', tree = 'post')

aic_df_habitat_post <- rbind(fitBM1_habitat_broad_post_AICc_df, fitOU1_habitat_broad_post_AICc_df,
                             fitBMM_habitat_broad_post_AICc_df, fitOUM_habitat_broad_post_AICc_df)
saveRDS(aic_df_habitat_post, 'objects/mvMorph/aic_df_habitat_post.rds')
aic_df_habitat_post <- readRDS('objects/mvMorph/aic_df_habitat_post.rds')

aic_plot_shape_habitat_broad_post <- ggplot(aic_df_habitat_post) +
  geom_density(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.8) +
  scale_fill_viridis_d()






