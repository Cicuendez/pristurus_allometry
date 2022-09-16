# Rates of morphological evolution ----
# 

version

libs <- c('tidyverse', 'treeio', 'phytools', 'geiger', 'ggtree', 'OUwie', 
          'patchwork', 'mvMORPH')
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

saveRDS(pd_land, 'objects/anc_rec/pd_land_cons.rds', version = 2)
saveRDS(pd_land_post, 'objects/anc_rec/pd_land_post.rds', version = 2)
saveRDS(pd_habitat_broad, 'objects/anc_rec/pd_habitat_broad_cons.rds', version = 2)
saveRDS(pd_habitat_broad_post, 'objects/anc_rec/pd_habitat_broad_post.rds', version = 2)
saveRDS(pd_habitat, 'objects/anc_rec/pd_habitat_cons.rds', version = 2)
saveRDS(pd_habitat_post, 'objects/anc_rec/pd_habitat_post.rds', version = 2)



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
          plot.title = element_text(size = 14, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 10, face = "plain"),          
          legend.position = "right")
}
theme_set(theme.clean())


# OUwie fit models ----
?OUwie
# We fit three models: BM1, BMS, and OU1. 
# We don't fit the rest of models in OUwie 

# 1. SVL ----
# 1.1. SVL land ----
rownames(morpho)
colnames(morpho)
svl_land_df <- morpho %>%
  select(species, land, SVL)

# 1.1.1. SVL land consensus ----
length(pd_land$tree)

# SVL land consensus BM1 ----
svl_land_cons_BM1 <- list()
for (i in 1:length(pd_land$tree)){
  print(i)
  svl_land_cons_BM1[[i]] <- OUwie(pd_land$tree[[i]], svl_land_df, model="BM1", 
                            simmap.tree=TRUE, root.station=FALSE, 
                            algorithm = "three.point", scaleHeight = TRUE)
}
saveRDS(svl_land_cons_BM1, 'objects/OUwie/svl_land_cons_BM1.rds', version = 2)
svl_land_cons_BM1 <- readRDS('objects/OUwie/svl_land_cons_BM1.rds')

# Get AICc values
svl_land_cons_BM1_AICc <- c()
for (i in 1:length(pd_land$tree)){
  svl_land_cons_BM1_AICc[i] <- svl_land_cons_BM1[[i]]$AICc
}

hist(svl_land_cons_BM1_AICc)
svl_land_cons_BM1_AICc_mean <- mean(svl_land_cons_BM1_AICc)

# SVL land consensus BMS ----
svl_land_cons_BMS <- list()
for (i in 1:length(pd_land$tree)){
  print(i)
  svl_land_cons_BMS[[i]] <- OUwie(pd_land$tree[[i]], svl_land_df, model="BMS", 
                                  simmap.tree=TRUE, root.station=FALSE, 
                                  algorithm = "three.point", scaleHeight = TRUE)
}
saveRDS(svl_land_cons_BMS, 'objects/OUwie/svl_land_cons_BMS.rds', version = 2)
svl_land_cons_BMS <- readRDS('objects/OUwie/svl_land_cons_BMS.rds')

svl_land_cons_BMS[[2]]

# Get AICc values
svl_land_cons_BMS_AICc <- c()
for (i in 1:length(pd_land$tree)){
  svl_land_cons_BMS_AICc[i] <- svl_land_cons_BMS[[i]]$AICc
}

hist(svl_land_cons_BMS_AICc)
svl_land_cons_BMS_AICc_mean <- mean(svl_land_cons_BMS_AICc)

# ...Plot rates BMS (sigmasq) ----
sigmasq_island_cons_BMS <- c()
sigmasq_mainland_cons_BMS <- c()
for (i in 1:length(svl_land_cons_BMS)){
  #print(i)
  solution <- as.data.frame(svl_land_cons_BMS[[i]]$solution)
  s <- solution['sigma.sq', 'island']
  m <- solution['sigma.sq', 'mainland']
  sigmasq_island_cons_BMS <- c(sigmasq_island_cons_BMS, s)
  sigmasq_mainland_cons_BMS <- c(sigmasq_mainland_cons_BMS, m)
}

sigmasq_island_cons_BMS_df <- data.frame(sigmasq = sigmasq_island_cons_BMS, 
                                         trait = 'land', state = 'island')
sigmasq_mainland_cons_BMS_df <- data.frame(sigmasq = sigmasq_mainland_cons_BMS, 
                                           trait = 'land', state = 'mainland')
sigmasq_land_cons_BMS_df <- rbind(sigmasq_island_cons_BMS_df, sigmasq_mainland_cons_BMS_df)

saveRDS(sigmasq_land_cons_BMS_df, 'objects/OUwie/svl_land_cons_BMS_rate_df.rds', version = 2)
sigmasq_land_cons_BMS_df <- readRDS('objects/OUwie/svl_land_cons_BMS_rate_df.rds')

rate_plot_land_cons <- ggplot(sigmasq_land_cons_BMS_df) +
  geom_density(aes(x = sigmasq, fill = land), color = 'transparent') +
  scale_fill_manual(values = land_colors) +
  labs(title = 'SVL rates consensus tree')

# SVL land consensus OU1 ----
svl_land_cons_OU1 <- list()
for (i in 1:length(pd_land$tree)){
  print(i)
  svl_land_cons_OU1[[i]] <- OUwie(pd_land$tree[[i]], svl_land_df, model="OU1", 
                                  simmap.tree=TRUE, root.station=TRUE, 
                                  algorithm = "three.point", scaleHeight = TRUE)
}
saveRDS(svl_land_cons_OU1, 'objects/OUwie/svl_land_cons_OU1.rds', version = 2)
svl_land_cons_OU1 <- readRDS('objects/OUwie/svl_land_cons_OU1.rds')

# Get AICc values
svl_land_cons_OU1_AICc <- c()
for (i in 1:length(pd_land$tree)){
  svl_land_cons_OU1_AICc[i] <- svl_land_cons_OU1[[i]]$AICc
}

hist(svl_land_cons_OU1_AICc)
svl_land_cons_OU1_AICc_mean <- mean(svl_land_cons_OU1_AICc)


# SVL land consensus model selection ----
svl_land_cons_BM1_AICc_mean
svl_land_cons_BMS_AICc_mean
svl_land_cons_OU1_AICc_mean
meanAICc_svl_land_cons <- data.frame(model = c('BM1', 'BMS', 'OU1'), 
           AICc = c(svl_land_cons_BM1_AICc_mean, 
                    svl_land_cons_BMS_AICc_mean,
                    svl_land_cons_OU1_AICc_mean))
saveRDS(meanAICc_svl_land_cons, 'objects/OUwie/meanAICc_svl_land_cons.rds')
meanAICc_svl_land_cons <- readRDS('objects/OUwie/meanAICc_svl_land_cons.rds')

meanAICc_svl_land_cons$model[meanAICc_svl_land_cons$AICc == min(meanAICc_svl_land_cons$AICc)]
# The minimum AICc value corresponds to the best fitting model: BMS

# SVL land consensus density plot AICc ----
svl_land_cons_BM1_AICc_df <- data.frame(AICc = svl_land_cons_BM1_AICc, model = 'BM1', 
                                        trait = 'land', tree = 'cons')
svl_land_cons_BMS_AICc_df <- data.frame(AICc = svl_land_cons_BMS_AICc, model = 'BMS', 
                                        trait = 'land', tree = 'cons')
svl_land_cons_OU1_AICc_df <- data.frame(AICc = svl_land_cons_OU1_AICc, model = 'OU1', 
                                        trait = 'land', tree = 'cons')

aic_df_svl_land_cons <- rbind(svl_land_cons_BM1_AICc_df, svl_land_cons_BMS_AICc_df, svl_land_cons_OU1_AICc_df)
head(aic_df_svl_land_cons)
saveRDS(aic_df_svl_land_cons, 'objects/OUwie/aic_df_svl_land_cons.rds')
aic_df_svl_land_cons <- readRDS('objects/OUwie/aic_df_svl_land_cons.rds')

aic_plot_svl_land_cons <- ggplot(aic_df_svl_land_cons) +
  geom_histogram(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.7) +
  scale_fill_viridis_d()

aic_df_svl_land_cons %>%
  filter(model %in% c('BMS', 'OU1')) %>%
  ggplot() +
  geom_histogram(mapping = aes(x = AICc, fill = model), color = 'transparent') +
  scale_fill_viridis_d()

# 1.1.2. SVL land posterior ----
# 100 trees, with 100 simulations each. 
ntrees <- 100
nsim <- 100

# SVL land posterior BM1 ----
length(pd_land_post)
length(pd_land_post[[1]]$tree)
pd_land_post[[1]]$tree[[1]]


svl_land_post_BM1 <- list()
for (i in 1:ntrees){
  svl_land_post_BM1[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_land_post_BM1[[i]][[j]] <- OUwie(pd_land_post[[i]]$tree[[j]], svl_land_df, 
                                    model="BM1", 
                                    simmap.tree=TRUE, root.station=FALSE, 
                                    algorithm = "three.point", scaleHeight = TRUE)
    
  }
}
saveRDS(svl_land_post_BM1, 'objects/OUwie/svl_land_post_BM1.rds', version = 2)
svl_land_post_BM1 <- readRDS('objects/OUwie/svl_land_post_BM1.rds')

# Get AICc values
svl_land_post_BM1_AICc <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_land_post_BM1_AICc <- c(svl_land_post_BM1_AICc, svl_land_post_BM1[[i]][[j]]$AICc)
  }
}


hist(svl_land_post_BM1_AICc)
svl_land_post_BM1_AICc_mean <- mean(svl_land_post_BM1_AICc)

# SVL land posterior BMS ----
svl_land_post_BMS <- list()
for (i in 1:ntrees){
  svl_land_post_BMS[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_land_post_BMS[[i]][[j]] <- OUwie(pd_land_post[[i]]$tree[[j]], svl_land_df, 
                                         model="BMS", 
                                         simmap.tree=TRUE, root.station=TRUE, 
                                         algorithm = "three.point", scaleHeight = TRUE)
    
  }
}
saveRDS(svl_land_post_BMS, 'objects/OUwie/svl_land_post_BMS.rds', version = 2)
svl_land_post_BMS <- readRDS('objects/OUwie/svl_land_post_BMS.rds')

# Get AICc values
svl_land_post_BMS_AICc <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_land_post_BMS_AICc <- c(svl_land_post_BMS_AICc, svl_land_post_BMS[[i]][[j]]$AICc)
  }
}

hist(svl_land_post_BMS_AICc)
svl_land_post_BMS_AICc_mean <- mean(svl_land_post_BMS_AICc)

# ...Plot rates BMS (sigmasq) ----
sigmasq_island_post_BMS <- c()
sigmasq_mainland_post_BMS <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    solution <- as.data.frame(svl_land_post_BMS[[i]][[j]]$solution)
    s <- solution['sigma.sq', 'island']
    m <- solution['sigma.sq', 'mainland']
    sigmasq_island_post_BMS <- c(sigmasq_island_post_BMS, s)
    sigmasq_mainland_post_BMS <- c(sigmasq_mainland_post_BMS, m)
  }
}

sigmasq_island_post_BMS_df <- data.frame(sigmasq = sigmasq_island_post_BMS, 
                                         trait = 'land', state = 'island')
sigmasq_mainland_post_BMS_df <- data.frame(sigmasq = sigmasq_mainland_post_BMS, 
                                           trait = 'land', state = 'mainland')
sigmasq_land_post_BMS_df <- rbind(sigmasq_island_post_BMS_df, sigmasq_mainland_post_BMS_df)

saveRDS(sigmasq_land_post_BMS_df, 'objects/OUwie/svl_land_post_BMS_rate_df.rds')
size_rate_land_post_BMS_df <- readRDS('objects/OUwie/svl_land_post_BMS_rate_df.rds')

rate_plot_land_post <- ggplot(size_rate_land_post_BMS_df) +
  geom_density(aes(x = sigmasq, fill = land), color = 'transparent') +
  scale_fill_manual(values = land_colors) +
  labs(title = 'SVL rates posterior trees')


# SVL land posterior OU1 ----
svl_land_post_OU1 <- list()
for (i in 1:ntrees){
  svl_land_post_OU1[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_land_post_OU1[[i]][[j]] <- OUwie(pd_land_post[[i]]$tree[[j]], svl_land_df, 
                                         model="OU1", 
                                         simmap.tree=TRUE, root.station=TRUE, 
                                         algorithm = "three.point", scaleHeight = TRUE)
    
  }
}
saveRDS(svl_land_post_OU1, 'objects/OUwie/svl_land_post_OU1.rds', version = 2)
svl_land_post_OU1 <- readRDS('objects/OUwie/svl_land_post_OU1.rds')

# Get AICc values
svl_land_post_OU1_AICc <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_land_post_OU1_AICc <- c(svl_land_post_OU1_AICc, svl_land_post_OU1[[i]][[j]]$AICc)
  }
}

hist(svl_land_post_OU1_AICc)
svl_land_post_OU1_AICc_mean <- mean(svl_land_post_OU1_AICc)

# SVL land consensus model selection ----
svl_land_post_BM1_AICc_mean
svl_land_post_BMS_AICc_mean
svl_land_post_OU1_AICc_mean
meanAICc_svl_land_post <- data.frame(model = c('BM1', 'BMS', 'OU1'), 
                                     AICc = c(svl_land_post_BM1_AICc_mean, 
                                              svl_land_post_BMS_AICc_mean,
                                              svl_land_post_OU1_AICc_mean))
saveRDS(meanAICc_svl_land_post, 'objects/OUwie/meanAICc_svl_land_post.rds')
meanAICc_svl_land_post <- readRDS('objects/OUwie/meanAICc_svl_land_post.rds')
meanAICc_svl_land_post$model[meanAICc_svl_land_post$AICc == min(meanAICc_svl_land_post$AICc)]
# The minimum AICc value corresponds to the best fitting model: 

# SVL land posterior density plot AICc ----
svl_land_post_BM1_AICc_df <- data.frame(AICc = svl_land_post_BM1_AICc, model = 'BM1',
                                        trait = 'land', tree = 'post')
svl_land_post_BMS_AICc_df <- data.frame(AICc = svl_land_post_BMS_AICc, model = 'BMS',
                                        trait = 'land', tree = 'post')
svl_land_post_OU1_AICc_df <- data.frame(AICc = svl_land_post_OU1_AICc, model = 'OU1',
                                        trait = 'land', tree = 'post')

aic_df_svl_land_post <- rbind(svl_land_post_BM1_AICc_df, svl_land_post_BMS_AICc_df, svl_land_post_OU1_AICc_df)
saveRDS(aic_df_svl_land_post, 'objects/OUwie/aic_df_svl_land_post.rds')
aic_df_svl_land_post <- readRDS('objects/OUwie/aic_df_svl_land_post.rds')

aic_plot_svl_land_post <- ggplot(aic_df_svl_land_post) +
  geom_density(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.7) +
  scale_fill_viridis_d()


# 1.2. SVL habitat_broad ----
# Three habitat categories: ground, rock, tree.
rownames(morpho)
colnames(morpho)
svl_habitat_broad_df <- morpho %>%
  select(species, habitat_broad, SVL)

# 1.2.1. SVL habitat_broad consensus ----
length(pd_habitat_broad$tree)

# SVL habitat_broad consensus BM1 ----
svl_habitat_broad_cons_BM1 <- list()
for (i in 1:length(pd_habitat_broad$tree)){
  print(i)
  svl_habitat_broad_cons_BM1[[i]] <- OUwie(pd_habitat_broad$tree[[i]], svl_habitat_broad_df, model="BM1", 
                                  simmap.tree=TRUE, root.station=FALSE, 
                                  algorithm = "three.point", scaleHeight = TRUE)
}
saveRDS(svl_habitat_broad_cons_BM1, 'objects/OUwie/svl_habitat_broad_cons_BM1.rds', version = 2)
svl_habitat_broad_cons_BM1 <- readRDS('objects/OUwie/svl_habitat_broad_cons_BM1.rds')

# Get AICc values
svl_habitat_broad_cons_BM1_AICc <- c()
for (i in 1:length(pd_habitat_broad$tree)){
  svl_habitat_broad_cons_BM1_AICc[i] <- svl_habitat_broad_cons_BM1[[i]]$AICc
}

hist(svl_habitat_broad_cons_BM1_AICc)
svl_habitat_broad_cons_BM1_AICc_mean <- mean(svl_habitat_broad_cons_BM1_AICc)

# SVL habitat_broad consensus BMS ----
svl_habitat_broad_cons_BMS <- list()
for (i in 1:length(pd_habitat_broad$tree)){
  print(i)
  svl_habitat_broad_cons_BMS[[i]] <- OUwie(pd_habitat_broad$tree[[i]], svl_habitat_broad_df, model="BMS", 
                                  simmap.tree=TRUE, root.station=FALSE, 
                                  algorithm = "three.point", scaleHeight = TRUE)
}
saveRDS(svl_habitat_broad_cons_BMS, 'objects/OUwie/svl_habitat_broad_cons_BMS.rds', version = 2)
svl_habitat_broad_cons_BMS <- readRDS('objects/OUwie/svl_habitat_broad_cons_BMS.rds')

# Get AICc values
svl_habitat_broad_cons_BMS_AICc <- c()
for (i in 1:length(pd_habitat_broad$tree)){
  svl_habitat_broad_cons_BMS_AICc[i] <- svl_habitat_broad_cons_BMS[[i]]$AICc
}

hist(svl_habitat_broad_cons_BMS_AICc)
svl_habitat_broad_cons_BMS_AICc_mean <- mean(svl_habitat_broad_cons_BMS_AICc)

# ...Plot rates BMS (sigmasq) ----
sigmasq_ground_cons_BMS <- c()
sigmasq_rock_cons_BMS <- c()
sigmasq_tree_cons_BMS <- c()
for (i in 1:length(svl_habitat_broad_cons_BMS)){
  #print(i)
  solution <- as.data.frame(svl_habitat_broad_cons_BMS[[i]]$solution)
  g <- solution['sigma.sq', 'ground']
  r <- solution['sigma.sq', 'rock']
  t <- solution['sigma.sq', 'tree']
  sigmasq_ground_cons_BMS <- c(sigmasq_ground_cons_BMS, g)
  sigmasq_rock_cons_BMS <- c(sigmasq_rock_cons_BMS, r)
  sigmasq_tree_cons_BMS <- c(sigmasq_tree_cons_BMS, t)
}

sigmasq_ground_cons_BMS_df <- data.frame(sigmasq = sigmasq_ground_cons_BMS, 
                                         trait = 'habitat', state = 'ground')
sigmasq_rock_cons_BMS_df <- data.frame(sigmasq = sigmasq_rock_cons_BMS, 
                                       trait = 'habitat', state = 'rock')
sigmasq_tree_cons_BMS_df <- data.frame(sigmasq = sigmasq_tree_cons_BMS, 
                                       trait = 'habitat', state = 'tree')
sigmasq_habitat_broad_cons_BMS_df <- rbind(sigmasq_ground_cons_BMS_df, 
                                           sigmasq_rock_cons_BMS_df,
                                           sigmasq_tree_cons_BMS_df)
saveRDS(sigmasq_habitat_broad_cons_BMS_df, 'objects/OUwie/svl_habitat_cons_BMS_rate_df.rds')
sigmasq_habitat_broad_cons_BMS_df <- readRDS('objects/OUwie/svl_habitat_cons_BMS_rate_df.rds')


rate_plot_habitat_broad_cons <- ggplot(sigmasq_habitat_broad_cons_BMS_df) +
  geom_density(aes(x = sigmasq, fill = habitat_broad), color = 'transparent') +
  scale_fill_manual(values = habitat_broad_colors) +
  labs(title = 'SVL rates consensus tree')

# SVL habitat_broad consensus OU1 ----
svl_habitat_broad_cons_OU1 <- list()
for (i in 1:length(pd_habitat_broad$tree)){
  print(i)
  svl_habitat_broad_cons_OU1[[i]] <- OUwie(pd_habitat_broad$tree[[i]], svl_habitat_broad_df, model="OU1", 
                                  simmap.tree=TRUE, root.station=TRUE, 
                                  algorithm = "three.point", scaleHeight = TRUE)
}
saveRDS(svl_habitat_broad_cons_OU1, 'objects/OUwie/svl_habitat_broad_cons_OU1.rds', version = 2)
svl_habitat_broad_cons_OU1 <- readRDS('objects/OUwie/svl_habitat_broad_cons_OU1.rds')

# Get AICc values
svl_habitat_broad_cons_OU1_AICc <- c()
for (i in 1:length(pd_habitat_broad$tree)){
  svl_habitat_broad_cons_OU1_AICc[i] <- svl_habitat_broad_cons_OU1[[i]]$AICc
}

hist(svl_habitat_broad_cons_OU1_AICc)
svl_habitat_broad_cons_OU1_AICc_mean <- mean(svl_habitat_broad_cons_OU1_AICc)


# SVL habitat_broad consensus model selection ----
svl_habitat_broad_cons_BM1_AICc_mean
svl_habitat_broad_cons_BMS_AICc_mean
svl_habitat_broad_cons_OU1_AICc_mean
meanAICc_svl_habitat_broad_cons <- data.frame(model = c('BM1', 'BMS', 'OU1'), 
                                     AICc = c(svl_habitat_broad_cons_BM1_AICc_mean, 
                                              svl_habitat_broad_cons_BMS_AICc_mean,
                                              svl_habitat_broad_cons_OU1_AICc_mean))
saveRDS(meanAICc_svl_habitat_broad_cons, 'objects/OUwie/meanAICc_svl_habitat_cons.rds')
meanAICc_svl_habitat_cons <- readRDS('objects/OUwie/meanAICc_svl_habitat_cons.rds')

meanAICc_svl_habitat_broad_cons$model[meanAICc_svl_habitat_broad_cons$AICc == min(meanAICc_svl_habitat_broad_cons$AICc)]
# The minimum AICc value corresponds to the best fitting model: BMS

# SVL habitat_broad consensus density plot AICc ----
svl_habitat_broad_cons_BM1_AICc_df <- data.frame(AICc = svl_habitat_broad_cons_BM1_AICc, model = 'BM1', 
                                                 trait = 'habitat', tree = 'cons')
svl_habitat_broad_cons_BMS_AICc_df <- data.frame(AICc = svl_habitat_broad_cons_BMS_AICc, model = 'BMS',
                                                 trait = 'habitat', tree = 'cons')
svl_habitat_broad_cons_OU1_AICc_df <- data.frame(AICc = svl_habitat_broad_cons_OU1_AICc, model = 'OU1',
                                                 trait = 'habitat', tree = 'cons')

aic_df_svl_habitat_cons <- rbind(svl_habitat_broad_cons_BM1_AICc_df, 
                svl_habitat_broad_cons_BMS_AICc_df, 
                svl_habitat_broad_cons_OU1_AICc_df
                )
saveRDS(aic_df_svl_habitat_cons, 'objects/OUwie/aic_df_svl_habitat_cons.rds')
aic_df_svl_habitat_cons <- readRDS('objects/OUwie/aic_df_svl_habitat_cons.rds')


aic_plot_svl_habitat_broad_cons <- ggplot(aic_df) +
  geom_histogram(mapping = aes(x = AICc, fill = model), color = 'transparent', alpha = 0.7) +
  scale_fill_viridis_d()

aic_df %>%
  filter(model %in% c('BMS', 'OU1')) %>%
  ggplot() +
  geom_histogram(mapping = aes(x = AICc, fill = model), color = 'transparent') +
  scale_fill_viridis_d()


# 1.2.2. SVL habitat_broad posterior ----
# 100 trees, with 100 simulations each. 
ntrees <- 100
nsim <- 100


# SVL habitat_broad posterior BM1 ----
length(pd_habitat_broad_post)
length(pd_habitat_broad_post[[1]]$tree)
pd_habitat_broad_post[[1]]$tree[[1]]


svl_habitat_broad_post_BM1 <- list()
for (i in 1:ntrees){
  svl_habitat_broad_post_BM1[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_habitat_broad_post_BM1[[i]][[j]] <- OUwie(pd_habitat_broad_post[[i]]$tree[[j]], svl_habitat_broad_df, 
                                         model="BM1", 
                                         simmap.tree=TRUE, root.station=FALSE, 
                                         algorithm = "three.point", scaleHeight = TRUE)
    
  }
}
saveRDS(svl_habitat_broad_post_BM1, 'objects/OUwie/svl_habitat_broad_post_BM1.rds', version = 2)
svl_habitat_broad_post_BM1 <- readRDS('objects/OUwie/svl_habitat_broad_post_BM1.rds')

# Get AICc values
svl_habitat_broad_post_BM1_AICc <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_habitat_broad_post_BM1_AICc <- c(svl_habitat_broad_post_BM1_AICc, svl_habitat_broad_post_BM1[[i]][[j]]$AICc)
  }
}

hist(svl_habitat_broad_post_BM1_AICc)
svl_habitat_broad_post_BM1_AICc_mean <- mean(svl_habitat_broad_post_BM1_AICc)


# SVL habitat_broad posterior BMS ----
svl_habitat_broad_post_BMS <- list()
for (i in 1:ntrees){
  svl_habitat_broad_post_BMS[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_habitat_broad_post_BMS[[i]][[j]] <- OUwie(pd_habitat_broad_post[[i]]$tree[[j]], svl_habitat_broad_df, 
                                         model="BMS", 
                                         simmap.tree=TRUE, root.station=FALSE, 
                                         algorithm = "three.point", scaleHeight = TRUE)
    
  }
}
saveRDS(svl_habitat_broad_post_BMS, 'objects/OUwie/svl_habitat_broad_post_BMS.rds', version = 2)
svl_habitat_broad_post_BMS <- readRDS('objects/OUwie/svl_habitat_broad_post_BMS.rds')

# Get AICc values
svl_habitat_broad_post_BMS_AICc <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_habitat_broad_post_BMS_AICc <- c(svl_habitat_broad_post_BMS_AICc, svl_habitat_broad_post_BMS[[i]][[j]]$AICc)
  }
}

hist(svl_habitat_broad_post_BMS_AICc)
svl_habitat_broad_post_BMS_AICc_mean <- mean(svl_habitat_broad_post_BMS_AICc)


# ...Plot rates BMS (sigmasq) ----
sigmasq_ground_post_BMS <- c()
sigmasq_rock_post_BMS <- c()
sigmasq_tree_post_BMS <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    solution <- as.data.frame(svl_habitat_broad_post_BMS[[i]][[j]]$solution)
    g <- solution['sigma.sq', 'ground']
    r <- solution['sigma.sq', 'rock']
    t <- solution['sigma.sq', 'tree']
    sigmasq_ground_post_BMS <- c(sigmasq_ground_post_BMS, g)
    sigmasq_rock_post_BMS <- c(sigmasq_rock_post_BMS, r)
    sigmasq_tree_post_BMS <- c(sigmasq_tree_post_BMS, t)
  }
}

sigmasq_ground_post_BMS_df <- data.frame(sigmasq = sigmasq_ground_post_BMS, 
                                         trait = 'habitat', state = 'ground')
sigmasq_rock_post_BMS_df <- data.frame(sigmasq = sigmasq_rock_post_BMS, 
                                       trait = 'habitat', state = 'rock')
sigmasq_tree_post_BMS_df <- data.frame(sigmasq = sigmasq_tree_post_BMS, 
                                       trait = 'habitat', state = 'tree')

sigmasq_habitat_post_BMS_df <- rbind(sigmasq_ground_post_BMS_df,
                                     sigmasq_rock_post_BMS_df,
                                     sigmasq_tree_post_BMS_df)
saveRDS(sigmasq_habitat_post_BMS_df, 'objects/OUwie/svl_habitat_post_BMS_rate_df.rds')
sigmasq_habitat_post_BMS_df <- readRDS('objects/OUwie/svl_habitat_post_BMS_rate_df.rds')


rate_plot_habitat_broad_post <- ggplot(sigmasq_habitat_post_BMS_df) +
  geom_density(aes(x = sigmasq, fill = habitat_broad), color = 'transparent') +
  scale_fill_manual(values = habitat_broad_colors) +
  labs(title = 'SVL rates posterior trees')


# SVL habitat_broad posterior OU1 ----
svl_habitat_broad_post_OU1 <- list()
for (i in 1:ntrees){
  svl_habitat_broad_post_OU1[[i]] <- list()
}

for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_habitat_broad_post_OU1[[i]][[j]] <- OUwie(pd_habitat_broad_post[[i]]$tree[[j]], svl_habitat_broad_df, 
                                         model="OU1", 
                                         simmap.tree=TRUE, root.station=TRUE, 
                                         algorithm = "three.point", scaleHeight = TRUE)
    
  }
}
saveRDS(svl_habitat_broad_post_OU1, 'objects/OUwie/svl_habitat_broad_post_OU1.rds', version = 2)
svl_habitat_broad_post_OU1 <- readRDS('objects/OUwie/svl_habitat_broad_post_OU1.rds')

# Get AICc values
svl_habitat_broad_post_OU1_AICc <- c()
for (i in 1:ntrees){
  for (j in 1:nsim){
    print(paste0('tree ', i, '; sim ', j))
    svl_habitat_broad_post_OU1_AICc <- c(svl_habitat_broad_post_OU1_AICc, svl_habitat_broad_post_OU1[[i]][[j]]$AICc)
  }
}

hist(svl_habitat_broad_post_OU1_AICc)
svl_habitat_broad_post_OU1_AICc_mean <- mean(svl_habitat_broad_post_OU1_AICc)

# SVL habitat_broad posterior model selection ----
svl_habitat_broad_post_BM1_AICc_mean
svl_habitat_broad_post_BMS_AICc_mean
svl_habitat_broad_post_OU1_AICc_mean
meanAICc_svl_habitat_broad_post <- data.frame(model = c('BM1', 'BMS', 'OU1'), 
                                     AICc = c(svl_habitat_broad_post_BM1_AICc_mean, 
                                              svl_habitat_broad_post_BMS_AICc_mean,
                                              svl_habitat_broad_post_OU1_AICc_mean))

saveRDS(meanAICc_svl_habitat_broad_post, 'objects/OUwie/meanAICc_svl_habitat_post.rds')
meanAICc_svl_habitat_post <- readRDS('objects/OUwie/meanAICc_svl_habitat_post.rds')

meanAICc_svl_habitat_broad_post$model[meanAICc_svl_habitat_broad_post$AICc == min(meanAICc_svl_habitat_broad_post$AICc)]
# The minimum AICc value corresponds to the best fitting model: BMS.

# SVL habitat_broad posterior density plot AICc ----
svl_habitat_broad_post_BM1_AICc_df <- data.frame(AICc = svl_habitat_broad_post_BM1_AICc, model = 'BM1',
                                                 trait = 'habitat', tree = 'post')
svl_habitat_broad_post_BMS_AICc_df <- data.frame(AICc = svl_habitat_broad_post_BMS_AICc, model = 'BMS',
                                                 trait = 'habitat', tree = 'post')
svl_habitat_broad_post_OU1_AICc_df <- data.frame(AICc = svl_habitat_broad_post_OU1_AICc, model = 'OU1',
                                                 trait = 'habitat', tree = 'post')

aic_df_svl_habitat_post <- rbind(svl_habitat_broad_post_BM1_AICc_df, 
                                 svl_habitat_broad_post_BMS_AICc_df, 
                                 svl_habitat_broad_post_OU1_AICc_df)
saveRDS(aic_df_svl_habitat_post, 'objects/OUwie/aic_df_svl_habitat_post.rds')
aic_df_svl_habitat_post <- readRDS('objects/OUwie/aic_df_svl_habitat_post.rds')

aic_plot_svl_habitat_broad_post <- ggplot(aic_df_svl_habitat_post) +
  geom_density(mapping = aes(x = AICc, fill = model), 
               color = 'transparent', alpha = 0.7) +
  scale_fill_viridis_d()

##############################################
##############################################
##############################################

# PLOTS - PLOTS - PLOTS - PLOTS - PLOTS ----
# Playing with plots ----
# https://patchwork.data-imaginist.com/
aic_plot_svl_land_cons
aic_plot_svl_land_post

rate_plot_land_cons +
  annotate(geom = 'text', x = 0, y = 500, 
           label = 'BMS model on 1,000\nstochastic character maps',
           hjust = -0.1, 
           vjust = 0.5)


rate_plot_land_cons

rate_plot_land_post +
  labs(subtitle = 'From the fit of a BMS model\n
       on 100 posterior trees, each with 100 stochastic\n
       character mapping') +
  theme(plot.subtitle = element_text(size = 12, vjust = 1, hjust = 0.5))


(aic_plot_svl_land_cons + aic_plot_svl_land_post) /
  (rate_plot_land_cons + rate_plot_land_post)

aic_plot_svl_habitat_broad_cons
aic_plot_svl_habitat_broad_post$layers[[1]]$aes_params$alpha = 0.9

rate_plot_habitat_broad_cons
rate_plot_habitat_broad_post

install.packages('patchwork')


rate_plot_habitat_broad_cons + rate_plot_habitat_broad_post



