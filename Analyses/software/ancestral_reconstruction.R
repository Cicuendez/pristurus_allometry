
# Ancestral reconstruction


## Load packages
libs <- c('tidyverse', 'treeio', 'phytools', 'geiger', 'ggtree', 'rBEAST')
lapply(libs, require, character.only = TRUE)

## Import data

#### Import phylogeny
tree <- read.nexus('data/phylogeny/pristurus_tree_final.nex')
posterior_pristurus <- readRDS('data/phylogeny/posterior_pristurus_final.rds')

#### Import morpho data
morpho <- read.table('data/morphology/morpho_sp_final.csv', sep=';', dec='.', 
                     header = TRUE)
rownames(morpho) <- morpho$species

# log10-transform morpho data
morpho_log <- morpho %>%
  mutate(across(where(is.numeric), log10))
rownames(morpho_log) <- morpho_log$species

## Ancestral reconstruction
# We will study the evolution of different traits on the phylogeny: 
# 1) Insularity, with two different states: mainland and island.
# 2) habitat (broad classification), three states: rock, tree, ground.
# 3) habitat, four states: rock, tree, soft-ground, hard-ground. 
# 
# First we fit different models of discrete character evolution: equal-rates (ER), 
# symmetrical (SYM), and all-rates-different (ARD). Then we will reconstruct the 
# discrete character across the phylogeny with the best-fit model. 
# We will do it for the consensus tree and then for 100 trees from the posterior.
# We will do 1,000 simulations for the consensus ancestral reconstruction, and 
# 100 simulations for the posterior trees. 

ntrees <- 100
nsim <- 100

#### 1) Insularity ancestral reconstruction
# Set colors and create trait variable
#land_colors0 <- c(mainland = 'coral2', island = 'seagreen')
#land_colors1 <- c(mainland = '#A54E29', island = '#42978A')
#land_colors2 <- c(mainland = '#C56542', island = '#67D1B7')
land_colors <- c(mainland = '#C56542', island = '#42978A')
plot(1:2, 1:2, col = land_colors, cex = 5, pch = 16)
land <- morpho$land
names(land) <- row.names(morpho)

# Fit models of discrete character evolution
fit_ER_land <- fitDiscrete(tree, land, model = 'ER')
fit_SYM_land <- fitDiscrete(tree, land, model = 'SYM')
fit_ARD_land <- fitDiscrete(tree, land, model = 'ARD')
aicw(c(fit_ER_land$opt$aicc, fit_SYM_land$opt$aicc, 
       fit_ARD_land$opt$aicc))
# Best fit model for land: ER

# Make simmap with consensus tree
simmap_land <- make.simmap(tree = tree, x = land, model = 'ER', nsim = 1000)
plotSimmap(simmap_land, colors = land_colors)
pd_land <- describe.simmap(simmap_land, plot=FALSE)

# Save simmap object
saveRDS(pd_land, 'objects/anc_rec/pd_land_cons.rds')

pd_land <- readRDS("objects/anc_rec/pd_land_cons.rds")
pd_land$ace
# Root: island 0.351, mainland 0.649.
plot(tree)
nodelabels()

# Plot reconstruction
pdf('plots/ancestral_reconstructions/anc_insularity.pdf', paper = 'a4r')
plot(pd_land, fsize=1, ftype="i", colors = land_colors)
add.simmap.legend(colors=land_colors, prompt=FALSE, fsize=1, x = 0, y = 10)
dev.off()

# Make simmap with posterior trees
#sampling_trees <- sample(1:length(posterior_pristurus), size = ntrees, replace = FALSE)
simmap_land_post <- make.simmap(tree = posterior_pristurus[1:ntrees], 
                                x = land, model = 'ER', nsim = nsim)

plotSimmap(simmap_land_post, colors = land_colors)
pd_land_post <- list()
n <- 1
m <- nsim
for (i in 1:ntrees){
  print(paste(n, m, sep = ' to '))
  pd_land_post[[i]] <-  describe.simmap(simmap_land_post[n:m], plot=FALSE)
  n <- n+nsim
  m <- m+nsim
}
plot(pd_land_post[[i]])

# pd_land_post is a list with ntrees simmap summaries of nsim simulations each.
saveRDS(pd_land_post, 'objects/anc_rec/pd_land_post.rds')

#### 2) Habitat broad ancestral reconstruction
# Set colors and create trait variable
#habitat_broad_colors0 <- c(ground = 'brown', rock = 'gray', tree = 'darkgreen')
#habitat_broad_colors <- c(ground = '#e0710b', rock = '#6C586E', tree = '#119616')
habitat_broad_colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")
plot(1:3, 1:3, col = habitat_broad_colors, cex = 10, pch = 16)
habitat_broad <- morpho$habitat_broad
names(habitat_broad) <- row.names(morpho)

# Fit models of discrete character evolution
fit_ER_habitat_broad <- fitDiscrete(tree, habitat_broad, model = 'ER')
fit_SYM_habitat_broad <- fitDiscrete(tree, habitat_broad, model = 'SYM')
fit_ARD_habitat_broad <- fitDiscrete(tree, habitat_broad, model = 'ARD')
aicw(c(fit_ER_habitat_broad$opt$aicc, fit_SYM_habitat_broad$opt$aicc, 
       fit_ARD_habitat_broad$opt$aicc))
# Best fit model for habitat broad: ER

# Make simmap with consensus tree
simmap_habitat_broad <- make.simmap(tree = tree, x = habitat_broad, model = 'ER', nsim = 1000)
plotSimmap(simmap_habitat_broad, colors = habitat_broad_colors)
pd_habitat_broad <- describe.simmap(simmap_habitat_broad, plot=FALSE)

# Save simmap object
saveRDS(pd_habitat_broad, 'objects/anc_rec/pd_habitat_broad_cons.rds')

pd_habitat_broad <- readRDS("objects/anc_rec/pd_habitat_broad_cons.rds")
pd_habitat_broad$ace
# Root: ground 0.012, rock 0.978, tree 0.010

# Plot reconstruction
pdf('plots/ancestral_reconstructions/anc_habitat.pdf', paper = 'a4r')
plot(pd_habitat_broad, fsize=1, 
     ftype="i", 
#     ftype = 'off',
     colors = habitat_broad_colors, 
     direction = 'leftwards')
add.simmap.legend(colors=habitat_broad_colors, prompt=FALSE, fsize=1, 
#                  x = 0, y = 10,
                  x = 35, y = 10)
dev.off()

# Plot insularity and habitat reconstructions face to face ----
# Arrange plot layout
layout.matrix <- matrix(c(2, 1, 0, 3), nrow = 2, ncol = 2)
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout.matrix
graphics::layout(mat = layout.matrix,
       heights = c(1, 2), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns
graphics::layout(mat = layout.matrix,
                 heights = 1, # Heights of the two rows
                 widths = c(2, 1)) # Widths of the two columns
layout.show(2)
?layout


pdf('plots/ancestral_reconstructions/ancestral_FaceToFace.pdf', paper = 'a4r')
par(mfrow = c(1, 2))
par(mar = c(5, 4, 0, 0))
plot(pd_land, fsize=0.6, ftype="i", colors = land_colors)
add.simmap.legend(colors=land_colors, prompt=FALSE, fsize=0.7, x = 0, y = 10)

par(mar = c(0,0,0,0))
plot(pd_habitat_broad, fsize=0.6, 
#     ftype="i", 
     ftype = 'off',
     colors = habitat_broad_colors, 
     direction = 'leftwards')
add.simmap.legend(colors=habitat_broad_colors, prompt=FALSE, fsize=0.7, x = 35, y = 10)

dev.off()

par(mfrow = c(1, 1))



# Make simmap with posterior trees
#sampling_trees <- sample(1:length(posterior_pristurus), size = ntrees, replace = FALSE)
simmap_habitat_broad_post <- make.simmap(tree = posterior_pristurus[1:ntrees], 
                                         x = habitat_broad, model = 'ER', nsim = nsim)

plotSimmap(simmap_habitat_broad_post, colors = habitat_broad_colors)
pd_habitat_broad_post <- list()
n <- 1
m <- nsim
for (i in 1:ntrees){
  print(paste(n, m, sep = ' to '))
  pd_habitat_broad_post[[i]] <-  describe.simmap(simmap_habitat_broad_post[n:m], plot=FALSE)
  n <- n+nsim
  m <- m+nsim
}
plot(pd_habitat_broad_post[[i]])
# pd_habitat_broad_post is a list with ntrees simmap summaries of nsim simulations each.

saveRDS(pd_habitat_broad_post, 'objects/anc_rec/pd_habitat_broad_post.rds')

#### 3) Habitat ancestral reconstruction
# Set colors and create trait variable
habitat_colors0 <- c('hard-ground' = 'brown', 'soft-ground' = 'orange', 
                     rock = 'gray', tree = 'darkgreen')
habitat_colors <- c('hard-ground' = '#D63916', 'soft-ground' = '#E9A800', 
                    rock = '#6C586E', tree = '#119616')
habitat <- morpho$habitat
names(habitat) <- row.names(morpho)

# Fit models of discrete character evolution
fit_ER_habitat <- fitDiscrete(tree, habitat, model = 'ER')
fit_SYM_habitat <- fitDiscrete(tree, habitat, model = 'SYM')
fit_ARD_habitat <- fitDiscrete(tree, habitat, model = 'ARD')
aicw(c(fit_ER_habitat$opt$aicc, fit_SYM_habitat$opt$aicc, 
       fit_ARD_habitat$opt$aicc))
# Best fit model for habitat: ER

# Make simmap with consensus tree
simmap_habitat <- make.simmap(tree = tree, x = habitat, model = 'ER', nsim = 1000)
plotSimmap(simmap_habitat, colors = habitat_colors)
pd_habitat <- describe.simmap(simmap_habitat, plot=FALSE)

# Save simmap object
saveRDS(pd_habitat, 'objects/anc_rec/pd_habitat_cons.rds')

pd_habitat <- readRDS("objects/anc_rec/pd_habitat_cons.rds")

# Plot reconstruction
plot(pd_habitat, fsize=0.6, ftype="i", colors = habitat_colors)
add.simmap.legend(colors=habitat_colors, prompt=FALSE, fsize=0.7, x = 0, y = 10)

# Make simmap with posterior trees
#sampling_trees <- sample(1:length(posterior_pristurus), size = ntrees, replace = FALSE)
nsim <- 100
ntrees <- 100
#sampling_trees <- sample(1:length(posterior_pristurus), size = ntrees, replace = FALSE)
simmap_habitat_post <- make.simmap(tree = posterior_pristurus[1:ntrees], 
                                   x = habitat, model = 'ER', nsim = nsim)

plotSimmap(simmap_habitat_post, colors = habitat_colors)
pd_habitat_post <- list()
n <- 1
m <- nsim
for (i in 1:ntrees){
  print(paste(n, m, sep = ' to '))
  pd_habitat_post[[i]] <-  describe.simmap(simmap_habitat_post[n:m], plot=FALSE)
  n <- n+nsim
  m <- m+nsim
}
plot(pd_habitat_post[[i]])
# pd_habitat_post is a list with ntrees simmap summaries of nsim simulations each.

saveRDS(pd_habitat_post, 'objects/anc_rec/pd_habitat_post.rds')

## Conclusion
# Now we have the *Pristurus* phylogeny (and the posterior trees) mapped with the discrete traits. 
# We will use them for subsequent analyses of evolutionary rates of continuous variables 
# (size and shape), exploring how the insularity or the habitat specialization have 
# affected morphological evolution in this genus. 











