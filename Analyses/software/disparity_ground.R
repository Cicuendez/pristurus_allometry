# GROUND DISPARITY 
# Explanation: ----
# Is morphological disparity higher in the ground species than in the rest? 
# We compute the actual disparity in ground and no-ground, and we simulate datasets 
# by simulating the evolution of the trait, to see if the disparity difference observed 
# between ground and no-ground is different from expected. 
# We calculate the ratio disp_ground/disp_noground, for the actual data (empiric part) 
# and for the simulations. If the observed ratio is significantly higher than 
# expected based on the simulations, then we'll assume that disparity in ground is 
# effectively higher, and that there is a 'habitat effect' for the trait under consideration 
# (either SVL, PC1 or PC1)

# Set seed
set.seed(19042021)

# Load packages ----
libs <- c('phytools', 'geiger', 'tidyverse')
lapply(libs, require, character.only = TRUE)

# Import data ----
tree <- read.nexus('data/phylogeny/pristurus_tree_final.nex')
morpho <- read.table('objects/phypca/phypca_scores.csv', sep = ";", dec = '.', 
                     header = TRUE, row.names = 1)
name.check(tree, morpho)

# Create a column 'hab_bin' (habitat_binary) in the morpho dataset, 
# with two categories: ground and noground.
morpho <- morpho %>% 
  mutate(hab_bin = case_when(habitat_broad == 'ground' ~ 'ground', 
                             T ~ 'noground'))

# 1. SVL ----
## 1.1. Empiric part ---- 
### Split ground and no-ground ----
svl_ground <- morpho$SVL[morpho$hab_bin=="ground"]
names(svl_ground) <- rownames(morpho)[morpho$hab_bin=="ground"]

svl_noground <- morpho$SVL[morpho$hab_bin=="noground"]
names(svl_noground) <- rownames(morpho)[morpho$hab_bin=="noground"]

### Calculate disparity of island and mainland separately ----
dispSVL_ground <- disparity(data=svl_ground, index="avg.sq")
dispSVL_noground <- disparity(data=svl_noground, index="avg.sq")

### Calculate disparity ratio island/mainland ----
dispSVL_ratio_emp <- dispSVL_ground/dispSVL_noground
dispSVL_ratio_emp
## If the value is greater than 1, the disparity is higher in the ground ## 
## If the value is less than 1, the disparity is higher in no-ground ##
### Ratio SVL disparity ground/no-ground = 2.25 -> Disparity is higher in ground
## END OF THE EMPIRIC PART ##

## 1.2. Simulation part ----
### Question: Is the difference between island and mainland disparity significant?
### Model of trait evolution ----
### Fit models of trait evolution on the tree
svl <- morpho$SVL
names(svl) <- rownames(morpho)
svl_bm <- fitContinuous(phy = tree, dat = svl, model = 'BM')
svl_ou <- fitContinuous(phy = tree, dat = svl, model = 'OU')
svl_eb <- fitContinuous(phy = tree, dat = svl, model = 'EB')
svl_rt <- fitContinuous(phy = tree, dat = svl, model = 'rate_trend')
svl_lambda <- fitContinuous(phy = tree, dat = svl, model = 'lambda')
svl_kappa <- fitContinuous(phy = tree, dat = svl, model = 'kappa')
svl_delta <- fitContinuous(phy = tree, dat = svl, model = 'delta')
svl_mt <- fitContinuous(phy = tree, dat = svl, model = 'mean_trend')
svl_white <- fitContinuous(phy = tree, dat = svl, model = 'white')

### Model selection ----
### Minimum AIC
?AIC
svl_aic <- AIC(svl_bm, svl_ou, svl_eb)
#AIC(svl_bm, svl_ou, svl_eb, svl_rt, svl_lambda, svl_kappa, svl_delta, svl_mt, svl_white)
rownames(svl_aic)[svl_aic$AIC == min(svl_aic$AIC)]
# The best fit model for SVL is BM

### Get a rate value (sigma square) ----
### Save sigma square from best model and put it later in 'par' option in the function 'sim.char'
sigmasq_svl <- svl_bm$opt$sigsq

### Set the number of simulations ----
nsim=10000

### Simulate the evolution of a trait with the sigma obtained in the best fit model ----
sim_svl <- sim.char(tree, par = sigmasq_svl, model = "BM", nsim = nsim)
dim(sim_svl)
### The result is a 3 dimensional dataset [x, y, z]:
### x) Number of taxon; 
### y) Number of variable (in this case only one); 
### z) Number of simulation. 

### Assign ground or no ground to the species in each simulation ----
sim_habitat_svl <- list()
for (i in 1:nsim){
  sim_habitat_svl[[i]] <- cbind(as.factor(morpho$hab_bin), sim_svl[,1,i])
}
## The result is a list of nsim datasets (one for each simulation) with the values of the trait 
## and the assignment to ground (1) and no-ground (2) ##

### Split ground (1) and no-ground (2) in each simulation ##
sim_ground_svl <- list()
for (i in 1:length(sim_habitat_svl)){
  sim_ground_svl[[i]] <- sim_habitat_svl[[i]][,2][sim_habitat_svl[[i]]==1]
}

sim_noground_svl <- list()
for (i in 1:length(sim_habitat_svl)){
  sim_noground_svl[[i]] <- sim_habitat_svl[[i]][,2][sim_habitat_svl[[i]]==2]
}
## The result is two lists of nsim simulated datasets, 
## one with ground species and another one with no-ground species, 
## and their simulated values of trait (svl)

## CALCULATE THE DISPARITY IN GROUND AND NO-GROUND SEPARATELY ##
disp_sim_ground_svl <- list()
for (i in 1:length(sim_ground_svl)){
  disp_sim_ground_svl[[i]] <- disparity(data=sim_ground_svl[[i]], 
                                        index="avg.sq")
}

disp_sim_noground_svl <- list()
for (i in 1:length(sim_noground_svl)){
  disp_sim_noground_svl[[i]] <- disparity(data=sim_noground_svl[[i]], 
                                          index="avg.sq")
}
## The result is two lists of values of disparity, one for island and one for mainland ##

### Calculate disparity ratio (ground/no-ground) for each pair of datasets ----
disp_ratio_sim_svl <- c()
for (i in 1:length(disp_sim_noground_svl)){
  disp_ratio_sim_svl[i] <- disp_sim_ground_svl[[i]]/disp_sim_noground_svl[[i]]
}
disp_ratio_sim_svl
## The result is a vector with n values of disparity ratios ##
## End of simulation part

## 1.3. Histogram of ratio values ----
hist(disp_ratio_sim_svl, col="gray", border="white")

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

simratios_svl <- data.frame(sim = c(1:nsim), ratio = disp_ratio_sim_svl)
head(simratios_svl)
range(simratios_svl$ratio)

ggplot(simratios_svl) +
  geom_histogram(aes(x = ratio), fill = 'gray80', color = 'white') +
  theme.clean() +
  xlim(0, 4)



## p-value for the empiric disparity ratio ----
## To do this, sum all the values of ratio equal or greater than the empirical ratio, 
## and divide by the number of simulations
pval_svl <- sum(disp_ratio_sim_svl >= dispSVL_ratio_emp)/nsim
pval_svl # 0.0291
## If the value is less than 0.05, the difference between ground and no-ground is significant (greater than 
## expected at random circumstances)

# Plot disparity ratio SVL ----
(hist_svl_ground <- ggplot(simratios_svl) +
  geom_histogram(aes(x = ratio), fill = 'gray80', color = 'white', 
                 ) +
  xlim(NA, 3.5) +
  geom_segment(mapping = aes(x = dispSVL_ratio_emp, y = 200, xend = dispSVL_ratio_emp, yend = 0), 
               arrow = arrow(angle=30, length = unit(0.2, "cm"), type = 'open'), 
               size = 0.3, color = 'red') +
  annotate(geom = 'text', x = 2, y = 1150, color = 'black', 
           label = 'Observed disparity ratio:') +
  annotate(geom = 'text', x = 2, y = 1000, color = 'black', 
           label = paste0('r = ', round(dispSVL_ratio_emp, digits = 2), 
                          ', p = ', round(pval_svl, digits = 2))) +
  xlab('Disparity ratio (r)') +
  ylab('Frequency') +
  ggtitle(label = 'Disparity in BODY SIZE\nground vs no-ground') +
  theme.clean()
)
#



#### 2. pPC1 ####
## 2.1. EMPIRIC PART ##
## SPLIT GROUND AND NO-GROUND ##

pPC1_ground <- morpho$PC1[morpho$hab_bin=="ground"]
names(pPC1_ground) <- rownames(morpho)[morpho$hab_bin=="ground"]

pPC1_noground <- morpho$PC1[morpho$hab_bin=="noground"]
names(pPC1_noground) <- rownames(morpho)[morpho$hab_bin=="noground"]

## CALCULATE THE DISPARITY OF GROUND AND NO-GROUND SEPARATELY ##
disppPC1_ground <- disparity(data=pPC1_ground, index="avg.sq")
disppPC1_noground <- disparity(data=pPC1_noground, index="avg.sq")

## CALCULATE THE DISPARITY RATIO ISLAND/MAINLAND ##
disppPC1_ratio_emp <- disppPC1_ground/disppPC1_noground
disppPC1_ratio_emp
## If the value is greater than 1, the disparity is higher in the ground ## 
## If the value is less than 1, the disparity is higher in no-ground ##
## END OF THE EMPIRIC PART ##
## The value is not greater than 1, so the disparity for pPC1 is higher in no-ground ##
## There's no need to do the simulation part then ##

#### 3. pPC2 ####
## 3.1. EMPIRIC PART ##
## SPLIT GROUND AND NO-GROUND ##

pPC2_ground <- morpho$PC2[morpho$hab_bin=="ground"]
names(pPC2_ground) <- rownames(morpho)[morpho$hab_bin=="ground"]

pPC2_noground <- morpho$PC2[morpho$hab_bin=="noground"]
names(pPC2_noground) <- rownames(morpho)[morpho$hab_bin=="noground"]

## CALCULATE THE DISPARITY OF GROUND AND NO-GROUND SEPARATELY ##
disppPC2_ground <- disparity(data=pPC2_ground, index="avg.sq")
disppPC2_noground <- disparity(data=pPC2_noground, index="avg.sq")

## CALCULATE THE DISPARITY RATIO GROUND/NO-GROUND ##
disppPC2_ratio_emp <- disppPC2_ground/disppPC2_noground
disppPC2_ratio_emp
## If the value is greater than 1, the disparity is higher in the ground ## 
## If the value is less than 1, the disparity is higher in no-ground ##
## END OF THE EMPIRIC PART ##
## The value is greater than 1 ##
## Ratio PC2 disparity ground/no-ground = 2.47 -> PC2 disparity is higher in ground.

## 3.2. Simulation part ----
### Question: Is the difference between island and mainland disparity significant?
### Model of trait evolution ----
### Fit models of trait evolution on the tree
pPC2 <- morpho$PC2
names(pPC2) <- rownames(morpho)
pPC2_bm <- fitContinuous(phy = tree, dat = pPC2, model = 'BM')
pPC2_ou <- fitContinuous(phy = tree, dat = pPC2, model = 'OU')
pPC2_eb <- fitContinuous(phy = tree, dat = pPC2, model = 'EB')
pPC2_rt <- fitContinuous(phy = tree, dat = pPC2, model = 'rate_trend')
pPC2_lambda <- fitContinuous(phy = tree, dat = pPC2, model = 'lambda')
pPC2_kappa <- fitContinuous(phy = tree, dat = pPC2, model = 'kappa')
pPC2_delta <- fitContinuous(phy = tree, dat = pPC2, model = 'delta')
pPC2_mt <- fitContinuous(phy = tree, dat = pPC2, model = 'mean_trend')
pPC2_white <- fitContinuous(phy = tree, dat = pPC2, model = 'white')

### Model selection ----
### Minimum AIC
?AIC
pPC2_aic <- AIC(pPC2_bm, pPC2_ou, pPC2_eb)
#AIC(pPC2_bm, pPC2_ou, pPC2_eb, pPC2_rt, pPC2_lambda, pPC2_kappa, pPC2_delta, pPC2_mt, pPC2_white)
rownames(pPC2_aic)[pPC2_aic$AIC == min(pPC2_aic$AIC)]
# The best fit model for pPC2 is BM


### Get a rate value (sigma square) ----
### Save sigma square from best model and put it later in 'par' option in the function 'sim.char'
sigmasq_pPC2 <- pPC2_bm$opt$sigsq

### Set the number of simulations ----
nsim=10000

### Simulate the evolution of a trait with the sigma obtained in the best fit model ----
sim_pPC2 <- sim.char(tree, par = sigmasq_pPC2, model = "BM", nsim = nsim)
dim(sim_pPC2)
### The result is a 3 dimensional dataset [x, y, z]:
### x) Number of taxon; 
### y) Number of variable (in this case only one); 
### z) Number of simulation. 

### Assign island or mainland to the species in each simulation ----
sim_habitat_pPC2 <- list()
for (i in 1:nsim){
  sim_habitat_pPC2[[i]] <- cbind(morpho$habitat, sim_pPC2[,1,i])
}
## The result is a list of nsim datasets (one for each simulation) with the values of the trait 
## and the assignment to mainland (2) and island (1) ##

### Split ground (1) and no-ground (2) in each simulation ##
sim_ground_pPC2 <- list()
for (i in 1:length(sim_habitat_pPC2)){
  sim_ground_pPC2[[i]] <- sim_habitat_pPC2[[i]][,2][sim_habitat_pPC2[[i]]==1]
}

sim_noground_pPC2 <- list()
for (i in 1:length(sim_habitat_pPC2)){
  sim_noground_pPC2[[i]] <- sim_habitat_pPC2[[i]][,2][sim_habitat_pPC2[[i]]==2]
}
## The result is two lists of nsim simulated datasets, 
## one with ground species and another one with no-ground species, 
## and their simulated values of trait (pc2)

## CALCULATE THE DISPARITY IN GROUND AND NO-GROUND SEPARATELY ##
disp_sim_ground_pPC2 <- list()
for (i in 1:length(sim_ground_pPC2)){
  disp_sim_ground_pPC2[[i]] <- disparity(data=sim_ground_pPC2[[i]], 
                                         index="avg.sq")
}

disp_sim_noground_pPC2 <- list()
for (i in 1:length(sim_noground_pPC2)){
  disp_sim_noground_pPC2[[i]] <- disparity(data=sim_noground_pPC2[[i]], 
                                           index="avg.sq")
}
## The result is two lists of values of disparity, one for ground and one for no-ground ##

### Calculate disparity ratio (ground/no-ground) for each pair of datasets ----
disp_ratio_sim_pPC2 <- c()
for (i in 1:length(disp_sim_noground_pPC2)){
  disp_ratio_sim_pPC2[i] <- disp_sim_ground_pPC2[[i]]/disp_sim_noground_pPC2[[i]]
}
disp_ratio_sim_pPC2
## The result is a vector with n values of disparity ratios ##
## End of simulation part

## 3.3. Histogram of ratio values ----
par(mar=c(5,5,2,2))
hist(disp_ratio_sim_pPC2, col="gray", border="white")

simratios_pPC2 <- data.frame(sim = c(1:nsim), ratio = disp_ratio_sim_pPC2)
head(simratios_pPC2)
range(simratios_pPC2$ratio)

ggplot(simratios_pPC2) +
  geom_histogram(aes(x = ratio), fill = 'gray80', color = 'white') +
  theme.clean() +
  xlim(0, 3)

## p-value for the empiric disparity ratio ----
## To do this, sum all the values of ratio equal or greater than the empirical ratio, 
## and divide by the number of simulations
p_val_pPC2 <- sum(disp_ratio_sim_pPC2 >= disppPC2_ratio_emp)/nsim
p_val_pPC2

## If the value is less than 0.05, the difference between island and mainland 
## is significant (greater than expected at random circumstances)
## p-value = 0.0107 for PC2 -> significantly higher ratio ground/no-ground than expected 
## for PC2 (head proportions) 

# Plot disparity ratio pPC2 ----
(hist_pPC2_ground <- ggplot(simratios_pPC2) +
  geom_histogram(aes(x = ratio), fill = 'gray80', color = 'white') +
  xlim(NA, 3) +
  geom_segment(mapping = aes(x = disppPC2_ratio_emp, y = 200, xend = disppPC2_ratio_emp, yend = 0), 
               arrow = arrow(angle=30, length = unit(0.2, "cm"), type = 'open'), 
               size = 0.3, color = 'red') +
  annotate(geom = 'text', x = 1.5, y = 1150, color = 'black', 
           label = 'Observed disparity ratio:') +
  annotate(geom = 'text', x = 1.5, y = 1000, color = 'black', 
           label = paste0('r = ', round(disppPC2_ratio_emp, digits = 2), 
                          ', p = ', round(p_val_pPC2, digits = 2))) +
  xlab('Disparity ratio (r)') +
  ylab('Frequency') +
  ggtitle(label = 'Disparity in HEAD PROPORTIONS\nground vs no-ground') +
  theme.clean()
)



# 4. Plot together ----
(hist_svl_ground <- ggplot(simratios_svl) +
   geom_histogram(aes(x = ratio), fill = 'gray80', color = 'white', 
   ) +
   xlim(NA, 3.5) +
   geom_segment(mapping = aes(x = dispSVL_ratio_emp, y = 170, xend = dispSVL_ratio_emp, yend = 0), 
                arrow = arrow(angle=30, length = unit(0.2, "cm"), type = 'open'), 
                size = 0.3, color = 'red') +
   annotate(geom = 'text', x = 2, y = 800, color = 'black', 
            label = 'Observed disparity ratio:') +
   annotate(geom = 'text', x = 2, y = 700, color = 'black', 
            label = paste0('r = ', round(dispSVL_ratio_emp, digits = 2), 
                           ', p = ', round(pval_svl, digits = 2))) +
   xlab('Disparity ratio (r)') +
   ylab('Frequency') +
   ggtitle(label = 'Disparity in BODY SIZE\nground vs no-ground') +
   theme.clean()
)


histogram_list <- list(hist_svl_ground, hist_pPC2_ground)
histogram_grid <- cowplot::plot_grid(plotlist = histogram_list, ncol = 2, labels = 'auto')

ggsave('plots/disparity/histograms_disparity.pdf', histogram_grid,
       height = 5, width = 10)

dispSVL_ratio_emp
disppPC1_ratio_emp
disppPC2_ratio_emp

# 5. DISPARITY THROUGH TIME ----

svl <- morpho$SVL
names(svl) <- rownames(morpho)
pPC1 <- morpho$PC1
names(pPC1) <- rownames(morpho)
pPC2 <- morpho$PC2
names(pPC2) <- rownames(morpho)

pdf('plots/disparity/dtt.pdf', paper = 'a4', height = 10, width = 5)
par(mfrow=c(3,1))
dtt_svl <- dtt(tree, svl, nsim = 1000, calculateMDIp = TRUE)
mtext(text = 'SVL disparity through time', 
      side = 3, line = -2, cex = 1, font = 2)
dtt_pc1 <- dtt(tree, pPC1, nsim = 1000, calculateMDIp = TRUE)
mtext(text = 'PC1 disparity through time', 
      side = 3, line = -2, cex = 1, font = 2)
dtt_pc2 <- dtt(tree, pPC2, nsim = 1000, calculateMDIp = TRUE)
mtext(text = 'PC2 disparity through time', 
      side = 3, line = -2, cex = 1, font = 2)
dev.off()






