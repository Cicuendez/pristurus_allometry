
# LAND DISPARITY 
# Explanation: ----
# Is morphological disparity higher in the island species than in the mainland? 
# We compute the actual disparity in mainland and island, and we simulate datasets 
# by simulating the evolution of the trait, to see if the disparity difference observed 
# between mainland and island is different from expected. 
# We calculate the ratio disp_island/disp_mainland, for the actual data (empiric part) 
# and for the simulations. If the observed ratio is significantly higher than 
# expected based on the simulations, then we'll assume that disparity in the island is 
# effectively higher, and that there is an 'island effect' for the trait under consideration. 


# Load packages ----
libs <- c('phytools', 'geiger', 'tidyverse')
lapply(libs, require, character.only = TRUE)

# Import data ----
tree <- read.nexus('data/phylogeny/pristurus_tree_final.nex')
morpho <- read.table('objects/phypca/phypca_scores.csv', sep = ";", dec = '.', 
                     header = TRUE, row.names = 1)
name.check(tree, morpho)
land <- morpho$land
names(land) <- rownames(morpho)


# 1. SVL ----
## 1.1. Empiric part ---- 
### Split island and mainland ----
svl_mainland <- morpho$SVL[morpho$land=="mainland"]
names(svl_mainland) <- rownames(morpho)[morpho$land=="mainland"]

svl_island <- morpho$SVL[morpho$land=="island"]
names(svl_island) <- rownames(morpho)[morpho$land=="island"]

### Calculate disparity of island and mainland separately ----
dispSVL_island <- disparity(data=svl_island, index="avg.sq")
dispSVL_mainland <- disparity(data=svl_mainland, index="avg.sq")
?disparity

### Calculate disparity ratio island/mainland ----
dispSVL_ratio_emp <- dispSVL_island/dispSVL_mainland
dispSVL_ratio_emp
### If the value is greater than 1, the disparity is higher in the island 
### If the value is less than 1, the disparity is higher in mainland
### Ratio SVL disparity island/mainland = 1.98 -> Disparity is higher in island.
## END OF THE EMPIRIC PART

## 1.2. Simulation part ----
### Question: Is the difference between island and mainland disparity significant?
### Model of trait evolution ----
### Fit models of trait evolution on the tree
?fitContinuous
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



### Assign island or mainland to the species in each simulation ----
sim_land_svl <- list()
for (i in 1:nsim){
  sim_land_svl[[i]] <- cbind(morpho$land, sim_svl[,1,i])
}
## The result is a list of nsim datasets (one for each simulation) with the values of the trait 
## and the assignment to mainland (2) and island (1) ##

### Split island (1) and mainland (2) in each simulation ##
sim_island_svl <- list()
for (i in 1:length(sim_land_svl)){
  sim_island_svl[[i]] <- sim_land_svl[[i]][,2][sim_land_svl[[i]]==1]
}

sim_mainland_svl <- list()
for (i in 1:length(sim_land_svl)){
  sim_mainland_svl[[i]] <- sim_land_svl[[i]][,2][sim_land_svl[[i]]==2]
}
## The result is two lists of nsim simulated datasets, 
## one with island species and another one with mainland species, 
## and their simulated values of trait (svl)

## CALCULATE THE DISPARITY IN ISLAND AND MAINLAND SEPARATELY ##
disp_sim_island_svl <- list()
for (i in 1:length(sim_island_svl)){
  disp_sim_island_svl[[i]] <- disparity(data=sim_island_svl[[i]],
                                        index="avg.sq")
}

disp_sim_mainland_svl <- list()
for (i in 1:length(sim_mainland_svl)){
  disp_sim_mainland_svl[[i]] <- disparity(data=sim_mainland_svl[[i]],
                                          index="avg.sq")
}

## The result is two lists of values of disparity, one for island and one for mainland ##

### Calculate disparity ratio (island/mainland) for each pair of datasets ----
### (each pair of elements from the lists)
disp_ratio_sim_svl <- c()
for (i in 1:length(disp_sim_mainland_svl)){
  disp_ratio_sim_svl[i] <- disp_sim_island_svl[[i]]/disp_sim_mainland_svl[[i]]
}
disp_ratio_sim_svl
## The result is a vector with n values of disparity ratios ##
## End of simulation part

## 1.3. Histogram of ratio values ----
par(mar=c(5,5,2,2))
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
  xlim(0, 8)



## p-value for the empiric disparity ratio ----
## To do this, sum all the values of ratio equal or greater than the empirical ratio, 
## and divide by the number of simulations
pval_svl <- sum(disp_ratio_sim_svl >= dispSVL_ratio_emp)/nsim
pval_svl

# Plot disparity ratio SVL ----
hist_svl_land <- ggplot(simratios_svl) +
  geom_histogram(aes(x = ratio), fill = 'gray80', color = 'white') +
  xlim(0, 8) +
  geom_segment(mapping = aes(x = dispSVL_ratio_emp, y = 200, xend = dispSVL_ratio_emp, yend = 0), 
               arrow = arrow(angle=30, length = unit(0.2, "cm"), type = 'open'), 
               size = 0.3, color = 'red') +
  annotate(geom = 'text', x = 4, y = 1150, color = 'black', 
           label = 'Observed disparity ratio:') +
  annotate(geom = 'text', x = 4, y = 1000, color = 'black', 
           label = paste0('r = ', round(dispSVL_ratio_emp, digits = 2), 
                          ', p = ', round(pval_svl, digits = 2))) +
  xlab('Disparity ratio (r)') +
  ylab('Frequency') +
  ggtitle(label = 'Body size disparity\nisland vs mainland') +
  theme.clean()

saveRDS(hist_svl_land, 'objects/disparity/hist_svl_land.rds')

## If the value is less than 0.05, the difference between island and mainland 
## is significant (greater than expected at random circumstances)
## p-value = 0.5249 for SVL -> no significant differences. 



#### 2. pPC1 ####
## 2.1. EMPIRIC PART ----
## SPLIT ISLAND AND MAINLAND ##
pPC1_mainland <- morpho$PC1[morpho$land=="mainland"]
names(pPC1_mainland) <- rownames(morpho)[morpho$land=="mainland"]

pPC1_island <- morpho$PC1[morpho$land=="island"]
names(pPC1_island) <- rownames(morpho)[morpho$land=="island"]

## CALCULATE THE DISPARITY OF ISLAND AND MAINLAND SEPARATELY ##
disppPC1_island <- disparity(data=pPC1_island, index="avg.sq")
disppPC1_mainland <- disparity(data=pPC1_mainland, index="avg.sq")

## CALCULATE THE DISPARITY RATIO ISLAND/MAINLAND ##
disppPC1_ratio_emp <- disppPC1_island/disppPC1_mainland
disppPC1_ratio_emp
## If the value is greater than 1, the disparity is higher in the island ## 
## If the value is less than 1, the disparity is higher in mainland ##
## END OF THE EMPIRIC PART ##
## The value is not greater than 1, so the disparity for pPC1 is higher in mainland ##
## There's no need to do the simulation part then ##

#### 3. pPC2 ####
## 3.1. EMPIRIC PART ----
## SPLIT ISLAND AND MAINLAND ##
pPC2_mainland <- morpho$PC2[morpho$land=="mainland"]
names(pPC2_mainland) <- rownames(morpho)[morpho$land=="mainland"]

pPC2_island <- morpho$PC2[morpho$land=="island"]
names(pPC2_island) <- rownames(morpho)[morpho$land=="island"]

## CALCULATE THE DISPARITY OF ISLAND AND MAINLAND SEPARATELY ##
disppPC2_island <- disparity(data=pPC2_island, index="avg.sq")
disppPC2_mainland <- disparity(data=pPC2_mainland, index="avg.sq")

## CALCULATE THE DISPARITY RATIO ISLAND/MAINLAND ##
disppPC2_ratio_emp <- disppPC2_island/disppPC2_mainland
disppPC2_ratio_emp
## If the value is greater than 1, the disparity is higher in the island ## 
## If the value is less than 1, the disparity is higher in mainland ##
## END OF THE EMPIRIC PART ##
## The value is not greater than 1, so the disparity for pPC2 is higher in mainland ##
## There's no need to do the simulation part then ##

dispSVL_ratio_emp
disppPC1_ratio_emp
disppPC2_ratio_emp

