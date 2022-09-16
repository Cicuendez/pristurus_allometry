library(RRPP)

# Packages ----
#libs <- c('geomorph', 'RRPP')
#easypackages::libraries(libs)

# Morphological data ----
data <- read.table('data/morpho/morpho_pristurus.csv', sep = ';', dec = '.', header = TRUE)

# Check the difference between first calculating the per species mean and then log-transforming 
# the data, and the other way around (log-transform the original data and calculate the mean)

species <- data$species
morpho <- data[, 7:ncol(data)]
log_morpho <- log(morpho)

# Get the number of specimens per species
table(species)

# Per-species mean of original data
morpho.mean <- rowsum(morpho, species)/as.vector(table(species))
morpho.dist <- dist(morpho.mean) # distance between species

# Per-species mean of log-transformed data
logmorpho.mean <- rowsum(log_morpho, species)/as.vector(table(species))
logmorpho.dist <- dist(logmorpho.mean) # distance between species

# Correlation between the two
cor(morpho.dist, logmorpho.dist)


# ALLOMETRY ----

# Allometry grouped by species ----
svl <- log(data$SVL)
shape <- as.matrix(log(data[, 8:ncol(data)]))
species.fctr <- as.factor(species)
rdf.sp <- rrpp.data.frame(svl = svl, shape = shape, species = species.fctr)

# Multivariate linear model
fit.sp <- lm.rrpp(shape~svl*species, data = rdf.sp)
anova(fit.sp)

plot(fit.sp, predictor = rdf.sp$svl, type = 'regression', pch = 16, col = as.numeric(rdf.sp$species))

# Pairwise differences in the angle 
pw.sp <- pairwise(fit.sp, groups = rdf.sp$species, covariate = rdf.sp$svl)
pw.sp_df <- summary(pw.sp, type = 'VC',stat.table = FALSE)

#write.csv(df_pw$pairwise.tables$P, 'pw_pvalue.csv', quote = FALSE, sep = ';', dec = '.')
# Check species with more than 5 specimens ----

# Allometry grouped by habitat (3 states) ----
habitat <- data$habitat_broad
habitat.fctr <- as.factor(habitat)
rdf.hab <- rrpp.data.frame(svl = svl, shape = shape, habitat = habitat.fctr)

# Multivariate linear model
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf.hab)
anova(fit.hab)
fit.hab$LM$coefficients

##### DCA: regression coefficients
fit.coef <- fit.hab$LM$coefficients

rbind(fit.coef[1,], fit.coef[2,]) #ground
rbind(fit.coef[1,]+fit.coef[3,], fit.coef[2,]+fit.coef[5,]) #rock
rbind(fit.coef[1,]+fit.coef[4,], fit.coef[2,]+fit.coef[6,]) #tree 


hab.colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")
plot(1:3, 1:3, col = hab.colors, cex = 10, pch = 16)

plot(fit.hab, predictor = rdf.hab$svl, type = 'regression', pch = 16, 
     col = hab.colors[as.numeric(rdf.hab$habitat)])
legend('topleft', levels(habitat.fctr), pt.bg = hab.colors, pch = 22)

# Pairwise differences in the angle ----
pw.hab <- pairwise(fit.hab, groups = rdf.hab$habitat, covariate = rdf.hab$svl)
pw.hab_df <- summary(pw.hab, type = 'VC', stat.table = FALSE)

pw.hab_df








#######


res <- two.b.pls(rdf$svl, rdf$shape)
res$right.pls.vectors

# Expected isometry for nvar variables
nvar <- 8
1/sqrt(nvar)

hab_list <- vector('list', length = 3)
names(hab_list) <- levels(habitat_fctr)

size.group <- split(rdf.habitat$svl, rdf.habitat$habitat)
shape.tmp <- split(rdf.habitat$shape, rdf.habitat$habitat)

shape.group <- list()
names(shape.group) <- names(size.group)
for (i in 1:length(shape.tmp)){
  shape.group[[i]] <- matrix(shape.tmp[[i]], ncol = nvar, byrow = TRUE)
  colnames(shape.group[[i]]) <- colnames(rdf.habitat$shape)
}

length(shape.group$ground)
?split


habs <- Map(function(x, y) two.b.pls(x, y), size.group, shape.group)
data.frame(ground = habs$ground$right.pls.vectors, rock = habs$rock$right.pls.vectors, 
            tree = habs$tree$right.pls.vectors)

log_shape - 



res <- two.b.pls(rdf$svl, rdf$shape)
res$right.pls.vectors


















