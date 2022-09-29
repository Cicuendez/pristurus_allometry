# R-script for primary analyses

libs <- c('geomorph', 'RRPP', 'phytools', 'geiger', 'tidyverse')
easypackages::libraries(libs)

# 0: Data Prep
data0 <- read.table('data/morpho/morpho_pristurus.csv', sep = ';', dec = '.', header = TRUE, stringsAsFactors = TRUE)
  sp.to.keep <- names(which(table(data0$species) >= 5))
data <- data0[data0$species %in% sp.to.keep, ]
  data$species <- droplevels(data$species)
  data$SVL <- log(data$SVL)
shape <- as.matrix(log(data[, 8:ncol(data)]))
rdf <- rrpp.data.frame(svl = data$SVL, shape = shape, habitat = data$habitat_broad, species = data$species)

# 2: Comparison of multivariate allometry among habitat types
fit.hab <- lm.rrpp(shape~svl*habitat, data = rdf)
  anova(fit.hab)
pw.hab <- pairwise(fit.hab, groups = rdf$habitat, covariate = rdf$svl)
  summary(pw.hab, type = 'VC', stat.table = FALSE)

  ######################## NEW STUFF IS IN HERE
mn.sz <- mean(rdf$svl)
mn.shape <- apply(rdf$shape,2,mean)
slopes <- rep(1,ncol(rdf$shape))
intercepts <- mn.shape - (slopes * mn.sz)
 X <- cbind(1,rdf$svl)
 b <- rbind(intercepts,slopes)
preds <- X%*%b
iso.line <- prcomp(preds)$x[,1]
plot(rdf$svl,iso.line)  #HERE IT IS

plot(fit.hab,type = "regression", predictor = rdf$svl)
 lines(rdf$svl,iso.line, col = "red")
  ######################## 
