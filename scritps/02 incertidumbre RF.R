
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(tidyverse, raster, rgdal, gtools, cclust, dismo, gtools, sp, rgeos, FactoMineR, pROC, randomForest, Hmisc, rgeos)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
load(file = '../rData/run_1/clustereddata.rData')

pnts
clusteredpresdata

prob <- raster('../rf/output/run_1/results/raw/RF_5Prob_current.asc')
plot(prob)
clst <- raster('../rf/output/run_1/results/raw/RF_5Clust_current.asc')

vlrs.prob <- raster::extract(prob, pnts[,1:2])
vlrs.prob
qntl <- quantile(x = vlrs.prob, seq(0, 1, 0.01))
qntl

qntl <- rownames_to_column(as.data.frame(qntl))
qntl <- mutate(qntl, rowname = parse_number(rowname))

# A simple plot 
ggplot(data = qntl, aes(x = rowname, y = qntl)) + 
  geom_point()

thrs <- filter(qntl, rowname == 10)
thrs <- pull(thrs, 2)

# To make a binary raster
prob.binr <- raster::reclassify(prob, c(0, thrs, 0, thrs, 1, 1))
plot(prob.binr)

# Binary raster cluster 
clst.binr <- raster::reclassify(clst, c(0.5, 2.5, 0, 2.5, 5.5, 1))
prob.binr <- raster::reclassify(prob, c(0, thrs, 0, thrs, 1, 2))

par(mfrow = c(1, 2))
plot(clst.binr, main = 'Cluster')
plot(prob.binr, main = 'Probability')
par(mfrow = c(1, 1))

diff <- prob.binr - clst.binr
rslt <- clst
rslt[which(diff[] == -1)] <- 6
rslt[which(diff[] ==  2)] <- 6
plot(rslt)


library(terra)
plot(rast(rslt))

raster::crs(rslt) <- '+proj=longlat +datum=WGS84 +no_defs'
raster::writeRaster(rslt, filename = '../rf/output/run_1/results/raw/RF_limitations.tif')
