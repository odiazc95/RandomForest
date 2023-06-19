
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(tidyverse, raster, rgdal, gtools, cclust, dismo, gtools, sp, rgeos, FactoMineR, pROC, randomForest, Hmisc, rgeos)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

run <- 'run_1'
source('./FunctionsRFclustering.R')

# Functions ---------------------------------------------------------------
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  # occ = back_swd; nforest = 50; ntrees = 500; nVars = 8; nclasses = 2
  datRF_presences <- occ[,3:ncol(occ)] %>% as.data.frame()
  print(nrow(datRF_presences))
  attach(datRF_presences)
  no.forests <- nforest
  no.trees <- ntrees
  distRF_presences <- RFdist(datRF_presences, mtry1 = nVars, no.trees, no.forests, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  no.presencesclasses <- nclasses
  labelRF <- pamNew(distRF_presences$cl1, no.presencesclasses)
  print(table(labelRF))
  clusterdata <- hclust(as.dist(distRF_presences$cl1), method = 'ward.D2')
  return(list(labelRF, clusterdata))
}


# Load data ---------------------------------------------------------------
mask <- raster('../tif/worldclim/bio_1.tif') * 0 

load('../rData/run_1/clustereddata.rData')
pnts
vars <- colnames(pnts)
vars <- vars[3:length(vars)]

# Stacking -----------------------------------------------------------------
fles <- list.files('../tif/worldclim', full.names = TRUE, pattern = '.tif$') 
fles <- mixedsort(fles)
stck <- raster::stack(fles)

# Create background --------------------------------------------------------
clls <- raster::extract(mask, pnts[,1:2], cellnumbers = TRUE)
back <- mask
back[clls[,1]] <- NA
back <- rasterToPoints(back, spatial = FALSE)
back <- as_tibble(back)
back <- sample_n(tbl = back, size = nrow(pnts), replace = FALSE)

# A simple plot
plot(mask)
points(pnts$X, pnts$Y, pch = 16, col = 'red')
points(back$x, back$y, pch = 16, col = 'black')

# Extract the values for pseudo-absences ----------------------------------
back_swd <- raster::extract(stck, back[,1:2])
back_swd <- cbind(back[,1:2], back_swd)
write.csv(back_swd, '../tbl/presences/background_swd.csv', row.names = FALSE)

# Cluster analysis pseudo-absences ----------------------------------------
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)

datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50
no.trees <- 500 
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = TRUE, addcl2 = F, imp = F, oob.prox1 = TRUE)
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

presvalue_swd <- clusteredpresdata[,3:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)

presvalue_swd <- dplyr::select(presvalue_swd, pb, bio_1:bio_19)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) 

unique(presvalue_swd$pb)
unique(classdata_2$pb)

dim(classdata_2); dim(presvalue_swd)

allclasses_swd <- rbind(classdata_2, presvalue_swd)
table(allclasses_swd$pb)

write.csv(allclasses_swd, '../tbl/presences/allclasses_swd.csv', row.names = FALSE)

# To make the random forest model - Presences and pseudoabsences ----------
vrs <- vars
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+'), sep = ' '))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

samplesize <- min(table(presvalue_swd$pb))


for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  NumberOfClusters <- 3
  save(rfmodel, file = paste('../rf/output/run_1/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}


auc <- unlist(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

save(rflist, file = paste('../rData/run_1/', '/rflist_', NumberOfClusters, '.rdata', sep = ''))
save(importance, file = paste0('../rData/run_1/', '/importanceRF.rData'))
save(auc, file = paste0('../rData/run_1/', '/aucRF_dist.rData'))
save(rff, file = paste0('../rData/run_1/', '/rff_dist.rData'))

# Predict modell
climatevalues  <- data.frame(getValues(stck))
NumberOfClusters <- 3

rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- stck[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- stck[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- stck[[1]]
values(rasterRFclass) <- rasterRF

writeRaster(rasterRFclass, paste0('../rf/output/run_1/results/raw/RF_5Clust_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFprob, paste0('../rf/output/run_1/results/raw/RF_5Prob_current.asc'), format = 'ascii', overwrite = T)
writeRaster(rasterRFuncertainty, paste0('../rf/output/run_1/results/raw/RF_5Unc_current.asc'), format = 'ascii', overwrite = T)

plot(rasterRFclass)
plot(rasterRFprob)
plot(rasterRFuncertainty)

plot(1 - rasterRFuncertainty)


pnts <- as_tibble(cbind(raster::extract(rasterRFclass, pnts[,1:2]), pnts))
library(glue)
colnames(pnts) <- c('pb', 'x', 'y', glue('bio_{1:19}'))

# To make boxplot
tble <- gather(data = pnts, variable, valor, -pb, -x, -y)
tble <- mutate(tble, variable = factor(variable, levels = paste0('bio_', 1:19)))
tble <- tble %>% filter(!pb %in% 1:2)
tble <- mutate(tble, pb = factor(pb, levels = c('3', '4', '5')))

tble %>% filter(is.na(pb))

gplot <- ggplot(data = tble, aes(x = pb, y = valor)) + 
  geom_boxplot() + 
  facet_wrap(~variable, scales = 'free_y')

dir.create('../png')
ggsave(plot = gplot, filename = '../png/boxplot_clases.png', units = 'in', width = 12, height = 10, dpi = 300)

# Tipo 1 Fresco - Seco ##4ACD9F
# Tipo 2 Muy calido - Humedo #8E9656
# Tipo 3 Muy fresco - Muy humedo #227EC5




