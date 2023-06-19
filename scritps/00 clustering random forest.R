
# Load libraries ------------------------------------------------ 
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, fs, 
               randomForest, multcomp, dismo, ggpubr)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data -----------------------------------------------------
pnts <- read_csv('../tbl/presences/durangensis_rmOtl.csv')
smoc <- sf::st_read('E:/asesorias/SDM/shp/base/sierra_madre_occidental/SMO.shp')




