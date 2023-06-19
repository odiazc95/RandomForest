

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, sf, tidyverse, gtools)

g <- gc(reset = TRUE)
rm(list = ls())

# Load data ---------------------------------------------------------------
clst <- raster('../rf/output/run_1/results/raw/RF_limitations.tif')
uncr <- raster('../rf/output/run_1/results/raw/RF_5Unc_current.asc')
prob <- raster('../rf/output/run_1/results/raw/RF_5Prob_current.asc')

load(file = '../rData/run_1/clustereddata.rData')
thrs.prob <- 0.798712

# Get the values for the uncertainty --------------------------------------
vles <- raster::extract(uncr, pnts[,1:2])
is.na(vles) %>% table()
vles <- quantile(x = vles, seq(0, 1, 0.01))
vles <- as.data.frame(vles)
vles
vles <- rownames_to_column(vles)
colnames(vles) <- c('porcentaje', 'valor')
glimpse(vles)
vles$porcentaje <- gsub('%', '', vles$porcentaje)
vles$porcentaje <- as.numeric(vles$porcentaje)

# Graphic with ggplot2 
ggplot(data = vles, aes(x = porcentaje, y = valor)) +
  geom_point()
  
thrs <- filter(vles, porcentaje == 10)
thrs <- pull(thrs, 2)

saveRDS(object = thrs, file = '../rData/run_1/threshold_unc.rds')

# Make the reclassify 
plot(uncr)
uncr.binr <- uncr
uncr.binr <- reclassify(x = uncr.binr, rcl = c(0, thrs, 0, thrs, 1, 1))

# A simple plot (both, raw and binary)
par(mfrow = c(1, 2))
plot(uncr, main = 'Raw')
plot(uncr.binr, main = 'Binary')
par(mfrow = c(1, 1))

# To make the result 
rslt <- clst
rslt[which(uncr[] < thrs & prob[] > thrs.prob)] <- 7

writeRaster(rslt, filename = '../rf/output/run_1/results/raw/RF_mixed.tif')

# Make a simple map
tble <- rasterToPoints(x = rslt)
tble <- as_tibble(tble)

lbls <- data.frame(RF_limitations = 1:7, 
                   clase = c('No idoneo', 'No idoneo', 'Idoneo tipo 1', 'Idoneo tipo 2', 'Idoneo tipo 3', 'Limitaciones', 'Aptitud incierta'))
tble <- inner_join(tble, lbls,  by = 'RF_limitations')
freq <- table(tble$clase) %>% 
  as.data.frame() %>% 
  mutate(porc = Freq / sum(Freq) * 100, 
         porc = round(porc, 1)) %>% 
  arrange(desc(porc))

tble <- mutate(
  tble, 
  clase = factor(clase, levels = c('No idoneo', 'Idoneo tipo 1', 'Idoneo tipo 2', 'Idoneo tipo 3', 'Limitaciones', 'Aptitud incierta'))
)

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x, y, fill = clase)) +
  scale_fill_manual(values = c('white', 'violetred', 'steelblue', 'palegreen3', 'grey', 'yellow')) + 
  coord_sf() + 
  theme(legend.position = 'bottom') +
  labs(fill = '', x = 'Longitud', y = 'Latitud')

ggsave(plot = gmap, 
       filename = '../png/map_final_bsl.png', 
       units = 'in', width = 10, height = 8, dpi = 300)
