library("mbmdata")
library(RColorBrewer)
library(rgdal)
library(raster)
# pal <- brewer.pal(3, "Dark2")
# ff <- quartzFont(c("Raleway Light", "Raleway Bold", "Raleway Light Italic", "Raleway Bold Italic"))

pal <- brewer.pal(3, "Set1")
proj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"

rows <- readRDS("3_alps/res/selectedRows.rds")
data(alps)
rasDat <- cbind(alps$siteEnv[,c('x', 'y')],0)
rasDat[rows$fit,3] <- 1
rasDat[rows$valid,3] <- 2
ras <- rasterFromXYZ(rasDat)
proj4string(ras) <- proj
alpMask <- projectRaster(raster("../mbmdata/data-raw/alps/alpsMask.img"), ras)

bnd <- readOGR("3_alps/dat/ne_10m_admin_0_countries", layer="ne_10m_admin_0_countries")
natE <- stack("3_alps/dat/NE2_HR_LC_SR_W_DR/NE2_HR_LC_SR_W_DR.tif")
cropE <- extent(c(-10, 20, 35,56))
bndC <- crop(bnd, cropE)
natEC <- crop(natE, cropE)
bndP <- spTransform(bndC, proj4string(ras))
natEP <- projectRaster(natEC, crs=proj4string(ras))
crop2 <- extent(c(-1346391.6, 624983.5, 671160.1, 2555926.6))

tp <- 'png'
quartz(width=6.5, height=4, pointsize=11, type=tp, file=paste0("3_alps/img/map.", tp), bg='white', dpi=600)
par(mfrow=c(1,2), bty='n', mar=c(4,4,4,0))
plot(ras, col=c('#dddddd', pal[2:3]), legend=FALSE, xaxt='n', yaxt='n', bty='n')
legend('topleft', pch=15, col=pal[2:3], legend=c("Calibration", "Validation"), cex=0.8)

plotRGB(crop(natEP, crop2))
plot(crop(bndP, crop2), add=TRUE, lwd=0.5)
plot(alpMask, col=paste0(pal[1], 'aa'), add=T, legend=FALSE)
dev.off()