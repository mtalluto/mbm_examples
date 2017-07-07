library(raster)
taxModelMF <- readRDS('3_alps/res/taxModelMF.rds')
funModel <- readRDS('3_alps/res/funModel.rds')
phyModel <- readRDS('3_alps/res/phyModel.rds')

taxSpatial <- spatial_predict(taxModelMF, rasDat)
funSpatial <- spatial_predict(funModel, rasDat)
phySpatial <- spatial_predict(phyModel, rasDat)

saveMBMSP <- function(x, file)
{
	writeRaster(x[[1]], paste0(file, '_fit.grd'))
	writeRaster(x[[2]], paste0(file, '_sd.grd'))
	saveRDS(x[[3]], paste0(file, '_pca.rds'))
}

loadMBMSP <- function(file)
{
	x <- list()
	x$fits <- stack(paste0(file, '_fit.grd'))
	x$stdev <- raster(paste0(file, '_sd.grd'))
	x$pca <- readRDS(paste0(file, '_pca.rds'))
	class(x) <- c('mbmSP', class(x))
	x
}

saveMBMSP(taxSpatial, '3_alps/res/taxSpatial')
saveMBMSP(funSpatial, '3_alps/res/funSpatial')
saveMBMSP(phySpatial, '3_alps/res/phySpatial')


taxSpatial <- loadMBMSP('3_alps/res/taxSpatial')
funSpatial <- loadMBMSP('3_alps/res/funSpatial')
phySpatial <- loadMBMSP('3_alps/res/phySpatial')

quartz(w=12, h=6, type='png', file="3_alps/img/spatialPredict.png", dpi=300, bg='white')
par(mfrow=c(1,3), oma=c(0.5, 0.5, 5, 0.5))
raster::plotRGB(taxSpatial$fits[[c(1,3,2)]], scale=1, main="TD")
mtext("Taxonomic Similarity", side=3, line=2)
raster::plotRGB(funSpatial$fits[[c(1,3,2)]], scale=1, main="FD")
mtext("Functional Similarity", side=3, line=2)
raster::plotRGB(phySpatial$fits[[c(2,3,1)]], scale=1, main="PD")
mtext("Phylogenetic Similarity", side=3, line=2)
dev.off()
