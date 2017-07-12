library(raster)
library(mbmdata)
library(mbm)
data(alps)

args <- commandArgs(trailingOnly = TRUE)

cat('preparing raster dataset\n')
rasDat <- as.data.frame(alps$envRaster)
coordinates(rasDat) <- c('x', 'y')
proj4string(rasDat) <- attr(alps$envRaster, 'proj4string')
gridded(rasDat) <- TRUE
rasDat <- stack(rasDat)
alpMask <- rasterFromXYZ(alps$alpMask)
proj4string(alpMask) <- attr(alps$alpMask, 'proj4string')
rasDat <- mask(rasDat, alpMask)

# for testing
# rasDatSm <- crop(rasDat, extent(c(-350000, -300000, 1600000, 1650000)))

saveMBMSP <- function(x, file)
{
	writeRaster(x[[1]], paste0(file, '_fit.grd'))
	writeRaster(x[[2]], paste0(file, '_sd.grd'))
	saveRDS(x[[3]], paste0(file, '_pca.rds'))
}

if('-t' %in% args || any(grepl('tax', args)))
{
	cat('doing taxo spatial prediction\n')
	taxModelMF <- readRDS('3_alps/res/taxModelMF.rds')
	taxSpatial <- spatial_predict(taxModelMF, rasDat)
	saveMBMSP(taxSpatial, '3_alps/res/spatial/taxSpatial')
}
if('-f' %in% args || any(grepl('func', args)))
{
	cat('doing functional spatial prediction\n')
	funModel <- readRDS('3_alps/res/funModel.rds')
	funSpatial <- spatial_predict(funModel, rasDat)
	saveMBMSP(funSpatial, '3_alps/res/spatial/funSpatial')
}
if('-p' %in% args || any(grepl('phylo', args)))
{
	cat('doing phylo spatial prediction\n')
	phyModel <- readRDS('3_alps/res/phyModel.rds')
	phySpatial <- spatial_predict(phyModel, rasDat)
	saveMBMSP(phySpatial, '3_alps/res/spatial/phySpatial')
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


if(any(grepl('draw', args)))
{
	taxSpatial <- loadMBMSP('3_alps/res/spatial/taxSpatial')
	funSpatial <- loadMBMSP('3_alps/res/spatial/funSpatial')
	phySpatial <- loadMBMSP('3_alps/res/spatial/phySpatial')

	xl <- extent(taxSpatial$fits)[1:2]
	yl <- extent(taxSpatial$fits)[3:4]
	quartz(w=9, h=6, type='png', file="3_alps/img/spatialPredict.png", dpi=300, bg='white')
	par(mfrow=c(2,3), oma=c(1,0,3,2), mar=c(0,0,2,1))

	raster::plot(taxSpatial$fits[[1]], box=FALSE, axes=FALSE, alpha=0, legend=FALSE)
	raster::plotRGB(taxSpatial$fits[[c(1,3,2)]], scale=1, add=TRUE)
	mtext("A. Taxonomic Diversity", side=3, line=2)

	raster::plot(taxSpatial$fits[[1]], box=FALSE, axes=FALSE, alpha=0, legend=FALSE)
	raster::plotRGB(funSpatial$fits[[c(1,3,2)]], scale=1, add=TRUE)
	mtext("B. Functional Diversity", side=3, line=2)

	raster::plot(taxSpatial$fits[[1]], box=FALSE, axes=FALSE, alpha=0, legend=FALSE)
	raster::plotRGB(phySpatial$fits[[c(2,3,1)]], scale=1, add=TRUE)
	mtext("C. Phylogenetic Diversity", side=3, line=2)

	# plot(xl,yl, type='n', axes=F, xlab='', ylab='')
	raster::plot(taxSpatial$stdev, col=heat.colors(100), box=FALSE, axes=FALSE)

	# plot(xl,yl, type='n', axes=F, xlab='', ylab='')
	raster::plot(funSpatial$stdev, col=heat.colors(100), box=FALSE, axes=FALSE)

	# plot(xl,yl, type='n', axes=F, xlab='', ylab='')
	raster::plot(phySpatial$stdev, col=heat.colors(100), box=FALSE, axes=FALSE)

	dev.off()	
}

