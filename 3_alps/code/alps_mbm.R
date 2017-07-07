# devtools::install_local("~/work/projects/mbm/mbm")

library("mbm")
library("mbmdata")
library("ape")
# library(reshape2)
data(alps)
library(raster)


# # site selection
# # keep all sites above 3000m, otherwise keep only sites with >5 releves
keep <- (alps$siteEnv[,'num_releves'] >= 5 | alps$siteEnv[,'elev'] > 3000)

# choose a maximum of 9 sites from 250m elevational bands
get_rows <- function(minElev, maxElev, elev, n, mask) {
	inds <- which(elev > minElev & elev <= maxElev & mask)
	if(length(inds) <= n)
	{
		return(inds)
	} else {
		return(sample(inds, n))
	}
}

elevBands <- seq(0,3500,250)
maxN <- ceiling(100 / (length(elevBands) - 1)) # 100 is the approx target final N; --> 9 per band

# 
# # update, rows have been chosen, so we just load them
# if(file.exists("3_alps/res/selectedRows.rds"))
# {
	rows <- readRDS("3_alps/res/selectedRows.rds")
# } else {
# 	# choose rows
# 	rows <- unlist(mapply(get_rows, minElev=elevBands[1:(length(elevBands)-1)], maxElev=elevBands[2:length(elevBands)], 
# 						  MoreArgs=list(elev=alps$siteEnv[,'elev'], n=maxN, mask=keep)))
# 	# validation rows
# 	rows <- list(fit = rows, valid = sample((1:nrow(alps$siteEnv))[-rows], length(rows)))
# 	saveRDS(rows, '3_alps/res/selectedRows.rds')
# }

#for testing
# rows$fit <- rows$fit[1:40]


envVars <- c('bio_4', 'bio_6', 'bio_7', 'bio_15')
envMat <- lapply(rows, function(r) as.matrix(alps$siteEnv[r,envVars]))
taxBeta <- lapply(rows, function(r) sorensen(alps$siteSpecies[r,]))
# saveRDS(taxBeta, "3_alps/res/gdmdat/taxBeta.rds")
# saveRDS(envMat, "3_alps/res/gdmdat/envMat.rds")
# saveRDS(lapply(rows, function(r) as.matrix(alps$siteEnv[r,c('x', 'y')])), "3_alps/res/gdmdat/coords.rds")

taxBeta <- readRDS("3_alps/res/gdmdat/taxBeta.rds")
envMat <- readRDS("3_alps/res/gdmdat/envMat.rds")
# functional MPD
trMat <- alps$spTrait
trMat <- trMat[,-2] # drop repro ht
trMat <- trMat[complete.cases(trMat),]
# trMat[,'SLA'] <- log(trMat[,'SLA'])
# trMat[,'SEEDM'] <- log(trMat[,'SEEDM'])
# trMat[,'PL_VEG_H'] <- log(trMat[,'PL_VEG_H'])
# center/scale traits before computing distance, then scale distance between 0 and 1
# trMat <- scale(trMat)
# trDis <- as.matrix(dist(trMat))
# trDis <- trDis / max(trDis)
# trMPD <- lapply(rows, function(r) mpd(alps$siteSpecies[r,], dis=trDis))
# saveRDS(trMPD, "3_alps/res/gdmdat/funBeta.rds")
trMPD <- readRDS("3_alps/res/gdmdat/funBeta.rds")

# # phylo MPD
# phDis <- cophenetic(alps$phylogeny)
# phMPD <- lapply(rows, function(r) mpd(alps$siteGenus[r,], dis=phDis))
# saveRDS(phMPD, "3_alps/res/gdmdat/phyBeta.rds")
phMPD <- readRDS("3_alps/res/gdmdat/phyBeta.rds")

# prep the raster dataset
rasDat <- as.data.frame(alps$envRaster)
coordinates(rasDat) <- c('x','y')
gridded(rasDat) <- TRUE
rasDat <- stack(rasDat)
proj4string(rasDat) <- attr(alps$envRaster, "proj4string")
alpMask <- raster("~/work/projects/mbm/mbmdata/data-raw/alps/alpsMask.img")
alpMask <- projectRaster(alpMask, rasDat)
rasDat <- mask(rasDat, alpMask)

names <- rasDat[[1]]
names[!is.na(names)] <- 1:sum(!is.na(getValues(names)))
rasDat$names <- names

taxModelMF <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo_mf', 
				  force_increasing = TRUE)
saveRDS(taxModelMF, '3_alps/res/taxModelMF.rds')


funModel <- mbm(trMPD$fit, envMat$fit, predictX = envMat['valid'], link='identity', response_curve = 'distance', y_name = 'functional')
saveRDS(funModel, '3_alps/res/funModel.rds')
# GDM is a bit different; could try both
# 1. perform hierarchical clustering on training cells
# 2. compute predicted dissimilarity for ALL cells (including training)
# 3. Allocate all cells to a class according to what training point is most similar
# 4. Compute a class-class dissimilarity matrix
# 5. perform an ordination on this matrix & proceed as above

phyModel <- mbm(phMPD$fit, envMat$fit, predictX = envMat['valid'], link='identity', response_curve = 'distance', y_name = 'phylogenetic')
saveRDS(phyModel, '3_alps/res/phyModel.rds')


taxSpatial <- spatial_predict(taxModelMF, rasDat))
funSpatial <- spatial_predict(funModel, rasDat))
phySpatial <- spatial_predict(phyModel, rasDat))

# 	}, error=function(e) print(e))
# tryCatch({
# 	}, error=function(e) print(e))
# tryCatch({
# taxModel <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo')
# 		saveRDS(taxModel, '3_alps/res/taxModel.rds')
# 	}, error=function(e) print(e))
# tryCatch({
# 		taxModel2 <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo_2',
# 				 lengthscale = c(round(taxModelMF$params['lengthscale.distance']*2, 1), rep(NA, ncol(envMat$fit))), force_increasing = TRUE, pyMsg = TRUE)
# 		saveRDS(taxModel2, '3_alps/res/taxModel2.rds')
# 	}, error=function(e) print(e))
# tryCatch({
# 		taxModel0.2 <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo_0.2',
# 			   lengthscale = c(round(taxModelMF$params['lengthscale.distance']*0.2, 1), rep(NA, ncol(envMat$fit))), force_increasing = TRUE)
# 		saveRDS(taxModel0.2, '3_alps/res/taxModel0_2.rds')
# 	}, error=function(e) print(e))












