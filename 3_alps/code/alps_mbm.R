# devtools::install_local("~/work/projects/mbm/mbm")

library("mbm")
library("mbmdata")
library("ape")
# library(reshape2)
data(alps)


# site selection
# keep all sites above 3000m, otherwise keep only sites with >5 releves
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
# choose rows
rows <- unlist(mapply(get_rows, minElev=elevBands[1:(length(elevBands)-1)], maxElev=elevBands[2:length(elevBands)], 
					  MoreArgs=list(elev=alps$siteEnv[,'elev'], n=maxN, mask=keep)))

## FOR TESTING: just 15 sites
# rows <- sample(rows, 15)

# validation rows
rows <- list(fit = rows, valid = sample((1:nrow(alps$siteEnv))[-rows], length(rows)))

envVars <- c('bio_4', 'bio_6', 'bio_7', 'bio_15')
envMat <- lapply(rows, function(r) as.matrix(alps$siteEnv[r,envVars]))
taxBeta <- lapply(rows, function(r) sorensen(alps$siteSpecies[r,]))

# functional MPD
trMat <- alps$spTrait
trMat <- trMat[,-2] # drop repro ht
trMat <- trMat[complete.cases(trMat),]
trMat[,'SLA'] <- log(trMat[,'SLA'])
trMat[,'SEEDM'] <- log(trMat[,'SEEDM'])
trMat[,'PL_VEG_H'] <- log(trMat[,'PL_VEG_H'])
# center/scale traits before computing distance, then scale distance between 0 and 1
trMat <- scale(trMat)
trDis <- as.matrix(dist(trMat))
trDis <- trDis / max(trDis)
trMPD <- lapply(rows, function(r) mpd(alps$siteSpecies[r,], dis=trDis))

# phylo MPD
phDis <- cophenetic(alps$phylogeny)
phMPD <- lapply(rows, function(r) mpd(alps$siteGenus[r,], dis=phDis))


tryCatch({
		funModel <- mbm(trMPD$fit, envMat$fit, predictX = envMat['valid'], link='identity', response_curve = 'distance', y_name = 'functional')
		saveRDS(funModel, '3_alps/res/funModel.rds')
	}, error=function(e) print(e))
tryCatch({
		phyModel <- mbm(phMPD$fit, envMat$fit, predictX = envMat['valid'], link='identity', response_curve = 'distance', y_name = 'phylogenetic')
		saveRDS(phyModel, '3_alps/res/phyModel.rds')
	}, error=function(e) print(e))
tryCatch({
taxModel <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo')
		saveRDS(taxModel, '3_alps/res/taxModel.rds')
	}, error=function(e) print(e))
tryCatch({
		taxModel2 <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo_2',
				 lengthscale = c(round(taxModel$params[2]*2, 1), rep(NA, ncol(envMat$fit))))
		saveRDS(taxModel2, '3_alps/res/taxModel2.rds')
	}, error=function(e) print(e))
tryCatch({
		taxModel0.2 <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo_0.2',
			   lengthscale = c(round(taxModel$params[2]*0.2, 1), rep(NA, ncol(envMat$fit))))
		saveRDS(taxModel0.2, '3_alps/res/taxModel0_2.rds')
	}, error=function(e) print(e))
tryCatch({
	taxModelMF <- mbm(taxBeta$fit, envMat$fit, predictX = envMat['valid'], link='probit', response_curve = 'distance', y_name = 'taxo_mf', 
					   mean_func='increasing')
	saveRDS(taxModelMF, '3_alps/res/taxModelMF.rds')
}, error=function(e) print(e))












