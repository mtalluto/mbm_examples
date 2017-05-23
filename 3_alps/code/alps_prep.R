#!/usr/bin/env Rscript
# devtools::install_local("~/work/projects/mbm/mbmdata")
# devtools::install_local("~/work/projects/mbm/mbmtools")

library("mbmtools")
library("mbmdata")
library("ape")
library(reshape2)
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


# covariate prep
# include seaonality (bio4), min temp of coldest month (6), annual temp range (7), and precip seasonality(15)
envVars <- c('bio_4', 'bio_6', 'bio_7', 'bio_15')
envMat <- env_dissim(alps$siteEnv[rows,envVars])

spBeta <- sorensen(alps$siteSpecies[rows,])
spBeta <- melt(spBeta,varnames=c('site1', 'site2'), value.name = 'sor')
taxBeta <- merge(spBeta, envMat)


# functional MPD and dissimilarity
trMat <- alps$spTrait
trMat <- trMat[,-2] # drop repro ht
trMat <- trMat[complete.cases(trMat),]
# try a transformed version as well
trMat_untransformed <- trMat
trMat[,'SLA'] <- log(trMat[,'SLA'])
trMat[,'SEEDM'] <- log(trMat[,'SEEDM'])
trMat[,'PL_VEG_H'] <- log(trMat[,'PL_VEG_H'])
trMat <- scale(trMat)
trDis <- as.matrix(dist(trMat))
trDis <- trDis / max(trDis)
trSor <- sorensen(alps$siteSpecies[rows,], trDis) ## this is still unacceptably slow, AND it is yielding negative values unless the distances are scaled to be between 0 and 1
trMPD <- mpd(alps$siteSpecies[rows,], dis=trDis)

## now melt and merge with environment
trSor_m <- melt(trSor,varnames=c('site1', 'site2'), value.name = 'f_sor')
trMPD_m <- melt(trMPD,varnames=c('site1', 'site2'), value.name = 'f_mpd')
funcBeta <- merge(taxBeta, trSor_m, all.x=TRUE)
funcBeta <- merge(funcBeta, trMPD_m, all.x=TRUE)

## now phylo
phDis <- cophenetic(alps$phylogeny)
phDis_sc <- phDis / max(phDis)
phSor <- sorensen(alps$siteGenus[rows,], phDis_sc)
phMPD <- mpd(alps$siteGenus[rows,], dis=phDis)
phSor_m <- melt(phSor,varnames=c('site1', 'site2'), value.name = 'p_sor')
phMPD_m <- melt(phMPD,varnames=c('site1', 'site2'), value.name = 'p_mpd')
phBeta <- merge(funcBeta, phSor_m, all.x=TRUE)
phBeta <- merge(phBeta, phMPD_m, all.x=TRUE)


# now try individual traits
trIndSor <- lapply(colnames(trMat), function(tr) {
	tdat <- trMat[,tr,drop=FALSE]
	tdis <- as.matrix(dist(tdat))
	tdis <- tdis/max(tdis)
	tsor <- sorensen(alps$siteSpecies[rows,], tdis)
	tmpd <- mpd(alps$siteSpecies[rows,], dis=tdis)
	tsor <- melt(tsor,varnames=c('site1', 'site2'), value.name = paste0(tr,'_sor'))
	merge(tsor, melt(tmpd,varnames=c('site1', 'site2'), value.name = paste0(tr,'_mpd')))
})
trIndSor <- Reduce(merge, trIndSor)
betaDiv <- merge(trIndSor, phBeta, all.y=TRUE)
saveRDS(betaDiv, "3_alps/dat/betaDiv.rds")
write.csv(betaDiv, "3_alps/dat/betaDiv.csv", row.names = FALSE)

##### TO TRY
# √ all traits (drop repro height) - we have 1054 species if we do this
# √ try log transforming sla, height, maybe seed mass (all quite skewed)
# √ multiple site sets, just to make sure this one isn't weird
# try adding soil data (maybe do env dist + soil dist separately)
#√√√√ fix my dissimilarity metric - it's busted right now; check sorensen_pd and try to fix it
# √ phylo (obvs)

