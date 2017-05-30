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
rows.valid <- sample((1:nrow(alps$siteEnv))[-rows], length(rows))
rows.predict <- (1:nrow(alps$siteEnv))[-c(rows, rows.valid)]

# covariate prep
# include seaonality (bio4), min temp of coldest month (6), annual temp range (7), and precip seasonality(15)
envVars <- c('bio_4', 'bio_6', 'bio_7', 'bio_15')
envMat <- env_dissim(alps$siteEnv[rows,envVars])
envMat.v <- env_dissim(alps$siteEnv[rows.valid,envVars])
envMat.p <- env_dissim(alps$siteEnv[rows.predict,envVars])

spBeta <- sorensen(alps$siteSpecies[rows,])
spBeta <- melt(spBeta,varnames=c('site1', 'site2'), value.name = 'sor')
taxBeta <- merge(spBeta, envMat)
spBeta.v <- sorensen(alps$siteSpecies[rows.valid,])
spBeta.v <- melt(spBeta.v,varnames=c('site1', 'site2'), value.name = 'sor')
taxBeta.v <- merge(spBeta.v, envMat.v)


# functional MPD
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
trMPD <- mpd(alps$siteSpecies[rows,], dis=trDis)
trMPD.v <- mpd(alps$siteSpecies[rows.valid,], dis=trDis)


## now melt and merge with environment
trMPD_m <- melt(trMPD,varnames=c('site1', 'site2'), value.name = 'f_mpd')
funcBeta <- merge(taxBeta, trMPD_m, all.x=TRUE)

trMPD_m.v <- melt(trMPD.v,varnames=c('site1', 'site2'), value.name = 'f_mpd')
funcBeta.v <- merge(taxBeta.v, trMPD_m.v, all.x=TRUE)

## now phylo
phDis <- cophenetic(alps$phylogeny)
phMPD <- mpd(alps$siteGenus[rows,], dis=phDis)
phMPD_m <- melt(phMPD,varnames=c('site1', 'site2'), value.name = 'p_mpd')
betaDiv <- merge(funcBeta, phMPD_m, all.x=TRUE)

phMPD.v <- mpd(alps$siteGenus[rows.valid,], dis=phDis)
phMPD_m.v <- melt(phMPD.v,varnames=c('site1', 'site2'), value.name = 'p_mpd')
betaDiv.v <- merge(funcBeta.v, phMPD_m.v, all.x=TRUE)

# make a response curve dataset for prediction
rcLim <- range(betaDiv$distance)
rcX <- data.frame(distance = seq(rcLim[1], rcLim[2], length.out = 200))
for(xx in envVars)
	rcX[,xx] <- 0
write.csv(rcX, "3_alps/dat/betaRC.csv", row.names=FALSE)

saveRDS(betaDiv, "3_alps/dat/betaDiv.rds")
write.csv(betaDiv, "3_alps/dat/betaDiv.csv", row.names = FALSE)
saveRDS(betaDiv.v, "3_alps/dat/betaDiv_valid.rds")
write.csv(betaDiv.v, "3_alps/dat/betaDiv_valid.csv", row.names = FALSE)
# write.csv(envMat.p, "3_alps/dat/betaDiv_predict.csv", row.names = FALSE)

##### TO TRY
# √ all traits (drop repro height) - we have 1054 species if we do this
# √ try log transforming sla, height, maybe seed mass (all quite skewed)
# √ multiple site sets, just to make sure this one isn't weird
# try adding soil data (maybe do env dist + soil dist separately)
#√√√√ fix my dissimilarity metric - it's busted right now; check sorensen_pd and try to fix it
# √ phylo (obvs)

