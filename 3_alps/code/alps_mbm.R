devtools::install_local("~/work/projects/mbm/mbmtools")

library("mbmtools")
library("mbmdata")
# library("ape")
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

## FOR TESTING: just 10 sites
rows <- rows[1:10]

envVars <- c('bio_4', 'bio_6', 'bio_7', 'bio_15')
envMat <- as.matrix(alps$siteEnv[rows,envVars])
taxBeta <- sorensen(alps$siteSpecies[rows,])

model <- mbm(taxBeta, envMat)
