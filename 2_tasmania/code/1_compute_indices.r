#!/usr/bin/env Rscript

# if not installed do it
# devtools::install_github("mtalluto/mbmtools", ref="develop")
# devtools::install_bitbucket("matt_talluto/mbmdata", username="matt_talluto", pw=readRDS("~/.bbpw.rds"))
## or if testing
# mbmdir <- "~/work/projects/mbm/mbmtools"
# devtools::document(mbmdir)
# devtools::install_local(mbmdir)
# devtools::install_local("~/work/projects/mbm/mbmdata")


library("mbmtools")
library("mbmdata")
library("geosphere")
library("reshape2")

data(tasmania)

# get distance matrix
geogDistance <- distm(tasmania$sites[,c('longitude', 'latitude')], fun=distVincentyEllipsoid)
rownames(geogDistance) <- colnames(geogDistance) <- tasmania$sites$site
# scale distances to be in 100s of km
geogDistance <- geogDistance/100000
geogDistance <- melt(geogDistance, varnames = c('site1', 'site2'), value.name = 'geogDistance100km')

# set up environmental variables
covarNames <- c('precip_PET_ratio', 'Temp_Isothermallity', 'Temp_ColdestPeriod_min', 'Rad_Jan_total')
covars <- tasmania$env[,covarNames]
covars <- scale(covars)
saveRDS(covars, "2_tasmania/dat/covars_scaled.rds")
cvDissim <- env_dissim(covars, scale=FALSE)
# covarDissim.list <- lapply(covarNames, function(nm) env_dissim(covars[,nm,drop=FALSE], scale=FALSE))
# covarDissim <- do.call(cbind, covarDissim.list)
# covarDissim <- covarDissim[,c(3,4,1,2,5,6,9,10,13,14)]
# colnames(covarDissim) <- c('site1', 'site2', 'dist_pToPET', 'pToPET', 'dist_isothermality', 
# 	'isothermality', 'dist_tempColdest', 'tempColdest', 'dist_radJan', 'radJan')

covarDissim <- merge(cvDissim, geogDistance, by=c('site1', 'site2'), all.x=TRUE, all.y=FALSE)

# compute taxonomic beta diversity
taxBeta <- sorensen(tasmania$species)
taxBeta <- melt(taxBeta, varnames = c('site1', 'site2'), value.name = 'taxSorensen')
# transform to deal with the ones:
taxBeta$taxSorensen <- taxBeta$taxSorensen - 1e-4

funBeta <- melt(tasmania$fun.dissim, varnames = c('site1', 'site2'), value.name = 'functionalDistance')

tasData <- merge(taxBeta, funBeta, by=c('site1', 'site2'))
tasData <- merge(tasData, covarDissim, by=c('site1', 'site2'), all.x=FALSE, all.y=TRUE)
saveRDS(tasData, "tasData_all.rds")

## withhold 20% for validation
validProp <- 0.2
validRows <- sample(nrow(tasData), as.integer(validProp*nrow(tasData)))
write.table(tasData[validRows,], '2_tasmania/dat/tasValid.csv', sep=',', row.names=FALSE)

# fit to random 20% subsets
numFits <- 5
fitProp <- 0.2
tasFit <- tasData[-validRows,]
fitRows <- lapply(1:numFits, function(x) sample(nrow(tasFit), as.integer(fitProp*nrow(tasFit))))
lapply(1:numFits, function(x) {
	write.table(tasFit[fitRows[[x]], ], paste0('2_tasmania/dat/tasFit_', x, '.csv'), sep=',', row.names=FALSE)
})

