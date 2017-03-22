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
library("reshape2")
library("ape")

data(guisane)

# unpack the list
coverSpecies <- guisane$species
coverGenus <- guisane$genera
phylogeny <- guisane$phylogenies[[guisane$bestPhylogeny]]
speciesTraits <- guisane$traits
covars <- guisane$env

# drop covariates we aren't using, then scale and save scaling
xycoords <- covars[,c('x', 'y')]
covars <- covars[,-c(1:2,23:28,31)]
covars <- scale(covars)
cvscale <- cbind(center=attr(covars, 'scaled:center'), scale=attr(covars, 'scaled:scale'))
saveRDS(cvscale, "1_guisane/dat/covariateScaling.rds")
covars <- cbind(covars, xycoords)

## taxonomic             
alphaT <- data.frame(simpson = simpson(coverSpecies), richness=richness(coverSpecies))
betaT <- sorensen(coverSpecies)
betaT <- melt(betaT, varnames=c('site1', 'site2'), value.name = 'sorensen_t')

## phylogenetic
## drop genera that are missing from the phylogeny; convert to PA; also drop nonflowering plants
dropGenera <- c('Selaginella', 'Botrychium', 'Equisetum', 'Polystichum', 'Dryopteris')
coverGenus <- coverGenus[,which(colnames(coverGenus) %in% phylogeny$tip.label)]
coverGenus <- coverGenus[,-which(colnames(coverGenus) %in% dropGenera)]
coverGenus <- 1.0 * (coverGenus > 0)

# compute MPD matrix
phyDivMat <- mbmtools::mpd(comm = coverGenus, phylogeny=phylogeny, dis.transform=sqrt)
phyloAlpha <- data.frame(sites = rownames(coverGenus), mpd_ap= diag(phyDivMat))
diag(phyDivMat) <- NA
phyDivMat[upper.tri(phyDivMat)] <- NA
phyloBeta <- melt(phyDivMat, varnames=c('site1', 'site2'), value.name = 'mpd_bp')
drop <- which(is.na(phyloBeta$mpd_bp))
phyloBeta <- phyloBeta[-drop,]


## functional
# clean up the trait data a bit
# cor(speciesTraits$h_rep, speciesTraits$h_veg, use="complete.obs") ## 0.864
# c(sum(is.na(speciesTraits$h_rep)), sum(is.na(speciesTraits$h_veg)))
# we are going to drop repro height, it has more NAs
speciesTraits <- speciesTraits[,-3]

# drop species that have 1 or 0 traits
dropSp <- which(rowSums(is.na(speciesTraits[3:6])) >= 3)
speciesTraits <- speciesTraits[-dropSp,]

traits <- melt(as.matrix(speciesTraits[,c(3:6)]), na.rm=TRUE, varnames =c('species', 'variable'))
funcDendro <- dendrogram(traits, data.frame(unique(traits$variable), 'Q'), species=1, traits=2, values=3)

# make the sp matrix into PA and drop species not in our trait database
fPASpecies <- 1.0 * (coverSpecies > 0)
fPASpecies <- fPASpecies[,which(colnames(fPASpecies) %in% traits$species)]
funDivMat <- mbmtools::mpd(comm = fPASpecies, phylogeny=funcDendro$dendrogram)
funcAlpha <- data.frame(sites=rownames(coverSpecies), mpd_af = diag(funDivMat))
diag(funDivMat) <- NA
funDivMat[upper.tri(funDivMat)] <- NA
funcBeta <- melt(funDivMat, varnames=c('site1', 'site2'), value.name = 'mpd_bf')
drop <- which(is.na(funcBeta$mpd_bf))
funcBeta <- funcBeta[-drop,]



## combine alpha and beta into single dataframes
alpha <- merge(alphaT, phyloAlpha, by.x=0, by.y='sites')
alpha <- merge(alpha, funcAlpha, by.x=1, by.y='sites')
beta <- merge(betaT, phyloBeta)
beta <- merge(beta, funcBeta)


## covariates
alpha <- merge(alpha, covars, by.x=1, by.y=0)
rownames(alpha) <- alpha[,1]
alpha <- alpha[,-1]


# view alpha covariates to choose
for(y in 1:4)
{
	for(x in 5:ncol(alpha))
	{
		fname = paste0('1_guisane/img/alpha/', colnames(alpha)[y], '_', colnames(alpha)[x], '.pdf')
		pdf(fname)
		plot(alpha[,x], alpha[,y], xlab=colnames(alpha)[x], ylab=colnames(alpha)[y],
				 pch=20, col='blue')
		legend('topright', pch=20, col='blue', legend=paste(round(cor(alpha[,x], alpha[,y]),3)))
		dev.off()
	}
}

# choose alpha covariates
# chosen because they have the highest correlations with the response and correlations
# with each other < 0.7
richCovars <- c('bio6', 'bio15')
simpCovars <- c('bio3', 'bio13')
mpdFCovars <- c('bio6','bio15')
mpdPCovars <- c('bio5', 'bio15')

# choose 20% of sites as holdouts and write data to disk
fit_proportion <- 0.8
fit_rows <- sample(nrow(alpha), as.integer(fit_proportion*nrow(alpha)))
saveRDS(fit_rows, '1_guisane/dat/selected_rows_alpha.rds')
write.table(alpha[fit_rows,], '1_guisane/dat/alpha_fit.csv', sep=',', row.names=FALSE)
write.table(alpha[-fit_rows,], '1_guisane/dat/alpha_valid.csv', sep=',', row.names=FALSE)





emp_logit <- function(x)
{
	eps <- 0
	if(any(x == 0 | x == 1))
	{
		eps <- min(x[x != 0])
	}
	log((x+eps) / (1-x+eps))
}


betaCovars <- env_dissim(covars, sites=0, scale=FALSE)
# for some reason, our site names are switched
cols <- which(colnames(betaCovars) %in% c('site1', 'site2'))
colnames(betaCovars)[cols] <- rev(colnames(betaCovars)[cols])

# merge covariates; we only keep the covariate rows, because it is just the lower diagonal and the matrices are symmetric
beta_explo <- merge(beta, betaCovars, all.x=FALSE, all.y=TRUE)
beta_explo$sorensen_t <- emp_logit(beta_explo$sorensen_t)

# look over beta covariates
# for(y in 3:5)
# {
# 	yn <- colnames(beta_explo)[y]
# 	yy <- beta_explo[,y]
# 	# if(yn == 'bray_curtis_t') {
# 	# 	kp <- which(yy < 8)
# 	# 	yy <- yy[kp]
# 	# }
# 	for(x in 6:ncol(beta_explo))
# 	{
# 		xx <- beta_explo[,x]
# 		xn <- colnames(beta_explo)[x]
# 		# if(yn == 'bray_curtis_t')	xx <- xx[kp]
# 		fname = paste0('img/beta/', yn, '_', xn, '.pdf')
# 		pdf(fname)
# 		plot(xx, yy, xlab=colnames(beta_explo)[x], ylab=colnames(beta_explo)[y],
# 				 pch=20, col='blue')
# 		legend('bottomright', pch=20, col='blue', legend=paste(round(cor(xx, yy),3)))
# 		dev.off()
# 	}
# }

# using the same vars for all ys for now
bVars <- c('bio3', 'bio6', 'bio12', 'bio15')

# compute separate beta sets (because distance metric is different with different variables)
finalBetaCovars <- env_dissim(covars[,bVars], sites=0, scale=FALSE)
# for some reason, our site names are switched
cols <- which(colnames(finalBetaCovars) %in% c('site1', 'site2'))
colnames(finalBetaCovars)[cols] <- rev(colnames(finalBetaCovars)[cols])

beta_final <- merge(beta, finalBetaCovars, all.x=FALSE, all.y=TRUE)


# for beta, the fit takes longer, so we fit to 50%
fit_proportion <- 0.5
fit_rows <- sample(nrow(beta_final), as.integer(fit_proportion*nrow(beta_final)))
saveRDS(fit_rows, '1_guisane/dat/selected_rows_beta.rds')
write.table(beta_final[fit_rows,], '1_guisane/dat/beta_fit.csv', sep=',', row.names=FALSE)
write.table(beta_final[-fit_rows,], '1_guisane/dat/beta_valid.csv', sep=',', row.names=FALSE)


