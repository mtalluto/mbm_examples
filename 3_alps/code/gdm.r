library("gdm")
library("mbm")
library("mbmdata")


# I checked this; there is no issue with the fit data not having reasonable values for sorensen
# even using all the data, we get more or less the same distribution
# data(alps)
# keep <- (alps$siteEnv[,'num_releves'] >= 5 | alps$siteEnv[,'elev'] > 3000)
# dat <- alps$siteSpecies[keep,]
# bioDat <- sorensen(dat)

# load the model
taxMBM <- readRDS("3_alps/res/taxModelMF.rds")
taxMBMnoMF <- readRDS("3_alps/res/taxModel.rds")
phyMBM <- readRDS("3_alps/res/phyModel.rds")
funMBM <- readRDS("3_alps/res/funModel.rds")
taxBeta <- readRDS("3_alps/res/gdmdat/taxBeta.rds")
phyBeta <- readRDS("3_alps/res/gdmdat/phyBeta.rds")
funBeta <- readRDS("3_alps/res/gdmdat/funBeta.rds")
coordDat <- readRDS("3_alps/res/gdmdat/coords.rds")
envDat <- readRDS("3_alps/res/gdmdat/envMat.rds")

bioDat <- cbind(as.numeric(rownames(taxBeta$fit)), taxBeta$fit)
colnames(bioDat)[1] <- "site"
envScale <- scale(envDat$fit)
env <- merge(coordDat$fit, envScale, by=0)
env[,1] <- as.integer(env[,1])
colnames(env)[1] <- 'site'


# try a variety of GDM models; first use only single covariates
cos <- c('site', 'x', 'y')
# gdm1Dat <- list(
# 	bio4 = formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_4')]),
# 	bio6 = formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_6')]),
# 	bio7 = formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_7')]),
# 	bio15 = formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_15')])
# 	)
# gdm1Models <- lapply(gdm1Dat, gdm, geo=FALSE)
# gdm1ModelsSp <- lapply(gdm1Dat, gdm, geo=TRUE)

## bio 6 is best

# best TD GDM - bio_7 and bio_15
gdmDat <- formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_6')])
gdmMod <- gdm(gdmDat, geo=FALSE)


# this was for testing all variables
# env2 <- cbind(site=as.integer(rownames(alps$siteEnv[rows$fit,])), alps$siteEnv[rows$fit,-4])
# env2[,4:ncol(env2)] <- scale(env2[,4:ncol(env2)])
env
# i <- 23  # i from 4 to 23
# gdmTestDat <- formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env2[,c(1:3, i)])
# gdmTestMod <- gdm(gdmTestDat, geo=FALSE)
# summary(gdmTestMod)
# plot(gdmTestMod)


# FIGURE 4: plot response curves against distance
tp <- "tiff"
cols <- c("#ff644b", "#005ab8")

quartz(width=6.2, height=8.5, pointsize=10, type=tp, file=paste0("3_alps/img/gdm.", tp), bg='white', dpi=300)
par(mfrow=c(3,2), bty='n', mar=c(4,4,0,0), bty='n', oma=c(0, 0, 0.5, 0.5))

rc(taxMBM, cex_pt=0.2, xlab="Climatic Distance", ylab="Sørensen Dissimilarity", ylim=c(0,1.1), lwd=2, pch=16, col_line=cols[1])
rc(taxMBMnoMF, add=TRUE, col_pt = NA, col_line = cols[2], lwd=2)

text(grconvertX(0.02, from='npc'), grconvertY(0.92, from='npc'), 'A', pos=4, cex=1.5)
# mtext("Sørensen Dissimilarity", side=2, line=3)
x <- gdmMod
plot(x$ecological, x$observed, xlab = "GDM Ecological Distance", 
	 ylab = "", type = "n", 
	 ylim = c(0, 1.1))
points(x$ecological, x$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(x$ecological), to = max(x$ecological), 
				length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[1])
text(grconvertX(0.02, from='npc'), grconvertY(0.92, from='npc'), 'B', pos=4, cex=1.5)
# dev.off()



#### fig S2 - GDM for PD and FD
# first need to scale from 0 to 1 for GDM
bioDat.p <- phyBeta$fit
bioDat.ps <- (bioDat.p - mean(bioDat.p))/sd(bioDat.p)
bioDat.pst <- pnorm(bioDat.ps)
# bioDat.pst <- (bioDat.p - min(bioDat.p))/max(bioDat.p)
bioDat.pstr <- cbind(as.numeric(rownames(bioDat.pst)), bioDat.pst)
colnames(bioDat.pstr)[1] <- "site"

# gdmPTest <- lapply(4:23, function(i) formatsitepair(bioDat.pstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env2[,c(1:3, i)]))
# gdmPMod <- lapply(gdmPTest, gdm, geo=FALSE)

# the only one that converged
gdmDat.p <- formatsitepair(bioDat.pstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_15')])
gdmMod.p <- gdm(gdmDat.p, geo=FALSE)


bioDat.f <- funBeta$fit
# bioDat.fs <- (bioDat.f - mean(bioDat.f))/sd(bioDat.f)
# bioDat.fst <- pnorm(bioDat.fs)
# bioDat.fstr <- cbind(as.numeric(rownames(bioDat.fst)), bioDat.fst)
# colnames(bioDat.fstr)[1] <- "site"

# for testing -- no response scaling -- I like this better
bioDat.fstr <- cbind(as.numeric(rownames(bioDat.f)), bioDat.f)
colnames(bioDat.fstr)[1] <- "site"

# gdmFTest <- lapply(4:23, function(i) formatsitepair(bioDat.fstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env2[,c(1:3, i)]))
# gdmFMod <- lapply(gdmFTest, gdm, geo=FALSE)

# these 2 are pretty good
gdmDat.f <- formatsitepair(bioDat.fstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env[,c(cos, 'bio_4', 'bio_6')])
gdmMod.f <- gdm(gdmDat.f, geo=FALSE)
# summary(gdmMod.f)



# quartz(width=6.5, height=6.5, pointsize=10, type=tp, file=paste0("3_alps/img/gdms2.", tp), bg='white', dpi=600)
# par(mfrow=c(2,2), bty='n', mar=c(4,4,0,0), bty='n', oma=c(0, 0, 0.5, 0.5))

rc(funMBM, cex_pt=0.2, xlab="Climatic Distance", ylab="Functional MPD", lwd=2, pch=16, col_line=cols[1])
# mtext("Sørensen Dissimilarity", side=2, line=3)
text(grconvertX(0, from='npc'), grconvertY(0.85, from='npc'), 'C', pos=4, cex=1.5)
x <- gdmMod.f
plot(x$ecological, x$observed, xlab = "GDM Ecological Distance", 
	 ylab = "", type = "n"
	 )
points(x$ecological, x$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(x$ecological), to = max(x$ecological), 
				length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[1])
text(grconvertX(0, from='npc'), grconvertY(0.85, from='npc'), 'D', pos=4, cex=1.5)



rc(phyMBM, cex_pt=0.2, xlab="Climatic Distance", ylab="Phylogenetic MPD", lwd=2, pch=16, col_line=cols[1])
text(grconvertX(0.05, from='npc'), grconvertY(0.9, from='npc'), 'E', pos=4, cex=1.5)
# mtext("Sørensen Dissimilarity", side=2, line=3)
x <- gdmMod.p
plot(x$ecological, x$observed, xlab = "GDM Ecological Distance", 
	 ylab = "", type = "n", 
	 ylim = c(0, 1.1))
points(x$ecological, x$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(x$ecological), to = max(x$ecological), 
				length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[1])
text(grconvertX(0.05, from='npc'), grconvertY(0.9, from='npc'), 'F', pos=4, cex=1.5)
dev.off()





# compute rmse for both using both fit and validation data
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

## fit data
rmse(gdmDat$distance, predict(gdmMod, gdmDat))
rmse(taxMBM$y_rev_transform(taxMBM$response), taxMBM$y_rev_transform(taxMBM$fitted.values))

bioDat.v <- cbind(as.numeric(rownames(taxBeta$valid)), taxBeta$valid)
colnames(bioDat.v)[1] <- "site"
envScale.v <- scale(envDat$valid, center=attr(envScale, "scaled:center"), scale=attr(envScale, "scaled:scale"))
env.v <- merge(coordDat$valid, envScale.v, by=0)
colnames(env.v)[1] <- 'site'
gdmValidDat <- formatsitepair(bioDat.v, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env.v)

rmse(gdmValidDat$distance, predict(gdmMod, gdmValidDat))
rmse(taxBeta$valid[lower.tri(taxBeta$valid)], taxMBM$y_rev_transform(taxMBM$predictions$valid[,3]))



## do rmse for FD and PD models as well

resp <- qnorm(gdmDat.p$distance) * sd(bioDat.p) + mean(bioDat.p)
pred <- qnorm(predict(gdmMod.p, gdmDat.p)) * sd(bioDat.p) + mean(bioDat.p)
rmse(resp, pred)
rmse(phyMBM$response, phyMBM$fitted.values)

bioDat.pv <-  phyBeta$valid
bioDat.pvs <- (bioDat.pv - mean(bioDat.p))/sd(bioDat.p)
bioDat.pvst <- pnorm(bioDat.pvs)
bioDat.pvstr <- cbind(as.numeric(rownames(bioDat.pvst)), bioDat.pvst)
colnames(bioDat.pvstr)[1] <- "site"
gdmValidDat.p <- formatsitepair(bioDat.pvstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env.v)

resp.v <- qnorm(gdmValidDat.p$distance) * sd(bioDat.p) + mean(bioDat.p)
pred.v <- qnorm(predict(gdmMod.p, gdmValidDat.p)) * sd(bioDat.p) + mean(bioDat.p)
rmse(resp.v, pred.v)
rmse(phyBeta$valid[lower.tri(phyBeta$valid)], phyMBM$predictions$valid[,3])




resp <- qnorm(gdmDat.f$distance) * sd(bioDat.f) + mean(bioDat.f)
pred <- qnorm(predict(gdmMod.f, gdmDat.f)) * sd(bioDat.f) + mean(bioDat.f)
rmse(resp, pred)
rmse(funMBM$response, funMBM$fitted.values)

bioDat.fv <-  funBeta$valid
bioDat.fvs <- (bioDat.fv - mean(bioDat.f))/sd(bioDat.f)
bioDat.fvst <- pnorm(bioDat.fvs)
bioDat.fvstr <- cbind(as.numeric(rownames(bioDat.fvst)), bioDat.fvst)
colnames(bioDat.fvstr)[1] <- "site"
gdmValidDat.f <- formatsitepair(bioDat.fvstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env.v)

resp.v <- qnorm(gdmValidDat.f$distance) * sd(bioDat.f) + mean(bioDat.f)
pred.v <- qnorm(predict(gdmMod.f, gdmValidDat.f)) * sd(bioDat.f) + mean(bioDat.f)
rmse(resp.v, pred.v)
rmse(funBeta$valid[lower.tri(funBeta$valid)], funMBM$predictions$valid[,3])
