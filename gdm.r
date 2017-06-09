library("gdm")
library("mbm")

# load the model
taxMBM <- readRDS("3_alps/res/taxModelMF.rds")
phyMBM <- readRDS("3_alps/res/phyModel.rds")
taxBeta <- readRDS("3_alps/res/gdmdat/taxBeta.rds")
phyBeta <- readRDS("3_alps/res/gdmdat/phyBeta.rds")
coordDat <- readRDS("3_alps/res/gdmdat/coords.rds")
envDat <- readRDS("3_alps/res/gdmdat/envMat.rds")

bioDat <- cbind(as.numeric(rownames(taxBeta$fit)), taxBeta$fit)
colnames(bioDat)[1] <- "site"
envScale <- scale(envDat$fit)
env <- merge(coordDat$fit, envScale, by=0)
colnames(env)[1] <- 'site'

gdmDat <- formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env)
gdmMod <- gdm(gdmDat, geo=FALSE)

# plot response curves against distance
tp <- "pdf"
cols <- c("#ff644b", "#005ab8")

quartz(width=6.5, height=3.5, pointsize=10, type=tp, file=paste0("3_alps/img/gdm.", tp), bg='white', dpi=600)
par(mfrow=c(1,2), bty='n', mar=c(4,4,0,0), bty='n', oma=c(0, 0, 0.5, 0.5))

rc(taxMBM, cex_pt=0.2, xlab="Environmental Distance", ylab="Sørensen Dissimilarity", ylim=c(0,1), lwd=2, pch=16, col_line=cols[1])
# mtext("Sørensen Dissimilarity", side=2, line=3)
x <- gdmMod
plot(x$ecological, x$observed, xlab = "GDM Ecological Distance", 
	 ylab = "", type = "n", 
	 ylim = c(0, 1))
points(x$ecological, x$observed, pch = 16, cex = 0.2, col = "#666666")
overlayX <- seq(from = min(x$ecological), to = max(x$ecological), 
				length = 200)
overlayY <- 1 - exp(-overlayX)
lines(overlayX, overlayY, lwd = 2, col=cols[1])

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
rmse(taxBeta$valid[lower.tri(taxBeta$valid)], taxMBM$y_rev_transform(taxMBM$y_rev_transform(taxMBM$predictions$valid[,1])))



## do rmse for a PD model as well
# first need to scale from 0 to 1 for GDM
bioDat.p <- phyBeta$fit
bioDat.ps <- (bioDat.p - mean(bioDat.p))/sd(bioDat.p)
bioDat.pst <- pnorm(bioDat.ps)
# bioDat.pst <- (bioDat.p - min(bioDat.p))/max(bioDat.p)
bioDat.pstr <- cbind(as.numeric(rownames(bioDat.pst)), bioDat.pst)

colnames(bioDat.pstr)[1] <- "site"

gdmDat.p <- formatsitepair(bioDat.pstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env)
gdmMod.p <- gdm(gdmDat.p, geo=FALSE)

resp <- qnorm(gdmDat.p$distance) * sd(bioDat.p) + mean(bioDat.p)
pred <- qnorm(predict(gdmMod.p, gdmDat.p)) * sd(bioDat.p) + mean(bioDat.p)
rmse(resp, pred)
rmse(phyMBM$y_rev_transform(phyMBM$response), phyMBM$y_rev_transform(phyMBM$fitted.values))

bioDat.pv <-  phyBeta$valid
bioDat.pvs <- (bioDat.pv - mean(bioDat.p))/sd(bioDat.p)
bioDat.pvst <- pnorm(bioDat.pvs)
bioDat.pvstr <- cbind(as.numeric(rownames(bioDat.pvst)), bioDat.pvst)
colnames(bioDat.pvstr)[1] <- "site"
gdmValidDat.p <- formatsitepair(bioDat.pvstr, 3, siteColumn='site', XColumn = 'x', YColumn = 'y', predData = env.v)

resp.v <- qnorm(gdmValidDat.p$distance) * sd(bioDat.p) + mean(bioDat.p)
pred.v <- qnorm(predict(gdmMod.p, gdmValidDat.p)) * sd(bioDat.p) + mean(bioDat.p)
rmse(resp.v, pred.v)
rmse(phyBeta$valid[lower.tri(phyBeta$valid)], phyMBM$y_rev_transform(phyMBM$y_rev_transform(phyMBM$predictions$valid[,1])))
