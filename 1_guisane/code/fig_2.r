#!/usr/bin/env Rscript
library(gam)

alphaFit <- read.csv('1_guisane/dat/alpha_fit.csv')
alphaValid <- read.csv('1_guisane/dat/alpha_valid.csv')
rcSimp <- read.csv(file.path('1_guisane', 'res', 'simpson', 'respCurve1D.csv'), header=TRUE)
rcFunc <- read.csv(file.path('1_guisane', 'res', 'funcAlpha', 'respCurve1D.csv'), header=TRUE)
rcPhy <- read.csv(file.path('1_guisane', 'res', 'phyloAlpha', 'respCurve1D.csv'), header=TRUE)
source('1_guisane/code/mbm_support.r')

cols <- c("#b24b70", "#6d81b1", "#8f4bbf", "#6c9d53", "#bb6f3d")
quartzFonts(ssp=c('Source Sans Pro', 'Source Sans Pro Bold', 'Source Sans Pro Italic', 'Source Sans Pro Bold Italic'))

# rescale vars
scaling <- readRDS("1_guisane/dat/covariateScaling.rds")
alphaFit <- descale(alphaFit, scaling)
alphaValid <- descale(alphaValid, scaling)
rcSimp <- descale(rcSimp, scaling)
rcFunc <- descale(rcFunc, scaling)
rcPhy <- descale(rcPhy, scaling)

# fig settings
lp=c(-10.5, 75)
l.cex=1.2
cex.l=0.7
lab.l = 2.5

quartz(file='1_guisane/img/fig2.pdf', h=9, w=6, type='pdf')
par(mfrow=c(3,2), bty='n', mar=c(3.4,3,0,0), oma=c(0.5,0.5,0.5,0.5), mgp=c(2, 0.5, 0), tcl=-0.2, family='ssp')
yl <- c(0,1)
mod <- gam(I(qnorm(simpson)) ~ s(bio3, 2) + s(bio13, 2), data=alphaFit)
with(rcSimp[1:250,], resp.curve(bio3, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio3', 'simpson')], xyValid = alphaValid[,c('bio3', 'simpson')], 
	link = pnorm, col=cols[1], ylim=yl, xlab='Isothermality', ylab=""))
mtext('Simpson Diversity', side=2, line=lab.l, cex=cex.l)
xx <- seq(38,42,.01)
yy <- pnorm(predict(mod, newdata=data.frame(bio3=xx, bio13=mean(alphaFit$bio13))))
lines(xx, yy, col=cols[2])

with(rcSimp[251:500,], resp.curve(bio13, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio13', 'simpson')], xyValid = alphaValid[,c('bio13', 'simpson')], 
	link = pnorm, col=cols[1], ylim=yl, xlab='Precipitation of the Wettest Month', ylab=""))
xx <- seq(100,150,.01)
yy <- pnorm(predict(mod, newdata=data.frame(bio13=xx, bio3=mean(alphaFit$bio3))))
lines(xx, yy, col=cols[2])



yl <- c(0,5)
mod <- gam(mpd_af ~ s(bio6, 2) + s(bio15, 2), data=alphaFit)
with(rcFunc[1:250,], resp.curve(bio6, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio6', 'mpd_af')], xyValid = alphaValid[,c('bio6', 'mpd_af')], 
	col=cols[1], ylim=yl, xlab='Min Temp of Coldest Month (°C)', ylab=""))
mtext('Functional MPD', side=2, line=lab.l, cex=cex.l)
xx <- seq(-13,-5,.01)
yy <- predict(mod, newdata=data.frame(bio6=xx, bio15=mean(alphaFit$bio15)))
lines(xx, yy, col=cols[2])

with(rcFunc[251:500,], resp.curve(bio15, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio15', 'mpd_af')], xyValid = alphaValid[,c('bio15', 'mpd_af')], 
	col=cols[1], ylim=yl, xlab='Precipitation Seasonality', ylab=""))
xx <- seq(15,22,.01)
yy <- predict(mod, newdata=data.frame(bio15=xx, bio6=mean(alphaFit$bio6)))
lines(xx, yy, col=cols[2])


yl <- c(0,20)
mod <- gam(mpd_ap ~ s(bio5, 2) + s(bio15, 2), data=alphaFit)
with(rcPhy[1:250,], resp.curve(bio5, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio5', 'mpd_ap')], xyValid = alphaValid[,c('bio5', 'mpd_ap')], 
	col=cols[1], ylim=yl, xlab='Max Temp of Warmest Month (°C)', ylab=""))
mtext('Phylogenetic MPD', side=2, line=lab.l, cex=cex.l)
xx <- seq(16,22,.01)
yy <- predict(mod, newdata=data.frame(bio5=xx, bio15=mean(alphaFit$bio15)))
lines(xx, yy, col=cols[2])

with(rcPhy[251:500,], resp.curve(bio15, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio15', 'mpd_ap')], xyValid = alphaValid[,c('bio15', 'mpd_ap')], 
	col=cols[1], ylim=yl, xlab='Precipitation Seasonality', ylab=""))
xx <- seq(15,22,.01)
yy <- predict(mod, newdata=data.frame(bio15=xx, bio5=mean(alphaFit$bio5)))
lines(xx, yy, col=cols[2])

dev.off()
