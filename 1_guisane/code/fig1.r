#!/usr/bin/env Rscript
res <- '1_guisane/res'
source('1_guisane/code/mbm_support.r')
alphaFit <- read.csv('1_guisane/dat/alpha_fit.csv')
alphaValid <- read.csv('1_guisane/dat/alpha_valid.csv')

rcRichness <- read.csv(file.path(res, 'richness', 'respCurve1D.csv'), header=TRUE)
rcRichness1 <- read.csv(file.path(res, 'richness_1', 'respCurve1D.csv'), header=TRUE)
rcRichness0.4 <- read.csv(file.path(res, 'richness_0.4', 'respCurve1D.csv'), header=TRUE)

cols <- c("#b24b70", "#6d81b1", "#8f4bbf", "#6c9d53", "#bb6f3d")

# rescale bio 6
scaling <- readRDS("1_guisane/dat/covariateScaling.rds")
rcRichness$bio6 <- (rcRichness$bio6 * scaling['bio6', 'scale']) + scaling['bio6', 'center']
rcRichness1$bio6 <- (rcRichness1$bio6 * scaling['bio6', 'scale']) + scaling['bio6', 'center']
rcRichness0.4$bio6 <- (rcRichness0.4$bio6 * scaling['bio6', 'scale']) + scaling['bio6', 'center']
alphaFit$bio6 <- (alphaFit$bio6 * scaling['bio6', 'scale']) + scaling['bio6', 'center']
alphaValid$bio6 <- (alphaValid$bio6 * scaling['bio6', 'scale']) + scaling['bio6', 'center']


quartzFonts(ssp=c('Source Sans Pro', 'Source Sans Pro Bold', 'Source Sans Pro Italic', 'Source Sans Pro Bold Italic'))
yl <- c(0,80)
lp=c(-10.5, 75)
l.cex=1.2
cex.l=0.7
lab.l = 2.5

quartz(file='1_guisane/img/fig1.pdf', h=7, w=4, type='pdf')
par(mfrow=c(3,1), bty='n', mar=c(2,3,0,0), oma=c(2,0.5,0.5,0.5), mgp=c(2, 0.5, 0), tcl=-0.2, family='ssp')
lensc <- round(read.csv('1_guisane/res/richness/params.csv', h=F)[2,1], 2)
with(rcRichness[1:250,], resp.curve(bio6, mean, lower=lower, upper=upper,
			xyDat=alphaFit[,c('bio6', 'richness')], xyValid = alphaValid[,c('bio6', 'richness')], 
			link = exp, col=cols[1], ylim=yl, xlab='', ylab=""))
text(lp[1], lp[2], paste("l =", lensc), cex=l.cex)
with(rcRichness1[1:250,], resp.curve(bio6, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio6', 'richness')], xyValid = alphaValid[,c('bio6', 'richness')], 
	link = exp, col=cols[1], ylim=yl, xlab='', ylab=""))
text(lp[1], lp[2], "l = 1.00", cex=l.cex)
mtext('Species Richness', side=2, line=lab.l, cex=cex.l)
with(rcRichness0.4[1:250,], resp.curve(bio6, mean, lower=lower, upper=upper,
	xyDat=alphaFit[,c('bio6', 'richness')], xyValid = alphaValid[,c('bio6', 'richness')], 
	link = exp, col=cols[1], ylim=yl, xlab='', ylab=""))
text(lp[1], lp[2], "l = 0.40", cex=l.cex)
mtext('Minnimum temperature of coldest month (°C)', side=1, line=lab.l, cex=cex.l)
dev.off()


