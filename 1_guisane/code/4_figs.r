#!/usr/bin/env Rscript
res <- '1_guisane/res'

rcRichness <- read.csv(file.path(res, 'richness', 'respCurve1D.csv'), header=TRUE)
rcRP <- read.csv(file.path(res, 'richnessPois', 'respCurve1D.csv'), header=TRUE)
rcRichnessS <- read.csv(file.path(res, 'richnessSamples', 'respCurve1D.csv'), header=TRUE)
rcRPS <- read.csv(file.path(res, 'richnessPoisSamples', 'respCurve1D.csv'), header=TRUE)
alphaFit <- read.csv('1_guisane/dat/alpha_fit.csv')
cols <- c("#b24b70", "#6c9d53", "#8f4bbf", "#bb6f3d", "#6d81b1")

rcRichness$lower = rcRichness$mean - 1.96 * rcRichness$sd
rcRichness$upper = rcRichness$mean + 1.96 * rcRichness$sd

resp.curve <- function(x, y, lower, upper, x0, y0, link, col='#6d81b1', alpha="66", ...)
{
	if(missing(link)) link <- function(x) x
	plot(x, link(y), col=col, type='l', ...)
	polygon(c(x, rev(x)), c(link(lower), rev(link(upper))), border=NA, col=paste0(col, alpha))
	points(x0, y0, pch='x', cex=0.7)
}

## RICHNESS is a bit f'ed
with(rcRichness[1:250,], {
	plot(bio6, exp(mean), col=cols[1], type='l', ylim=c(0, 100))
	polygon(c(bio6, rev(bio6)), c(exp(lower), rev(exp(upper))), border=NA, col=paste0(cols[1], '66'))
})
with(rcRP[1:250,], {
	lines(bio6, exp(mean), col=cols[2])
	polygon(c(bio6, rev(bio6)), c(exp(lower), rev(exp(upper))), border=NA, col=paste0(cols[2], '66'))
})
with(rcRichnessS[1:250,], {
	lines(bio6, exp(mean), col=cols[3])
	polygon(c(bio6, rev(bio6)), c(exp(lower), rev(exp(upper))), border=NA, col=paste0(cols[3], '66'))
})
with(rcRPS[1:250,], {
	lines(bio6, exp(mean), col=cols[4])
	polygon(c(bio6, rev(bio6)), c(exp(lower), rev(exp(upper))), border=NA, col=paste0(cols[4], '66'))
})
points(alphaFit$bio6, alphaFit$richness, pch='x', cex=0.7)



# rcSimpson <- read.csv(file.path(res, 'alphaSimpson', 'respCurve.csv'), header=TRUE)
# with(rcSimpson[1:250,], resp.curve(bio3, mean, lower, upper, alphaFit$bio3, alphaFit$simpson, 
# 	plogis, ylim=c(0,1), xlab='bio3', ylab='Simpson diversity'))
# with(rcSimpson[251:500,], resp.curve(bio15, mean, lower, upper, alphaFit$bio15, alphaFit$simpson, 
# 	plogis, ylim=c(0,1), xlab='bio15', ylab='Simpson diversity'))
# 
# rcFunc <- read.csv(file.path(res, 'alphaFunc', 'respCurve.csv'), header=TRUE)
# with(rcFunc[1:250,], resp.curve(bio6, mean, lower, upper, alphaFit$bio6, alphaFit$mpd_af, 
# 			ylim=c(0,6), xlab='bio6', ylab='Functional MPD (within)'))
# with(rcFunc[251:500,], resp.curve(bio15, mean, lower, upper, alphaFit$bio15, alphaFit$mpd_af, 
# 			ylim=c(0,6), xlab='bio15', ylab='Functional MPD (within)'))
# 
# 
# rcPhy <- read.csv(file.path(res, 'alphaPhylo', 'respCurve.csv'), header=TRUE)
# with(rcPhy[1:250,], resp.curve(bio5, mean, lower, upper, alphaFit$bio5, alphaFit$mpd_ap, 
# 	ylim=c(0,350), xlab='bio5', ylab='Phylogenetic MPD (within)'))
# with(rcPhy[251:500,], resp.curve(bio15, mean, lower, upper, alphaFit$bio15, alphaFit$mpd_ap, 
# 	ylim=c(0,350), xlab='bio15', ylab='Phylogenetic MPD (within)'))
# 

