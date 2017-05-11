#!/usr/bin/env Rscript

# tbEps <- readRDS('dat/betaUntransformed.rds')$sorensen_t
# tbEps <- min(tbEps[tbEps != 0])
# tb <- read.csv('res/taxBeta/respCurve1D.csv')
# emp_unlogit <- function(x, eps)
# 	((exp(x) * (1 + eps)) - eps)/(1 + exp(x))
library(reshape2)
library(RColorBrewer)
library(fields)

cols <- c("#b24b70", "#6d81b1", "#8f4bbf", "#6c9d53", "#bb6f3d")
quartzFonts(ssp=c('Open Sans', 'Open Sans Bold', 'Open Sans Italic', 'Open Sans Bold Italic'))
# jet.colors <- colorRampPalette( c("#110000", "#b04b48", "#ff7777") ) 
# nbcol <- 100
# color <- jet.colors(nbcol)
# 
# makecols <- function(z, cl) {
# 	nrz=nrow(z); ncz=ncol(z)
# 	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# 	# Recode facet z-values into color indices
# 	facetcol <- cut(zfacet, nbcol)
# 	cl[facetcol]
# }

quartz(file='img/fig3.pdf', h=9, w=6, type='pdf')
par(mfcol=c(3,2), bty='n', mar=c(3.4,3,0,0), oma=c(0.5,0.5,0.5,0.5), mgp=c(2, 0.5, 0), tcl=-0.2, family='ssp')
tb <- read.csv('res/taxBeta/respCurve1D.csv')
tbFit <- read.csv("dat/betaTaxo_fit.csv")
yl <- c(0,1)
with(tbFit, plot(distance, sorensen_t, cex=0.2, pch=20, ylim=yl, ylab=bquote(beta~Taxonomic), xlab="Env. Distance"))
with(tb[1:250,], {
	lines(distance, plogis(mean), col=cols[1], lwd=1.5)
	polygon(c(distance, rev(distance)), c(plogis(lower), rev(plogis(upper))), border=NA, col=paste0(cols[1], '88'))
})


fb <- read.csv('res/funBeta/respCurve1D.csv')
fbFit <- read.csv("dat/betaFun_fit.csv")
yl <- c(0,6)
with(fbFit, plot(distance, mpd_bf, cex=0.2, pch=20, ylim=yl, ylab=bquote(beta~Functional~MPD), xlab="Env. Distance"))
with(fb[1:250,], {
	lines(distance, mean, col=cols[1], lwd=1.5)
	polygon(c(distance, rev(distance)), c(lower, rev(upper)), border=NA, col=paste0(cols[1], '88'))
})

pb <- read.csv('res/phyBeta/respCurve1D.csv')
pbFit <- read.csv("dat/betaPhylo_fit.csv")
yl <- c(12,18)
with(pbFit, plot(distance, mpd_bp, cex=0.2, pch=20, ylim=yl, ylab=bquote(beta~Phylogenetic~MPD), xlab="Env. Distance"))
with(pb[1:250,], {
	lines(distance, mean, col=cols[1], lwd=1.5)
	polygon(c(distance, rev(distance)), c(lower, rev(upper)), border=NA, col=paste0(cols[1], '88'))
})

par(mar=c(3.4,3,3,1))
color = colorRampPalette(brewer.pal(7, "OrRd"))(100)

tb2 <- read.csv('res/taxBeta/responseCurve2D.csv')
z = plogis(acast(tb2, distance ~ bio6, value.var='mean'))
x = as.numeric(colnames(z))
y = as.numeric(rownames(z))
x = c(x, x[length(x)] + x[2] - x[1])
y = c(y, y[length(y)] + y[2] - y[1])
# persp(x=as.numeric(colnames(z)), y=as.numeric(rownames(z)), z=z, xlab='\nMin Temperature\nColdest Month', ylab="Env. Distance",
# 	zlab="\nβ Taxonomic", theta=60, phi=30, lwd=0.5, col=makecols(z, color), ticktype='detailed')
image.plot(x=x, y=y, z=z, col=color, ylab="Env Distance", xlab="\nMin Temperature\nColdest Month", main="β Taxonomic")
contour(x=x[1:100], y=y[1:100], z=z, add=TRUE, nlevels=5)

plot.new()
fb2 <- read.csv('res/funBeta/responseCurve2D.csv')
z = acast(fb2, distance ~ bio6, value.var='mean')
persp(x=as.numeric(colnames(z)), y=as.numeric(rownames(z)), z=z, xlab='\nMin Temperature\nColdest Month', ylab="Env. Distance",
	zlab="\nβ Functional MPD", theta=60, phi=30, border=NA, lwd=0.2, col=makecols(z, color), ticktype='detailed')

plot.new()
pb2 <- read.csv('res/phyBeta/responseCurve2D.csv')
z = acast(pb2, distance ~ bio15, value.var='mean')
x = as.numeric(colnames(z))
y = as.numeric(rownames(z))
x = c(x, x[length(x)] + x[2] - x[1])
y = c(y, y[length(y)] + y[2] - y[1])
image.plot(x=x, y=y, z=z, col=color, ylab="Env Distance", xlab="Precip Seasonality", main="β Phylogenetic")
# persp(x=as.numeric(colnames(z)), y=as.numeric(rownames(z)), z=z, xlab='\nPrecip Seasonality', ylab="Env. Distance",
# 	zlab="\nβ Phylogenetic MPD", theta=60, phi=30, lwd=0.5, col=makecols(z, color), ticktype='detailed')

dev.off()




