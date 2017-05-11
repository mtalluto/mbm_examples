#!/usr/bin/env Rscript
library(reshape2)
library(fields)
library(RColorBrewer)

res <- '2_tasmania/res'
rc1d_tax <- read.csv('2_tasmania/res/taxo/respCurve1D.csv')
rc1d_fun <- read.csv('2_tasmania/res/func/respCurve1D.csv')
dat <- read.csv("2_tasmania/dat/tasFit_1.csv")
cols <- c("#b24b70", "#6d81b1", "#8f4bbf", "#6c9d53", "#bb6f3d")


quartzFonts(ssp=c('Source Sans Pro', 'Source Sans Pro Bold', 'Source Sans Pro Italic', 'Source Sans Pro Bold Italic'))

get.window <- function(x, y, lower, upper, w.size = diff(range(x)) / 25, length.out=1000)
{
	do.window <- function(x, nx, y, w.size) {
		sapply(nx, function(xx) {
			i <- which(x >= xx - w.size & x <= xx + w.size)
			mean(y[i])
		})
	}
	newx <- seq(min(x), max(x), length.out = length.out)
	data.frame(x = newx, y = do.window(x, newx, y, w.size), lower = do.window(x, newx, lower, w.size), upper = do.window(x, newx, upper, w.size))
}

quartz(file='2_tasmania/img/fig4.png', h=4.5, w=7, type='png', dpi=150)
par(mfrow=c(1,2), bty='n', mar=c(4,3,0,0), oma=c(0.5,0.5,0.5,0.5), mgp=c(2, 0.5, 0), tcl=-0.2, family='ssp')

plot(dat$distance, dat$taxSorensen, cex=0.5, pch='+', ylim=c(0,1), ylab='Sorensen Taxonomic Distance', xlab="Environmental Distance")
with(rc1d_tax[1:250,], {
	lines(distance, pnorm(mean), col=cols[3])
	polygon(c(distance, rev(distance)), c(pnorm(lower), rev(pnorm(upper))), col=paste0(cols[3], '88'), border=NA)
})

plot(dat$distance, dat$functionalDistance, cex=0.5, pch='+', ylim=c(0,0.4), ylab='Functional Distance', xlab="Environmental Distance")
datPr_f <- read.csv("2_tasmania/res/func/datPredict.csv")
prPreds <- get.window(datPr_f$distance, datPr_f$mean, datPr_f$lower, datPr_f$upper)
lines(prPreds$x, pnorm(prPreds$y), col=cols[3])
polygon(c(prPreds$x, rev(prPreds$x)), c(pnorm(prPreds$lower), rev(pnorm(prPreds$upper))), col=paste0(cols[3], '88'), border=NA)
# with(rc1d_fun[1:250,], {
# 	lines(distance, pnorm(mean), col=cols[5])
# 	polygon(c(distance, rev(distance)), c(pnorm(lower), rev(pnorm(upper))), col=paste0(cols[5], '88'), border=NA)
# })

dev.off()
