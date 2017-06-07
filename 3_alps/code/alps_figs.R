## FIG 1 - look at response curves of 3 models
# plot(taxModel, sterr=TRUE)
# par(mfrow = c(3, 1))
# rc(taxModel, ylim=c(0,1), ylab="Sørensen Dissimilarity", xlab="Environmental Distance")
# rc(funModel, ylab = "Functional MPD", xlab = "Environmental Distance")
# rc(phyModel, ylab = "Phylogenetic MPD", xlab = "Environmental Distance")




####
# fig 3: compare response curves
rc(taxModel, col_line='#c6000b', col_polygon=NA, lwd=2, xlab="Environmental Distance", ylab="Sørensen Dissimilarity")
rc(taxModel0.5, add=TRUE, col_line='#7851f7', col_polygon=NA, lwd=2)
rc(taxModel2, add=TRUE, col_line='#00f6f9', col_polygon=NA, lwd=2)
legend('bottomright', c(paste('l =', round(taxModel$params[2], 1)), paste('l =', round(taxModel2$params[2], 1)), 
						paste('l =', round(taxModel0.5$params[2], 1))), lwd=1.5, col=c('#c6000b', '#00f6f9', '#7851f7'))



mods <- list.dirs("3_alps/res", full.names=FALSE, recursive=FALSE)
betaDiv <- read.csv("3_alps/dat/betaDiv.csv")
betaDiv.v <- read.csv("3_alps/dat/betaDiv_valid.csv")

ws <- 0.5
wstep <- 100

qt <- 'png'
quartz(w=12, h=20, type=qt, file= paste0('3_alps/img/rc.', qt), dpi=300, bg='white')
par(mfrow=c(5,3))
for(m in mods) {
	dpr <- read.table(file.path('3_alps', 'res', m, 'dat_predict.csv'), header=TRUE)
	# vpr <- read.table(file.path('3_alps', 'res', m, 'valid_predict.csv'))
	
	wmids <- seq(min(dpr$distance) + ws, max(dpr$distance) - ws, length.out=wstep)
	rc <- t(sapply(wmids, function(x) {
		with(dpr[dpr$distance >= x - ws & dpr$distance < x + ws,], c(x, mean(mean)))
	}))
	fun <- pnorm
	if(grepl('mpd', m)) fun <- function(x) 10^x
	plot(rc[,1], fun(rc[,2]), type='n', ylim=range(log(betaDiv[,m])), xlab='distance', ylab=m)
	points(betaDiv[,'distance'], log(betaDiv[,m]), pch='x', cex=0.4)
	lines(rc[,1], fun(rc[,2]), col='cyan', lwd=2)
}
dev.off()



# rcDat <- read.csv("3_alps/res/sor/predict0.csv")
rcDat <- read.table("3_alps/res/sor/predict0.csv", h=TRUE)
## transform predictions
rcDat$mean <- pnorm(rcDat$mean)
rcDat$lower <- pnorm(rcDat$lower)
rcDat$upper <- pnorm(rcDat$upper)

betaDat <- readRDS('3_alps/dat/betaDiv.rds')
rc(rcDat$distance, rcDat$mean, betaDat$distance, betaDat$sor, rcDat[,c('lower', 'upper')], xlab='Environmental Distance', ylab='Taxonomic Dissimilarity')
