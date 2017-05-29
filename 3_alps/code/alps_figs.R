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