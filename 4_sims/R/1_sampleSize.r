# RNG seed, selected at random
source("R/landscape.r")
library("mbm")
library("gdm")

args <- commandArgs(TRUE)
do.gdm <- ("--gdm" %in% args)
do.mbm <- ("--mbm" %in% args)

set.seed(419124096)

sample_grid <- function(n, lscape) 
{
	require(raster)
	xrange <- range(lscape$landscape[,1])
	yrange <- range(lscape$landscape[,2])
	xstep <- diff(xrange)/sqrt(n)
	ystep <- diff(yrange)/sqrt(n)
	xstart <- runif(1, xrange[1], xrange[1] + xstep)
	ystart <- runif(1, yrange[1], yrange[1] + ystep)
	grid <- expand.grid(env1=xstart + (xstep * (seq_len(sqrt(n))-1)),
				env2 = ystart + (ystep * (seq_len(sqrt(n))-1)))

	lsRas <- rasterFromXYZ(cbind(lscape$landscape, 1:nrow(lscape$landscape)))
	rows <- extract(lsRas, grid)
	return(list(X=lscape$landscape[rows,], Y=lscape$site_sp[rows,], rows=rows))
}


# sample sizes (number of sites) to test
# svgp parameters for runs that need it
simParams <- data.frame(
	sampSizes  = c(25, 49, 100, 225, 529, 1024, 2025),
	batchSizes = c(NA, NA, 40,  60,  80,  120,  120 ),
	indSizes   = c(NA, NA, 40,  60,  80,  120,  200 ),
	iterLen    = c(NA, NA, 5e3, 5e3, 1e4, 1e4,  1e5 ))

# how many sites for validation
validN <- 200

# make a landscape and stuff here
lsFname <- 'dat/landscape.rds'
validFname <- 'dat/vDat.rds'
if(file.exists(lsFname)) {
	landscape <- readRDS(lsFname)
} else {
	landscape <- simulate_landscape(300, dims=c(100,100))
	saveRDS(landscape, lsFname)
}
if(file.exists(validFname)) {
	grid_v <- readRDS(validFname)
} else {
	grid_v <- sample_grid(validN, landscape)
	saveRDS(grid_v, validFname)
}

dissim_v <- sorensen(grid_v$Y)
dissim_v_tall <- dissim_v
diag(dissim_v_tall) <- NA
dissim_v_tall[upper.tri(dissim_v_tall)] <- NA
dissim_v_tall <- melt(dissim_v_tall, na.rm=TRUE)
colnames(dissim_v_tall) <- c('site2', 'site1', 'obs')

# plot(landscape)
# quartz()
# density(landscape)



format_gdm <- function(X, Y)
{
	# set up data to match GDM expectations
	bioDat <- cbind(as.integer(sub("site([0-9]+)", "\\1", rownames(Y))), Y)
	colnames(bioDat)[1] <- 'site'
	X$site <- as.integer(sub("site([0-9]+)", "\\1", rownames(X)))
	X$x <- X$env1
	X$y <- X$env2
	formatsitepair(bioDat, 3, siteColumn='site', XColumn = 'x', 
		YColumn = 'y', predData = X)
}
run_gdm <- function(X, Y)
{
	gdmDat <- format_gdm(X, Y)
	gdm(gdmDat, geo=FALSE)
}
gdm_predict <- function(mod, newX, newY)
{
	gdmDat <- format_gdm(newX, newY)
	data.frame(observations = gdmDat$distance, fits =predict(mod, gdmDat))
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

statsRow <- function(df, mt, stat, n, val)
{
	if(missing(df)) {
		newdf <- data.frame(modelType = character(), statistic = character(), 
			n = numeric(), value = numeric(), stringsAsFactors = FALSE)
		return(newdf)
	} else {
		newdf <- data.frame(modelType = mt, statistic = stat, 
			n = n, value = val, stringsAsFactors = FALSE)
		return(rbind(df, newdf))
	}
}

stFile <- 'res/sim_stats.rds'
if(file.exists(stFile)) {
	stats <- readRDS(stFile)
} else {
	stats <- statsRow()
}

# eventually loop over sampSizes
for(i in 1:nrow(simParams))
{
	ss <- simParams$sampSizes[i]
	batchN <- simParams$batchSizes[i]
	indN <- simParams$indSizes[i]
	iterN <- simParams$iterLen[i]

	# set up data
	grFname <- paste0('dat/sim/sample_grid_', ss, '.rds')
	if(file.exists(grFname)) {
		grid <- readRDS(grFname)
	} else {
		grid <- sample_grid(ss, landscape)
		saveRDS(grid, grFname)
	}
	dissim <- sorensen(grid$Y)

	## GDM block
	if(do.gdm)
	{
		cat(paste0(as.character(Sys.time()), "     starting gdm n = ", ss, "\n"))
		gdmTime <- system.time(gdmMod <- run_gdm(grid$X, dissim))
		gdmTime <- unname(as.numeric(gdmTime[1] + gdmTime[3]))
		stats <- statsRow(stats, 'gdm', 'fit.time', ss, gdmTime)
		gdmPrTime <- system.time(gdmPr <- gdm_predict(gdmMod, grid_v$X, dissim_v))
		gdmPrTime <- unname(as.numeric(gdmPrTime[1] + gdmPrTime[3]))
		stats <- statsRow(stats, 'gdm', 'predict.time', ss, gdmPrTime)
		gdmValid <- rmse(gdmPr[[1]], gdmPr[[2]])
		stats <- statsRow(stats, 'gdm', 'rmse', ss, gdmValid)	
		saveRDS(gdmMod, paste0('res/gdmmod_', ss, '.rds'))	
		saveRDS(gdmPr, paste0('res/gdmpr_', ss, '.rds'))	
	}

	## MBM block
	if(do.mbm) 
	{
		if(ss <= 100)
		{
			cat(paste0(as.character(Sys.time()), "     starting mbm n = ", ss, "\n"))
			mbmTime <- system.time(mbmMod <- mbm(y=dissim, x=grid$X, link='probit', 
					scale=FALSE, force_increasing = TRUE))
			mbmTime <- unname(as.numeric(mbmTime[1] + mbmTime[3]))
			stats <- statsRow(stats, 'mbm', 'fit.time', ss, mbmTime)
			mbmPrTime <- system.time(mbmPr <- predict(mbmMod, grid_v$X))
			mbmPrTime <- unname(as.numeric(mbmPrTime[1] + mbmPrTime[3]))
			stats <- statsRow(stats, 'mbm', 'predict.time', ss, mbmPrTime)
			mbmPr <- merge(mbmPr, dissim_v_tall)
			mbmValid <- rmse(mbmPr$obs, mbmMod$y_rev_transform(mbmPr$fit))
			stats <- statsRow(stats, 'mbm', 'rmse', ss, mbmValid)

			## here - calculate conf interval accuracy

			saveRDS(mbmMod, paste0('res/mbmmod_', ss, '.rds'))	
			saveRDS(mbmPr, paste0('res/mbmpr_', ss, '.rds'))	

		} else {
			mbmTime <- mbmMod <- NA
		}

		if(ss >= 100)
		{
			cat(paste0(as.character(Sys.time()), "     starting svgp n = ", ss, "\n"))
			svgpTime <- system.time(svgpMod <- mbm(y = dissim, x = grid$X, link = 'probit', 
				scale = FALSE, force_increasing = TRUE, svgp = TRUE, svgp_inducing = indN, 
				svgp_batch = batchN, svgp_iter=iterN))
			svgpTime <- unname(as.numeric(svgpTime[1] + svgpTime[3]))
			stats <- statsRow(stats, 'svgp', 'fit.time', ss, svgpTime)
			svgpPrTime <- system.time(svgpPr <- predict(svgpMod, grid_v$X))
			svgpPrTime <- unname(as.numeric(svgpPrTime[1] + svgpPrTime[3]))
			stats <- statsRow(stats, 'svgp', 'predict.time', ss, svgpPrTime)
			svgpPr <- merge(svgpPr, dissim_v_tall)
			svgpValid <- rmse(svgpPr$obs, svgpMod$y_rev_transform(svgpPr$fit))
			stats <- statsRow(stats, 'svgp', 'rmse', ss, svgpValid)

			## here - calculate conf interval accuracy

			saveRDS(svgpMod, paste0('res/svgpmod_', ss, '.rds'))	
			saveRDS(svgpPr, paste0('res/svgppr_', ss, '.rds'))	
		} else {
			svgpTime <- svgpMod <- NA
		}


	}

	saveRDS(stats, stFile)
}






		
# run MBM
		# plot(gdmPr[,1], gdmPr[,2], pch=16, cex=0.5, bty='n', xlim=c(0,1), ylim=c(0,1))
# abline(0,1, lty=2, col='gray')
# plot(mbmPr$obs, mbmMod$y_rev_transform(mbmPr$fit), pch=16, cex=0.5, bty='n', xlim=c(0,1), ylim=c(0,1))
# abline(0,1, lty=2, col='gray')

		
## save results
# ideally, save all models along with the stats