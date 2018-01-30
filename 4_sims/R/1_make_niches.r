# script to generate a bunch of multivariate gaussian niches 
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(raster)
library(viridis)

make_niche <- function(center_means = c(0,0), center_sig = diag(2), 
	width_shapes = c(1.2, 1.2), width_rates = c(1.2, 1.2), rho_ab = c(0,0),
	scale_ab = c(1.2, 1.2))
# for width parameters, first entry in each describes width for env1, second
# for env2
# for rho params, rho ~ 2 * Beta(rho_ab[1], rho_ab[2]) - 1
# for the special case of rho = 0,0, covariance will be zero
{
	# vector of two, the niche in each dimension
	cntrs <- rmvnorm(1, center_means, center_sig)
	# niche widths
	widths <- rgamma(2, width_shapes, width_rates)
	# niche axis correlations
	if(all(rho_ab == 0)) {
		rho <- 0
	} else {
		rho <- 2 * rbeta(1, rho_ab[1], rho_ab[2]) - 1
	}
	# how high is the niche; controls abundance and detectability of the species
	scale <- rbeta(1, scale_ab[1], scale_ab[2])

	data.frame(center1 = cntrs[1], center2 = cntrs[2], width1 = widths[1], 
		width2=widths[2], rho=rho, scale=scale)
}

niche_height <- function(x, niche)
{
	# computes the niche height as the multivariate normal density times the scale
	mu <- niche[1:2]
	sig <- matrix(NA, nrow=2, ncol=2)
	sig[1,1] <- niche[3]
	sig[2,2] <- niche[4]
	sig[1,2] <- sig[2,1] <- niche[3]*niche[4]*niche[5]
	sc <- niche[6]
	# scale height to one
	ht <- dmvnorm(x, mu, sig) / dmvnorm(mu, mu, sig)
	sc * ht
}
# now use the function above to generate a list of n niches
# for covariance, best to use large params for rho; 100,100 makes for just a bit
# if the rho params are equal, mean correlation will be 0; the larger the numbers the
# smaller the variance
niches <- do.call(rbind, lapply(1:300, function(x) make_niche(center_sig=12*diag(2),
	width_shapes = c(4,7), width_rates = c(1.2,1.2), rho_ab=c(100,100), scale_ab=c(15,5))))

# make a figure with histogram (or density) of each parameter 
#		(i.e., 2 centers, 4 widths, 1 scale)
# niches_tall <- melt(niches)
# ggplot(niches_tall, aes(x = value)) + geom_density() + facet_wrap(~variable, ncol=4)

#generate a matrix as a test landscape
lscape <- expand.grid(env1=seq(-5,5,length.out=50), env2=seq(-5,5,length.out=50))
# for each cell the in the matrix, compute the sum of the niche heights for those env 
#		conditions; the value is the expected species richness
site_sp_pr <- apply(niches, 1, function(ni) niche_height(lscape, ni))
richness_exp <- rowSums(site_sp_pr)

site_sp <- matrix(rbinom(length(site_sp_pr), 1, site_sp_pr), ncol=ncol(site_sp_pr))
richness <- rowSums(site_sp)

# plot it using image, also examine mean, max, min, range of the matrix
richRasExp <- rasterFromXYZ(cbind(lscape, richness_exp))
richRas <- rasterFromXYZ(cbind(lscape, richness))
par(mfrow=c(2,3))
plot(richRasExp, col=magma(100), zlim=c(0,ceiling(maxValue(richRasExp))))
plot(lscape$env1, richness_exp, cex=0.5, pch=16)
plot(lscape$env2, richness_exp, cex=0.5, pch=16)

plot(richRas, col=magma(100), zlim=c(0,ceiling(maxValue(richRas))))
plot(lscape$env1, richness, cex=0.5, pch=16)
plot(lscape$env2, richness, cex=0.5, pch=16)


## if you want to view the relatinoship looking at all plots
# xDF <- env_dissim(lscape)
# sor <- sorensen(site_sp, FALSE)
# rownames(sor) <- colnames(sor) <- 1:2500
# yDF <- reshape2::melt(sor,varnames=c('site1', 'site2'), value.name = 'sor')
# yDF[,1] <- as.character(yDF[,1])
# yDF[,2] <- as.character(yDF[,2])
# dat <- merge(xDF, yDF, all.x = TRUE)
# plot(dat$distance, dat$sor, pch=16, cex=0.3)
