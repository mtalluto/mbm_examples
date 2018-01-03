# script to generate a bunch of multivariate gaussian niches 
library(MASS)

make_niche <- function(center_means = c(0,0), center_sig = diag(2), 
	width_shapes = c(1.2, 0, 0, 1.2), width_rates = c(1.2, 0, 0, 1.2),
	scale_alpha = 1.2, scale_beta = 1.2)
{
	list(
		center = mvrnorm(1, center_means, center_sig), # vector of two, the niche in each dimension
		widths = matrix(c(rgamma(1, width_shapes[1], width_rates[1]), 
			0, 0, rgamma(1, , width_shapes[2], width_rates[2])), nrow=2), # variance covariance matrix of this niche; for now assume no covariance
		scale = rbeta(1, scale_alpha, scale_beta)# how high is the niche; controls abundance and detectability of the species
	)
}


# now use the function above to generate a list of n niches
# make a figure with histogram (or density) of each parameter (i.e., 2 centers, 4 widths, 1 scale)
# also generate a matrix as a test landscape; maybe start from -15 to 15 in each dimension with a resolution of 1
# for each cell the in the matrix, compute the sum of the niche heights for those env conditions; the value is the expected species richness
# plot it using image, also examine mean, max, min, range of the matrix
# repeat using a finer and finer grid; the goal is to have no 'gaps' or other oddities in environmental space; species richness gradient could be smooth
# if not, increase the number of niches generated or tweak the hyperparameters