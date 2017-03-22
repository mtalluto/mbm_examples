rmse <- function(pr, dat)
{
	# compute root mean square error between a prediction set and a data set
	sqrt(mean((pr -dat)^2))
}

descale <- function(x, scaling)
{
	for(rn in rownames(scaling))
	{
		if(rn %in% colnames(x))
		{
			x[,rn] <- (x[,rn] * scaling[rn,'scale']) + scaling[rn,'center']
		} else warning("Variable ", rn, " not found in input data")
	}
	x
}

resp.curve <- function(x, y, lower, upper, xyDat, xyValid, link, col='#6d81b1', alpha="66", add=FALSE, ...)
{
	if(missing(link)) link <- function(x) x
	if(add) {
		lines(x, link(y), col=col, ...)
	} else {
		plot(x, link(y), col=col, type='l', ...)
	}
	if(! missing(lower) & !missing(upper))
		polygon(c(x, rev(x)), c(link(lower), rev(link(upper))), border=NA, col=paste0(col, alpha))
	if(! missing(xyDat))
		points(xyDat[,1], xyDat[,2], pch=20, cex=0.5)
	if(! missing(xyValid))
		points(xyValid[,1], xyValid[,2], pch='+', col='red', cex=0.7)
}
