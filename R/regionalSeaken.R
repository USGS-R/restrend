#' Trend Test
#'
#' Compute the regional seasonal Kendall trend test.
#'
#' @param series a matrix of the regular series for each site. Each
#'column is the regular series of observations for that site.
#'Missing values are permitted.
#' @param nseas the number of seasons per year.
#' @return An object of class "rsktest" containing the following components: 
#'  \item{method}{a description of the method.}
#'  \item{statistic}{the value of Kendall's tau.}
#'  \item{p.values}{the p-value. See \bold{Note}.}
#'  \item{estimate}{the Sen estimate of the 
#'slope in units per year--the median of the site medians.}
#'  \item{data.name}{a string containing the actual name of the input 
#'series with the number of years and seasons.}
#' @note The values of \code{p.values} are the attained p-values considering
#'the raw results, corrected for serial correlation, corrected for spatial
#'correlation using the method of Dietz and Killeen (1981), and the corrected 
#'value for spatial correlation using the method of Douglas and others (2000).
#' @references
#'Dietz, E.J., and Killeen, T.J., 1981, A nonparametric multivariate test for
#'monotone trend with pharmaceutical applications: Journal of the American
#'Statistical Association, v. 76, p 169--174.\cr
#'Douglas, E.M., Vogel, R.M., and Kroll, C.N., 2000, Trends in floods and low 
#'flows in the United States: impact of spatial correlation: Journal of 
#'Hydrology, v. 240, p. 90--105.
#'Sprague, L.A., Mueller, D.K., Schwarz, G.E., and Lorenz, D.L., 2009,
#'Nutrient trends in streams and rivers of the United States, 1993--2003:
#'U.S. Geolgical Survey Scientific Investigations Report 2008-5202, 196 p.
#' @examples
#'\dontrun{
#'data(RK3b)
#'# Build matrix, site 11 is missing 2004, would be row 60
#'Ammonia <- with(RK3b, matrix(c(value[1:59], NA, value[60:299]), 
#'  ncol=25))
#'regionalSeaken(Ammonia, 1)
#'}
#' @useDynLib restrend regionsk
#' @export
regionalSeaken <- function(series, nseas=12) {
  ## Coding history:
  ##    2013Oct25 DLLorenz Original Coding
  ##
  ## Preliminary processing
	Dname <- deparse(substitute(series))
	series <- as.matrix(series)
	nobs <- nrow(series)
	nseas <- as.integer(nseas)
  nyrs <- as.integer(nobs / nseas)
	nsite <- ncol(series)
  if(nyrs < 2L)
    stop("Must have at least 2 years of data")
	if(nobs %% nseas != 0L) {
		stop(paste("The original series is not an even multiple of the nseas."))
	}
	series <- na2miss(series, -99999)
  ## Do it
	results <- .Fortran("regionsk", as.single(series), nsite, nseas, nyrs,
											results=single(4), scomp=single(1), slope=single(1))
	method <- "Regional Seasonal Kendall Test for Trend"
	stats <- c(tau=results$results[1L]/results$scomp, S=results$results[1L])
	S <- results$results[1L] - sign(results$results[1L]) # Continuity corr.
	# Compute p.value using the method of Douglas, Vogel and Kroll
	x.cor <- cor(series,  method="spearman", use="pair")
	rhoxx <- 2*mean(x.cor[lower.tri(x.cor)])
	NewVarS <- sum(results$results[2L:3L])*(1 + rhoxx)
	zscores <- S/sqrt(c(cumsum(results$results[2L:4L]), NewVarS))
	p.values <- (1 - pnorm(abs(zscores)))*2 # two-sided test
	est = c(slope=results$slope)
	# Put it together
	data.name = paste(Dname, " (", nsite, " sites, ",
										nyrs, " years, and ", nseas, " seasons)", sep = "")
	retval <- list(method=method, data.name=data.name, 
								 statistic=stats,
								 p.value=p.values, estimate=est, results=results,
								 spatialcorr=rhoxx)
  oldClass(retval) <- c("rsktest")
  return(retval)
}
