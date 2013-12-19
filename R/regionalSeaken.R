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
#'  \item{p.value}{the p-value. See \bold{Note}.}
#'  \item{p.value.raw}{the p-value computed without correction for serial 
#'correlation. See \bold{Note}.}
#'  \item{p.value.corrected}{the p-value computed with correction for 
#'serial correlation. See \bold{Note}.}
#'  \item{estimate}{a named vector containing the Sen estimate of the 
#'slope in units per year, the median value of the data, and the median 
#'value of time.}
#'  \item{data.name}{a string containing the actual name of the input 
#'series with the number of years and seasons.}
#'  \item{alternative}{a character string describing alternative to the
#'test ("two.sided").}
#'  \item{null.value}{the value for the hypothesized slope (0).}
#' @note The value of \code{p.value} is \code{p.value.raw} if there are 
#'fewer than 10 years of data and is \code{p.value.corrected} otherwise.
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
	zscores <- S/sqrt(cumsum(results$results[2L:4L]))
	p.values <- (1 - pnorm(abs(zscores)))*2 # two-sided test
	est = c(slope=results$slope)
	data.name = paste(Dname, " (", nsite, " sites,",
										nyrs, " years, and ", nseas, " seasons)", sep = "")
	retval <- list(method=method, data.name=data.name, 
								 statistic=stats,
								 p.value=p.values, estimate=est, results=results)
  oldClass(retval) <- c("rsktest")
  return(retval)
}
