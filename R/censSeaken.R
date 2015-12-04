#' Trend Test
#'
#' Compute the seasonal Kendall trend test
#'with the Turnbull slope estimator for left-censored data.
#'
#' @param series any regularly spaced object that can be forced to 
#'class "lcens" to test for trend. Missing values are permitted.
#' @param nseas the number of seasons per year.
#' @return An object of class "htest" also inhereting class "seaken" 
#' containing the following components: 
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
#'  \item{nyears}{the number of years.}
#'  \item{nseasons}{the number of seasons.}
#'  \item{series}{the data that was analyzed converted to numeric values.}
#' @note The value of \code{p.value} is \code{p.value.raw} if there are 
#'fewer than 10 years of data and is \code{p.value.corrected} otherwise.
#' @references The approach used in \code{censSeaken} was used ooriginally in Sullivan 
#'and others (2009). The original version of the code was published in Lorenz 
#'and others (2011). It is based on principles for comparing censored values
#'in Helsel (2012) and the Turnbull slope estiamte is described by Turnbul(1974).
#'
#'Helsel, D.R. 2012, Statistics for censored environmental data 
#'using Minitab and R: New York, Wiley, 324 p.\cr
#'Lorenz, D.L., Ahearn, E.A., Carter, J.M., Cohn, T.A., Danchuk, W.J., 
#'Frey, J.W., Helsel, D.R., Lee, K.E., Leeth, D.C., Martin, J.D., 
#'McGuire, V.L., Neitzert, K.M., Robertson, D.M., Slack, J.R., Starn, J., 
#'Vecchia, A.V., Wilkison, D.H., and Williamson, J.E., 2011, 
#'USGS library for S-PLUS for Windows---Release 4.0: U.S. Geological 
#'Survey Open-File Report 2011-1130. 
#'(Available at \url{http://pubs.er.usgs.gov/publication/ofr20111130}).\cr
#'Sullivan, D.J., Vecchia, A.V., Lorenz, D.L., Gilliom, R.J., and 
#'Martin, J.D., 2009, Trends in pesticide concentrations in corn-belt streams,
#'1996???-2006: U.S. Geological Survey Scientific Investigations Report 
#'2009-5132, 75 p.\cr
#'Turnbull, B.W., 1974, Nonparametric estimation of a survivorship function 
#'with doubly censored data: Journal of the American Statistical Society, 
#'v. 69, p. 169--173.
#' @examples
#'\dontrun{
#'# Compare censored and uncensored to seaken
#'library(USGSwsData)
#'data(KlamathTP)
#'# Construct the regular series
#'KlamathTP.rs <- with(KlamathTP, regularSeries(TP_ss, sample_dt, 
#'  begin="1972-01-01", end="1980-01-01"))
#'# Uncensored, differences due to rounding
#'with(KlamathTP.rs, seaken(Value, 12))
#'with(KlamathTP.rs, censSeaken(Value, 12))
#'# About 30 percent censoring, censSeaken closer to uncensored slope
#'with(KlamathTP.rs, seaken(ifelse(Value < 0.05, 0.025, Value), 12))
#'with(KlamathTP.rs, censSeaken(as.lcens(Value, 0.05), 12))
#'}
#' @useDynLib restrend seakenlcd kendalllcd
#' @export
censSeaken <- function(series, nseas=12) {
  ## Coding history:
  ##    2013Sep06 DLLorenz Original Coding from S+
  ##
  ## Error checking.
  ## Make sure we have complete data
	sernam <- deparse(substitute(series))
	series <- as.lcens(series)
	nobs <- length(series)
	nseas <- as.integer(nseas)
  nyrs <- as.integer(nobs / nseas)
  if(nyrs < 2L)
    stop("Must have at least 2 years of data")
	if(nobs %% nseas != 0L) {
		nfull <- ceiling(nobs/nseas)*nseas
		series <- c(series, rep(NA_real_, nfull - nobs))
		warning(paste("The original series of", nobs,
									"values was padded with missing values to a length of", nfull,
									"so that it contains an integral number of years."))
		nobs <- as.integer(nfull)
		nyrs <- as.integer(nobs / nseas)
	}
  ## construct t (xx)
  xx <- rep(1:nyrs, each=nseas)
  sel <- !is.na(series)
  ## Median of the data values.
  median.data <- median(series, na.rm=TRUE)
  ## Median of the time values
  median.time <- nyrs / 2.0 # as defined in seaken
  ## Construct matrices of yy, cy, xx, cx, sel
  yy <- matrix(series@.Data[, 1L], ncol=nseas, byrow=TRUE)
  cy <- na2miss(as.double(series@censor.codes), -1)
  cy <- matrix(cy, ncol=nseas, byrow=TRUE)
  ycadj <- min(diff(unique(sort(yy))), na.rm=TRUE)/1000.
  yy <- yy - ycadj * cy # y-values adjusted a small amount down
  xx <- matrix(xx, ncol=nseas, byrow=TRUE)
  cx <- matrix(0, nrow=nyrs, ncol=nseas)
  sel <- matrix(sel, ncol=nseas, byrow=TRUE)
  ## compute s, vars, tau, tau-b (if desired)
  results <- rep(0.0, 6L)
  slopemaxy <- double(0)
  slopeminy <- double(0)
  for(i in 1:nseas) {
    selin <- sel[, i]
    yyin <- yy[selin, i]
    cyin <- cy[selin, i]
    xxin <- xx[selin, i]
    cxin <- cx[selin, i]
    nin <- length(yyin)
    if(nin > 1) {
      nsl <- max(nin,nin*(nin-1)/2)
      out <- .Fortran("kendalllcd",
                      xx = as.double(xxin),
                      cx = as.double(cxin),
                      yy = as.double(yyin),
                      cy = as.double(cyin),
                      n = as.integer(nin),
                      results = double(6),
                      slopel = double(nsl),
                      slopeh = double(nsl),
                      nsl = nsl)
      results <- results + out$results
      slopemaxy <- c(slopemaxy, out$slopeh[1:out$nsl])
      slopeminy <- c(slopeminy, out$slopel[1:out$nsl])
    } # end if nin > 1
  } # end loop
  ## compute the adjustment for serial dependence
  adj <- c(0, 0)
  yy <- na2miss(yy, -99)
  cy <- ifelse(yy < -98, -1, cy) # flag for missing values
  if(nseas > 1L) { # adjust only if number of seasons > 1
    for(i in 1:(nseas-1)) {
      for(j in (i+1):nseas) {
        yyin <- yy[,i]
        cyin <- cy[,i]
        xxin <- yy[,j]
        cxin <- cy[,j]
        nin <- length(yyin)
        out <- .Fortran("seakenlcd",
                        xx = as.double(xxin),
                        cx = as.integer(cxin),
                        yy = as.double(yyin),
                        cy = as.integer(cyin),
                        n = as.integer(nin),
                        results = double(2))
        adj <- adj + out$results * 2
      }}
  } # end if
  losl <- pmin(slopemaxy, slopeminy)
  hisl <- pmax(slopemaxy, slopeminy)
  ## compute median slope
	helselslope <- median(as.mcens(losl, hisl))
  method <- "Seasonal Kendall's tau with the Turnbull slope estimator"
  tau <- results[1L]/results[2L]* 2.0
  kenS <- results[1L]
  if(is.na(results[6L])) # protect against certain unusual conditions
    varS <- results[5L]
  else
    varS <- results[5L] - results[6L]
  p.val <- 2 * (1 - pnorm((abs(kenS - sign(kenS)))/ sqrt(varS)))
  varS <- varS + sum(adj)/3.
  p.val.c <- 2 * (1 - pnorm((abs(kenS - sign(kenS)))/ sqrt(varS)))
  names(tau) <- "tau"
  zero <- 0
  names(zero) <- "slope"
  est <- c(helselslope, median.data, median.time)
  names(est) <- c("slope", "median data", "median time")
  if(nyrs < 10)
    p.value <- p.val
  else
    p.value <- p.val.c
	# Modify series to be numeric, needed for plotting 
	# mostly for compatibility with seaken!
	# Find the smallest integer value less than the smallest value
	# and set the simple substitution value to that
	minser <- min(series@.Data[,1L], na.rm=TRUE)
	offset <- (minser - trunc(minser - 0.0001))/2
	series <- ifelse(series@censor.codes, series@.Data[,1L] - offset,
									 series@.Data[,1L])
  retval <- list(method=method,  
  							 data.name=paste(sernam, " (",nyrs, " years and ", nseas, " seasons)", sep=""),
                 nyears=nyrs, nseasons=nseas, series=series,
                 statistic=tau, p.value=p.value,
                 p.value.raw=p.val,
                 p.value.corrected=p.val.c,
                 estimate=est, alternative="two.sided",
                 null.value = zero)
  oldClass(retval) <- c("htest", "seaken")
  return(retval)
}
