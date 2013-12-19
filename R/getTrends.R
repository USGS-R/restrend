#' Get Trends
#' 
#' Extract the trend results.
#' 
#' @param Stations a vector of the station identifiers on which 
#'to do the flow adjustment.
#' @param Snames a vector of the response variables on which to do 
#'the flow adjustment.
#' @param sig.level the alpha level of the test. If the attained p-value
#'of the test is less than \code{sig.level}, then the null hypothesis of
#'no trend is rejected and the trend is classes as "up" or "down."
#' @return A data frame of the requested trend tests having these columns:
#'\item{Station}{ the station identifier.}
#'\item{Response}{ the response variable name.}
#'\item{Type}{ the type of trend test.}
#'\item{NumYears}{ the number of years for the trend test.}
#'\item{NumSeasons}{ the number of seasons. Set to \code{NA} for the Tobit
#'trend test.}
#'\item{Nobs}{ the number of observations used in the test.}
#'\item{RepValue}{ the representative value of the response variable. See
#'\code{Note}.}
#'\item{Trend}{ the trend expressed as the average rate in Response units
#'per year}
#'\item{Trend.pct}{ the trend expressed as the average rate in percent
#'change per year}
#'\item{P.value}{ the attained p-value of the test.}
#'\item{Trend.dir}{ an indicator of the trend. \code{Trend.dir} is "none" if the
#'attained p-value is greater than \code{sig.level}, "up" if the
#'attained p-value is less than \code{sig.level} and \code{Trend} is positive,
#'and "down" if the attained p-value is less than \code{sig.level} and 
#'\code{Trend} is negative. For the seasonal Kendall test, when there are many 
#'tied values, the results from the trend test and the trend may not agree,
#'in these cases \code{Trend.dir} will be "*."}
#' @note The representative value is the median from a sample of the data
#'taken to represent a reasonably uniform sampling over time and throughout
#'the year. For the seasonal Kendall test, it is the median of the data, for 
#'the Tobit test, it is the median of the data sampled as for a seasonal 
#'Kendall test.
#' @export
getTrends <- function(Stations="All", Snames="All", sig.level=0.05) {
	## Coding history:
	##    2013Oct21 DLLorenz Original Coding
	##
	## Preliminaries
	pos <- ckProj()
	estrend.in <- get("estrend.in", pos=pos)
	estrend.st <- get("estrend.st", pos=pos)
	if(Stations[1L] == "All")
		Stations <- estrend.in$stations
	if(Snames[1L] == "All")
		Snames <- estrend.in$snames
	retval <- list()
	i <- 0L
	## Seasonal Kendall
	if(estrend.in$type == "seasonal") {
		estrend.sk <- get("estrend.sk", pos=pos)
		for(station in Stations) {
			for(sname in Snames) {
				if(estrend.st[station, sname] == "OK") {
					results <- estrend.sk[station, sname][[1L]]
					restyps <- names(results)
					reslog <- substring(results[[1L]]$data.name, 1L, 4L) == "log("
					RepValue <- results[[1L]]$estimate[2L]
					if(reslog)
						RepValue <- exp(RepValue)
					Type <- c(SKu="uncensored seasonal Kendall",
										SKFAC="flow-adjusted seasonal Kendall",
										SKc="censored seasonal Kendall")[restyps]
					NumYears <- sapply(results, function(x) x$nyear)
					NumSeas <- sapply(results, function(x) x$nseasons)
					Nobs <- sapply(results, function(x) sum(!is.na(x$series)))
					Tnd <- sapply(results, function(x) x$estimate[1L])
					if(reslog) {
						## Tnd is the trend in logs
						Trend.pct <- 100*(exp(Tnd) - 1)
						Trend <- Trend.pct/100*RepValue
					} else {
						## Tnd is the trend in actual units
						Trend <- Tnd
						Trend.pct <- 100*Trend/RepValue
					}
					P.value <- sapply(results, function(x) x$p.value)
					Result <- ifelse(P.value > sig.level, "none",
													 ifelse(Trend > 0, "up", "down"))
					Cktnd <- sapply(results, function(x) x$statistic)
					## Error check on the direction of trend
					Result <- ifelse(sign(Trend) == sign(Cktnd), Result, "*")
					## Pack it up
					i <- i + 1L
					retval[[i]] <- data.frame(Station=station, Response=sname,
																	Type=Type, NumYears=NumYears,
																	NumSeas=NumSeas, Nobs=Nobs,
																	RepValue=RepValue, Trend=Trend,
																	Trend.pct=Trend.pct, P.value=P.value,
																	Trend.dir=Result, stringsAsFactors=FALSE)
				}
			}
		}
	} else if(estrend.in$type == "tobit") {
		estrend.tb <- get("estrend.tb", pos=pos)
		estrend.df <- get("estrend.df", pos=pos)
		DATES <- estrend.in$DATES
		for(station in Stations) {
			for(sname in Snames) {
				if(estrend.st[station, sname] == "OK") {
					results <- estrend.tb[station, sname][[1L]]$TB
					## determine the rep value
					temp.df <- estrend.df[station, sname][[1L]]
					peryear <- table(year(temp.df[[DATES]]))
					medper <- median(peryear)
					if(medper > 10) {
						k.per <- 1
					} else if(medper >= 5) {
						k.per <- 2
					} else if(medper >= 4) {
						k.per <- 3
					} else
						k.per <- 4
					## Create an index rather than the qw object
					temp.sel <- regularSeries(seq(nrow(temp.df)), temp.df[[DATES]],
																		k.period=k.per)$Value
					RepValue <- median(temp.df[[sname]][temp.sel], na.rm=TRUE)
					reslog <- substring(deparse(results$call$formula), 1L, 4L) == "log("
					NumYears <- diff(range(temp.df[[DATES]]))
					if(attr(NumYears, "units") == "days") {
						NumYears <- as.vector(NumYears)/365.25
					} else # seconds for POSIXt
						NumYears <- as.vector(NumYears)/22791600
					NumSeas <- NA_real_
					Nobs <- results$NOBSC
					NPAR <- results$NPAR
					Tnd <- results$PARAML[NPAR]
					if(reslog) {
						## Tnd is the trend in logs
						Trend.pct <- 100*(exp(Tnd) - 1)
						Trend <- Trend.pct/100*RepValue
					} else {
						## Tnd is the trend in actual units
						Trend <- Tnd
						Trend.pct <- 100*Trend/RepValue
					}
					P.value <- results$PVAL[NPAR]
					Result <- ifelse(P.value > sig.level, "none",
													 ifelse(Trend > 0, "up", "down"))
					## Pack it up
					i <- i + 1L
					retval[[i]] <- data.frame(Station=station, Response=sname,
																		Type="Tobit", NumYears=NumYears,
																		NumSeas=NumSeas, Nobs=Nobs,
																		RepValue=RepValue, Trend=Trend,
																		Trend.pct=Trend.pct, P.value=P.value,
																		Trend.dir=Result, stringsAsFactors=FALSE)
					
				}
			}
		}	
	} else { # Must be annual
		estrend.an <- get("estrend.an", pos=pos)
		estrend.df <- get("estrend.df", pos=pos)
		DATES <- estrend.in$DATES
		for(station in Stations) {
			for(sname in Snames) {
				if(estrend.st[station, sname] == "OK") {
					temp.df <- estrend.df[station, sname][[1L]]
					results <- estrend.an[station, sname][[1L]]
					restyps <- names(results)
					reslog <- substring(results[[1L]]$data.name, 1L, 4L) == "log("
					RepValue <- results[[1L]]$estimate[2L]
					if(reslog)
						RepValue <- exp(RepValue)
					Type <- c(aK="annual Kendall",
										aKFAC="flow-adjusted annual Kendall")[restyps]
					NumYears <- as.integer(diff(range(dectime(temp.df[[DATES]]))) + 1.1) # Best guess
					NumSeas <- 1L
					Nobs <- sapply(results, function(x) as.integer(x$estimate[4L]))
					Tnd <- sapply(results, function(x) x$estimate[1L])
					if(reslog) {
						## Tnd is the trend in logs
						Trend.pct <- 100*(exp(Tnd) - 1)
						Trend <- Trend.pct/100*RepValue
					} else {
						## Tnd is the trend in actual units
						Trend <- Tnd
						Trend.pct <- 100*Trend/RepValue
					}
					P.value <- sapply(results, function(x) x$p.value)
					Result <- ifelse(P.value > sig.level, "none",
													 ifelse(Trend > 0, "up", "down"))
					## Pack it up
					i <- i + 1L
					retval[[i]] <- data.frame(Station=station, Response=sname,
																		Type=Type, NumYears=NumYears,
																		NumSeas=NumSeas, Nobs=Nobs,
																		RepValue=RepValue, Trend=Trend,
																		Trend.pct=Trend.pct, P.value=P.value,
																		Trend.dir=Result, stringsAsFactors=FALSE)
				}
			}
		}
	}
	retval <- do.call(rbind, retval)
	row.names(retval) <- NULL
	return(retval)
}
