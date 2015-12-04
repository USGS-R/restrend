#' Regional Seasonal Kendall Trends
#' 
#' Perform the regional seasonal Kendall trend test. The test is only
#'appropriate for uncensored data and is performed after the trends have
#'been computed.
#' 
#' @param Stations a vector of the the station identifiers on which 
#'to do the trend test.
#' @param Sname the response variables on which to do 
#'the trend test.
#' @param FAC logical, if \code{TRUE}, then do the regional test on the
#'flow ajusted values, otherwise the analysis in done on the raw values.
#' @return An object of class "rsktest" from \code{regionalSeaken}.
#' @seealso \code{\link{regionalSeaken}}
#' @export
RSKTrends <- function(Stations="All", Sname, FAC=FALSE) {
	##
	## Preliminaries
	## 
	if(missing(Sname) || length(Sname) > 1L) {
		stop("must specify exactly one constituent")
	}
	pos <- ckProj()
	estrend.in <- get("estrend.in", pos=pos)
	if(estrend.in$analysis != "regular") {
		stop("the anaylsis must be regular")
	}
	estrend.st <- get("estrend.st", pos=pos)
	if(Stations[1L] == "All")
		Stations <- estrend.in$stations
	# Build the regional matrix
	ck <- FALSE
	if(estrend.in$type == "seasonal") {
		estrend.sk <- get("estrend.sk", pos=pos)
		for(station in Stations) {
			if(estrend.st[station, Sname] == "OK") {
				results <- estrend.sk[station, Sname][[1L]]
				results <- if(FAC) results$SKFAC else results$SKu
				if(is.null(results)) {
					warning("No results for station ", station)
				} else { # do it
					if(ck) { # after first OK
						if(nseas != results$nseasons) {
							stop("Mush use consistent season definition")
						}
						Region <- cbind(Region, results$series)
					} else {
						Region <- results$series
						nseas <- results$nseasons
						ck <- TRUE
					}
				}
			}
		} # end of station loop
	} else if(estrend.in$type == "annual") {
		estrend.sk <- get("estrend.sk", pos=pos)
		nseas <- 1L
		Region <- NULL
		for(station in Stations) {
			if(estrend.st[station, Sname] == "OK") {
				results <- estrend.sk[station, Sname][[1L]]
				results <- if(FAC) results$SKFAC else results$SKu
				if(is.null(results)) {
					warning("No results for station ", station)
				} else { # do it
					Region <- cbind(Region, results$series)
				}
			}
		} # end of station loop
	} else {
		stop("the type must be either \"seasonal\" or \"annual\"")
	}
	# Do it
	retval <- regionalSeaken(Region, nseas)
	return(retval)
}
