#' Summarize Samples
#' 
#' Produce a summary of sample data by station in a dataset.
#' 
#' @param data the dataset to summarize.
#' @param DATES the name of the column containing the sample dates.
#' @param STAID the name of the column containing the station identifiers.
#' @param by.numeric compute summaries for each numeric column in \code{data}?
#' @return A data frame contining the starting and ending dates of the 
#'samples and the number of samples by station identifier if \code{by.numeric}
#'is \code{FALSE}. If \code{by.numeric} is \code{TRUE}, then the returned 
#'data are by station and numeric column (Reponse) and an indicator of 
#'censoring is included.
#' @examples
#' # do something here
#' @export
sumSamp <- function(data, DATES="DATES", STAID="STAID", by.numeric=TRUE) {
	## Coding history:
	##    2013Aug19 DLLorenz Original Coding
	##
	if(by.numeric) {
		numerics <- sapply(data, is.numeric)
		numerics <- names(numerics)[numerics] # Get the column names
		retval <- lapply(numerics, function(i) {
			ret1 <- sumSamp(data[!is.na(data[[i]]), ], DATES=DATES, STAID=STAID, 
										 by.numeric=FALSE)
			ret2 <- tapply(data[[i]], data[[STAID]], censoring)
			## Need to subset ret2 for all missing
			ret2 <- ret2[names(ret2) %in% ret1$STAID]
			ret <- cbind(Response=i, ret1, Censoring=ret2)
		} )
		retval <- do.call(rbind, retval)
	} else {
		beg <- tapply(data[[DATES]], data[[STAID]], min, na.rm=TRUE)
		end <- tapply(data[[DATES]], data[[STAID]], max, na.rm=TRUE)
		class(beg) <- class(end) <- class(data[[DATES]]) # Recover that info
		Nsamp <- tapply(data[[DATES]], data[[STAID]], function(x) sum(!is.na(x)))
		retval <- data.frame(STAID=names(Nsamp), FirstSamp=beg, 
												 LastSamp=end, NumSamp=Nsamp, 
												 row.names=seq(length(Nsamp)))
	}
	return(retval)
}
