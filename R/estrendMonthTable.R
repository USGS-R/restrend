#' Evaluate Monthly Definitions
#' 
#' Create tables of observations for each season in both the first and last
#'halfs of the record. Internal use only.
#'
#' @param ssn12 the vector of data selected for the monthly
#'seaken analysis.
#' @param frctn limit for the fraction of possible observations in 
#'each half.
#' @return A list of the results for the percentage of observations in each
#'half, the selected months, and the number of selected months.
estrendMonthTable <- function(ssn12, frctn=.5) {
	## Coding history:
	##    2015Nov18 DLLorenz Original Coding
	##
	years.analyzed <- length(ssn12)/12
	## split the record into 1/2s, skipping middle if odd
	first.half <- seq(1, years.analyzed/2)
	last.half <- years.analyzed + 1 - rev(first.half)
	half.len <- length(first.half)
	target <- frctn * 100
	# Convert to matrix and sum number of samples per month
	ssn12 <- matrix(ssn12, nrow=12)
	months.list <- list(first=rowSums(!is.na(ssn12[, first.half])) / half.len * 100,
											last=rowSums(!is.na(ssn12[, last.half])) / half.len * 100)
	months.pick <- months.list$first >= target & months.list$last >= target
	return(list(months.list=months.list, months.pick=months.pick, nmons=sum(months.pick)))
}

