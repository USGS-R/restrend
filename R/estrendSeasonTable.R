#' Evaluate Seasonal Definitions
#' 
#' Create tables of observations for each season in both the BE and MI
#'records. Internal use only.
#'
#' @param ssn12 the vector of data selected for the 12 seasons per
#'year seaken analysis.
#' @param frctn limit for the fraction of possible observations in 
#'the BE.
#' @param pctg required percentage of seasons needed to exceed frctn.
#' @return A list of the results for 12, 6, 4, and 3 seasons per year
#'and the "best."
estrendSeasonTable <- function(ssn12, frctn=.5, pctg=.8) {
	## Coding history:
	##    2013Aug20 DLLorenz Original Coding
	##
	years.analyzed <- length(ssn12)/12
	years.be <- max(round((years.analyzed + 1)/5,0),2) # get at least 2 years
	years.mi <- years.analyzed - 2*years.be
	retval <- list()
	## 12 seasons per year
	seasons <- !is.na(ssn12)
  season.list <- rep(seq(12), years.analyzed)
  record.list <- c(rep("BE", 12*years.be), rep("MI", 12*years.mi),
                   rep("BE", 12*years.be))
  ## To absolutely guarantee a full table tack on a complete 
	## replication and then subtract 1
	season.list <- c(season.list[seasons],rep(seq(12), 2L))
	record.list <- c(record.list[seasons],rep(c("BE","MI"), each=12L))
	season.table <- table(record.list, season.list)
	season.table <- 100 * (season.table - 1) / 
		rep(c(2*years.be, years.mi), 12)
	retval[["s12"]] <- season.table
	g12 <- sum(season.table[1L,] > frctn * 100)/12 > pctg
	## 6 seasons per year
	ssn6 <- colSums(matrix(ssn12, nrow=2), na.rm=TRUE)
	seasons <- !is.na(ssn6)
	season.list <- rep(seq(6), years.analyzed)
	record.list <- c(rep("BE", 6*years.be), rep("MI", 6*years.mi),
									 rep("BE", 6*years.be))
	## Guarantee a full table
	season.list <- c(season.list[seasons],rep(seq(6), 2L))
	record.list <- c(record.list[seasons],rep(c("BE","MI"), each=6L))
	season.table <- table(record.list, season.list)
	season.table <- 100 * (season.table - 1) / 
		rep(c(2*years.be, years.mi), 6)
	retval[["s6"]] <- season.table
	g6 <- sum(season.table[1L,] > frctn * 100)/6 > pctg
	## 4 seasons per year
	ssn4 <- colSums(matrix(ssn12, nrow=3), na.rm=TRUE)
	seasons <- !is.na(ssn4)
	season.list <- rep(seq(4), years.analyzed)
	record.list <- c(rep("BE", 4*years.be), rep("MI", 4*years.mi),
									 rep("BE", 4*years.be))
	## Guarantee a full table
	season.list <- c(season.list[seasons],rep(seq(4), 2L))
	record.list <- c(record.list[seasons],rep(c("BE","MI"), each=4L))
	season.table <- table(record.list, season.list)
	season.table <- 100 * (season.table - 1) / 
		rep(c(2*years.be, years.mi), 4)
	retval[["s4"]] <- season.table
	g4 <- sum(season.table[1L,] > frctn * 100)/4 > pctg
	## 3 seasons per year
	ssn3 <- colSums(matrix(ssn12, nrow=4), na.rm=TRUE)
	seasons <- !is.na(ssn3)
	season.list <- rep(seq(3), years.analyzed)
	record.list <- c(rep("BE", 3*years.be), rep("MI", 3*years.mi),
									 rep("BE", 3*years.be))
	## Guarantee a full table
	season.list <- c(season.list[seasons],rep(seq(3), 2L))
	record.list <- c(record.list[seasons],rep(c("BE","MI"), each=3L))
	season.table <- table(record.list, season.list)
	season.table <- 100 * (season.table - 1) / 
		rep(c(2*years.be, years.mi), 3)
	retval[["s3"]] <- season.table
	g3 <- sum(season.table[1L,] > frctn * 100)/3 > pctg
	if(g12) {
		retval$best <- 12
	} else if(g6) {
		retval$best <- 6
	} else if(g4) {
		retval$best <- 4
	} else if(g3) {
		retval$best <- 3
	} else
		retval$best <- 0 # No can do!
	return(retval)
}

