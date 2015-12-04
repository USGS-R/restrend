#' Handle Missing Values in Data
#' 
#' This function is useful for dealing with \code{NA}s in data frames. Only 
#'rows in the data specified by \code{columns} that have only missing values
#'are removed.
#' 
#' @details The value for \code{columns} can be a character vector containing the
#'name of the columns to check, or "All," which checks all columns. Or it can be
#'a function that returns a logical value--\code{TRUE} checks the column and 
#'\code{FALSE} does not. Any additonal arguments the the function can be given
#'by \dots.
#'
#' @param data the dataset to subset for all missing values.
#' @param columns the columns to check for missing values. See \bold{Details}.
#' @param \dots any additional arguments to \code{columns} if it is a funciton.
#' @return The dataset \code{data} having rows with at least one nonmissing value
#'in the columns specified by \code{columns}.
#' @examples
#' # create a short test dataset
#' test.df <- data.frame(A=c("a", "b", "c", "d", "e"),
#'    B=c(1, 2, NA, NA, 5), # numeric values
#'    C=c(1L, 3L, NA, 4L, NA), # integer values
#'    D=c(1, 2, NA, NA, 5)) # more numeric values
#' # The default, no row has all missing values
#' dropAllMissing(test.df)
#' # Check 2 columns
#' dropAllMissing(test.df, columns=c("B", "C"))
#' # check all numeric, including both integer and numeric types
#' dropAllMissing(test.df, columns=is.numeric)
#' # Check only those that are type numeric
#' dropAllMissing(test.df, columns=inherits, what="numeric")
#' @export
dropAllMissing <- function(data, columns="All", ...) {
	##
	## identify the columns
	if(is.character(columns)) {
		if(length(columns) == 1L && columns=="All") {
			picks <- rep(TRUE, ncol(data))
		} else {
			picks <- names(data) %in% columns
		}
	} else { # better be a function
		picks <- sapply(data, columns, ...)		
	}
	if(all(!picks)) {
		warning("No columns selected, returning the original data")
	} else {
		test <- rowSums(!is.na(data[picks]))
		data <- data[test > 0, ]
	}
	return(data)
}
