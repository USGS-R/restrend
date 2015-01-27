#' Print an Object
#'
#' Print the results of the regional seasonal Kendall trend test.
#'
#' @param x the object to be printed, must be the output from 
#'\code{regionalSeaken}.
#' @param digits the number of digits to use when printing numeric values.
#' @param ... further arguments for other methods.
#' @return The object is returned invisibly.
#' @note Three p-values are printed for the analysis. The raw p-value is
#'based only on the computed variance of S. The p-value corrected for serial
#'correlation includes the adjustement described in Hirsh and Slack. The
#'p-value corected for serial and spatial correlation also includes the
#'adjustment describe in Sprague and others.
#' @method print rsktest
#' @export
print.rsktest <- function(x, digits=4, ...) {
  ## Coding history:
  ##    2013Oct26 DLLorenz Original Coding
  ##
	cat("\n")
	cat(strwrap(x$method, prefix = "\t"), sep = "\n")
	cat("\ndata:  ", x$data.name, "\n", sep = "")
	cat("\nThe value of S is ", x$statistic[2L], ".", sep="")
	cat("\nThe tau correlation coefficient is ", signif(x$statistic[1L], digits), ".\n", sep="")
	cat("\nThe attained p-values:")
	minp <- format(10^-digits, scientific=8)
	pval <- round(x$p.value[1L], digits)
	if(pval > 0) {
		cat("\n                                          Raw = ", pval, sep="")
	} else
		cat("\n                                          Raw < ", minp, sep="")
	pval <- round(x$p.value[2L], digits)
	if(pval > 0) {
		cat("\n             Corrected for serial correlation = ", pval, sep="")
	} else
		cat("\n             Corrected for serial correlation < ", minp, sep="")
	pval <- round(x$p.value[3L], digits)
	if(pval > 0) {
		cat("\n Corrected for serial and spatial correlation = ", pval, sep="")
	} else
		cat("\n Corrected for serial and spatial correlation < ", minp, sep="")
	cat("\n\nThe median of the median slope at each site = ", 
			signif(x$estimate, digits), "\n\n", sep="")
  invisible(x)
}
