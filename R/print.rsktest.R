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
#'adjustment based on Dietz and Killeen (1981) described in Sprague and others
#'(2009). An alternative corrected value based on the methods described 
#'Douglas and others (2000) is also printed.
#' @references
#'Dietz, E.J., and Killeen, T.J., 1981, A nonparametric multivariate test for
#'monotone trend with pharmaceutical applications: Journal of the American
#'Statistical Association, v. 76, p 169--174.\cr
#'Douglas, E.M., Vogel, R.M., and Kroll, C.N., 2000, Trends in floods and low 
#'flows in the United States: impact of spatial correlation: Journal of 
#'Hydrology, v. 240, p. 90--105.\cr
#'Sprague, L.A., Mueller, D.K., Schwarz, G.E., and Lorenz, D.L., 2009,
#'Nutrient trends in streams and rivers of the United States, 1993--2003:
#'U.S. Geolgical Survey Scientific Investigations Report 2008-5202, 196 p.
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
	prnp <- format(pval, scientific=8)
	if(pval > 0) {
		cat("\n                                            Raw = ", prnp, sep="")
	} else
		cat("\n                                            Raw < ", minp, sep="")
	pval <- round(x$p.value[2L], digits)
	prnp <- format(pval, scientific=8)
	if(pval > 0) {
		cat("\n               Corrected for serial correlation = ", prnp, sep="")
	} else
		cat("\n               Corrected for serial correlation < ", minp, sep="")
	pval <- round(x$p.value[3L], digits)
	prnp <- format(pval, scientific=8)
	if(pval > 0) {
		cat("\n   Corrected for serial and spatial correlation = ", prnp, sep="")
	} else
		cat("\n   Corrected for serial and spatial correlation < ", minp, sep="")
	pval <- round(x$p.value[4L], digits)
	prnp <- format(pval, scientific=8)
	if(pval > 0) {
		cat("\n Alternative for serial and spatial correlation = ", prnp, sep="")
	} else
		cat("\n Alternative for serial and spatial correlation < ", minp, sep="")
	cat("\n\nThe median of the median slope at each site = ", 
			signif(x$estimate, digits), "\n\n", sep="")
  invisible(x)
}
