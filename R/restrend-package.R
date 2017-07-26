#'Estimate Trends (ESTREND)
#'
#'This package has specialized functions for managing data to 
#'facilitate testing for linear or monotonic trends in hydrologic data.
#'
#'\tabular{ll}{ Package: \tab restrend\cr 
#'Type: \tab Package\cr 
#'License: \tab File LICENSE\cr 
#'Depends: \tab g.data, smwrBase, smwrGraphs, smwrStats, smwrQW\cr }
#'This package contains functions that facilitate testing for linear or monotonic trends 
#'in hydrologic data. Water-quality data or any other data collected on a
#'nearly regular basis can be uncensored, left censored, or multiply censored.
#'Data for annual analysis must be uncensored.
#'
#' @name restrend-package
#' @aliases restrend-package restrend
#' @docType package
#' @author Dave Lorenz 
#'
#' @seealso \code{\link[smwrQW:smwrQW-package]{smwrQW}}
#' @import g.data dataRetrieval lubridate smwrBase smwrGraphs smwrStats smwrQW
#' @references Lorenz, D.L., in preparation, restrend---an R package for trend
#'estimation in hydrologic data, version 0.4.1
#' @keywords package
NULL
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This information is preliminary or provisional and
is subject to revision. It is being provided to meet
the need for timely best science. The information
has not received final approval by the U.S. Geological
Survey (USGS) and is provided on the condition that
neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from the authorized
or unauthorized use of the information.

****Orphaned Package****
This package is looking for a new maintainer. For more information, 
see: https://owi.usgs.gov/R/packages.html#orphan")
}
