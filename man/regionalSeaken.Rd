\name{regionalSeaken}
\alias{regionalSeaken}
\title{Trend Test}
\usage{
  regionalSeaken(series, nseas = 12)
}
\arguments{
  \item{series}{a matrix of the regular series for each
  site. Each column is the regular series of observations
  for that site. Missing values are permitted.}

  \item{nseas}{the number of seasons per year.}
}
\value{
  An object of class "rsktest" containing the following
  components: \item{method}{a description of the method.}
  \item{statistic}{the value of Kendall's tau.}
  \item{p.value}{the p-value. See \bold{Note}.}
  \item{p.value.raw}{the p-value computed without
  correction for serial correlation. See \bold{Note}.}
  \item{p.value.corrected}{the p-value computed with
  correction for serial correlation. See \bold{Note}.}
  \item{estimate}{a named vector containing the Sen
  estimate of the slope in units per year, the median value
  of the data, and the median value of time.}
  \item{data.name}{a string containing the actual name of
  the input series with the number of years and seasons.}
  \item{alternative}{a character string describing
  alternative to the test ("two.sided").}
  \item{null.value}{the value for the hypothesized slope
  (0).}
}
\description{
  Compute the regional seasonal Kendall trend test.
}
\note{
  The value of \code{p.value} is \code{p.value.raw} if
  there are fewer than 10 years of data and is
  \code{p.value.corrected} otherwise.
}
\examples{
\dontrun{
data(RK3b)
# Build matrix, site 11 is missing 2004, would be row 60
Ammonia <- with(RK3b, matrix(c(value[1:59], NA, value[60:299]),
 ncol=25))
regionalSeaken(Ammonia, 1)
}
}

