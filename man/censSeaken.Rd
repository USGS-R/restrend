\name{censSeaken}
\alias{censSeaken}
\title{Trend Test}
\usage{
  censSeaken(series, nseas = 12)
}
\arguments{
  \item{series}{any regularly spaced object that can be
  forced to class "lcens" to test for trend. Missing values
  are permitted.}

  \item{nseas}{the number of seasons per year.}
}
\value{
  An object of class "htest" containing the following
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
  Compute the seasonal Kendall trend test with the Turnbull
  slope estimator for left-censored data.
}
\note{
  The value of \code{p.value} is \code{p.value.raw} if
  there are fewer than 10 years of data and is
  \code{p.value.corrected} otherwise.
}
\examples{
\dontrun{
# Compare censored and uncensored to seaken
library(USGSwsData)
data(KlamathTP)
# Construct the regular series
KlamathTP.rs <- with(KlamathTP, regularSeries(TP_ss, sample_dt,
 begin="1972-01-01", end="1980-01-01"))
# Uncensored, differences due to rounding
with(KlamathTP.rs, seaken(Value, 12))
with(KlamathTP.rs, censSeaken(Value, 12))
# About 30 percent censoring, censSeaken closer to uncensored slope
with(KlamathTP.rs, seaken(ifelse(Value < 0.05, 0.025, Value), 12))
with(KlamathTP.rs, censSeaken(as.lcens(Value, 0.05), 12))
}
}

