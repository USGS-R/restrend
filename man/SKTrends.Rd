\name{SKTrends}
\alias{SKTrends}
\title{Seasonal Kendall Trends}
\usage{
  SKTrends(Stations = "All", Snames = "All",
    use.logs = TRUE, max.cens = 5, nseas = NULL, report)
}
\arguments{
  \item{Stations}{a vector of the the station identifiers
  on which to do the flow adjustment.}

  \item{Snames}{a vector of the response variables on which
  to do the flow adjustment.}

  \item{use.logs}{log transform the data before the trend
  test, applies only to uncensored seasonal Kendall test?}

  \item{max.cens}{the maximum percent censoring permitted
  for the uncensored seasonal Kendall test. If the
  percentage od censoring exceeds this value, then the
  censored seasonal Kendall test is performed. Set to a
  negative value to force the censored seasonal Kendall
  test for all \code{Stations} and \code{Snames}.}

  \item{nseas}{the number of seasons to use for all of the
  tests. If \code{NULL} (default), then use the selected
  number of seasons defined in \code{setProj}.}

  \item{report}{the name of the PDF file that contains a
  report for each test. The default is to use the name of
  the project with "_sk" appended. If the PDF file exists,
  then it is not overwritten, but the name is appended with
  a sequence of numbers until one that is valid is
  created.}
}
\value{
  The name of the report file.
}
\description{
  Perform the seasonal Kendall trend test.
}

