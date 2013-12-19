\name{flowAdjust}
\alias{flowAdjust}
\title{Flow Adjusted Concentrations}
\usage{
  flowAdjust(Stations = "All", Snames = "All",
    use.logs = TRUE, span = 0.75, max.cens = 5, report)
}
\arguments{
  \item{Stations}{a vector of the the station identifiers
  on which to do the flow adjustment.}

  \item{Snames}{a vector of the response variables on which
  to do the flow adjustment.}

  \item{use.logs}{fit a log-log LOWESS curve? If
  \code{FALSE}, then do not use the log transforms.}

  \item{span}{the span to use for the LOWESS curve.}

  \item{max.cens}{the maximum percent censoring permitted.
  A warning is printed if there are any censored data and a
  warning is printed if any exceed the value.}

  \item{report}{the name of the PDF file that contains
  graphs of each fit. The default is to use the name of the
  project with "_fa" appended. If the PDF file exists, then
  it is not overwritten, but the name is appended with a
  sequence of numbers until one that is valid is created.}
}
\value{
  The name of the report file.
}
\description{
  Create flow adjusted concentrations for the seasonal
  Kendall trend test.
}

