\name{undoFA}
\alias{undoFA}
\title{Flow Adjusted Concentrations}
\usage{
  undoFA(Stations, Snames)
}
\arguments{
  \item{Stations}{a vector of the the station identifiers
  on which to do the flow adjustment.}

  \item{Snames}{a vector of the response variables on which
  to do the flow adjustment.}
}
\value{
  Nothing is returned.
}
\description{
  Remove flow adjusted concentrations for the seasonal
  Kendall trend test.
}
\details{
  The arguments for \code{Stations} and \code{Snames} must
  be valid station identifiers and response varaibles.
  "All" is not valid.
}
\note{
  The \code{undoFA} function is used to remove the flow
  adjusted concentration data and prevent that analysis for
  the uncensored seasonal Kendall test. It is intended for
  the rare case when a reasonable flow-adjustment model
  cannot be found.
}

