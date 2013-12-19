\name{plotTrends}
\alias{plotTrends}
\title{Plot Trends}
\usage{
  plotTrends(data, which = "Trend", device)
}
\arguments{
  \item{data}{the data set returned from \code{getTrends}.}

  \item{which}{which trend to plot, either "Trend" or
  "Trend.pct."}

  \item{device}{the name of the graphics device. Defaults
  to a valid graphics device on any platform. May be "pdf"
  to create a pdf file.}
}
\value{
  The name of the graphics device.
}
\description{
  Create a series of diagnostic plots for a single Tobit
  trend test.
}
\note{
  The project from which \code{data} was retrieved must be
  the current project.
}

