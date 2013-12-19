\name{sampReport}
\alias{sampReport}
\title{Summarize Samples}
\usage{
  sampReport(data, DATES = "DATES", STAID = "STAID", file)
}
\arguments{
  \item{data}{the dataset to summarize}

  \item{DATES}{the name of the column containing the sample
  dates}

  \item{STAID}{the name of the column containing the
  station identifiers}

  \item{file}{the output file base name; the .pdf suffix is
  appended to make the actual file name. If missing, then
  the name of \code{data} is used as the base name.}
}
\value{
  The actual file name is returned invisibly.
}
\description{
  Create a multi-page pdf file of sample data by station in
  a dataset. The first pages are a listing of the first and
  last sample and total numerof of samples at each station.
  The following pages are dot plots of the sample dates. No
  more than 40 stations per page are liested or plotted.
}

