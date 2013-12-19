\name{RK3b}
\alias{RK3b}
\docType{data}
\title{
Annual Series Data
}
\description{
  Ammonium (µeq/L) in snowpack, Colorado and New Mexico.
}
\usage{RK3b}
\format{
Data frame with 299 rows and 3 columns\cr
\tabular{lll}{
 Name \tab Type \tab Description\cr
 \cr
 year \tab integer \tab The year in which the sample was taken\cr
 site \tab integer \tab The site number\cr
 value \tab numeric \tab The ammonium concentration of the sample in µeq/L\cr
 }
}
\source{
Data retrieved from the supplemental data set for Helsel and others (
}
\references{
Helsel, D.R., Mueller, D.K., and Slack, J.R., 2006, 
Computer program for the Kendall family of trend tests: 
U.S. Geological Survey Scientific Investigations Report 2005–5275, 4 p.
}
\examples{
\dontrun{
data(RK3b)
# How many samples at each site?
with(RK3b, table(site))
}}
\keyword{datasets}
