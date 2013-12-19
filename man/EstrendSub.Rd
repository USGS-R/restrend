\name{EstrendSub}
\alias{EstrendSub}
\docType{data}
\title{
EstrendSub Data
}
\description{
A subset of stations and water-quality constituents from Schertz and others (1991). 
The water-quality constituents were selected to cover a range of censoring levels. 
The stations were selected to represent a range of sampling intensity and duration.
}
\usage{EstrendSub}
\format{
Data frame with 5428 rows and 18 columns\cr
\tabular{lll}{
 Name \tab Type \tab Description\cr
 \cr
 STAID \tab character \tab USGS station identifier\cr
 DATES \tab Date \tab Sample date\cr
 QI \tab numeric \tab Instantaneous streamflow at time of sample\cr
 QD \tab numeric \tab daily mean streamflow for sample\cr
 RN.organic \tab character \tab Remark code for organic nitrogen concentration\cr
 PN.organic \tab numeric \tab Organic nitrogen concentration in mg/L\cr
 RAmmonia \tab character \tab Remark code for ammonia concentration\cr
 PAmmonia \tab numeric \tab Ammonia concentration in mg/L as N\cr
 RKjeldahl \tab character \tab Remark code for Kjeldahl nitrogen concentration\cr
 PKjeldahl \tab numeric \tab Kjeldahl nitrogen concentration in mg/L\cr
 RTotal.P \tab character \tab Remark code for total (whole water) phosphorus concentration\cr
 PTotal.P \tab numeric \tab Total (whole water) phosphorus concentrationin mg/L\cr
 RCopper \tab character \tab Remark code for copper concentration\cr
 PCopper \tab numeric \tab Copper concentration in ug/L\cr
 RIron \tab character \tab Remark code for iron concentration\cr
 PIron \tab numeric \tab Iron concentration in ug/L\cr
 Calcium \tab numeric \tab Calcium concentration in mg/L\cr
 Chloride \tab numeric \tab Chloride concentration in mg/L\cr
 }
}
\source{
The data include 8 water-quality constituents from 19 stations from Schertz and others (1991).
}
\references{
Schertz, T.L., Alexander, R.B., and Ohe, D.J., 1991, The computer program EStimate TREND 
(ESTREND), a system for the detection of trends in water-quality data: U.S. Geological 
Survey Water Resources Investigations Report 91-4040, 72 p. 
}
\examples{
\dontrun{
data(EstrendSub)
# Sampling date ranges for each station
with(EstrendSub, tapply(DATES, STAID, range))
}}
\keyword{datasets}
