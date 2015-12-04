### R code from vignette source 'Regional.Rnw'

###################################################
### code chunk number 1: Regional.Rnw:33-39
###################################################
# Load restrend and other packages and the data
library(restrend)
library(smwrBase)
library(smwrQW)
data(EstrendSub)
head(EstrendSub)


###################################################
### code chunk number 2: Regional.Rnw:47-57
###################################################
# Convert to class qw
EstrendSub.qw <- convert2qw(EstrendSub)
# Create the subset
KN <- subset(EstrendSub.qw, select=c("STAID", "DATES",
  "PKjeldahl"))
# Rename to remove leading P, not required--just pretty
names(KN)[3] <- "Kjeldahl"
# The sampling for nutrients started later, so remove the samples that
# have no data.
KN <- dropAllMissing(KN, "Kjeldahl")


###################################################
### code chunk number 3: Regional.Rnw:62-64
###################################################
# Create the report
sampReport(KN, DATES="DATES", STAID="STAID", file="KNSampling")


###################################################
### code chunk number 4: Regional.Rnw:80-84
###################################################
# Set up the project
setProj("kn", KN, STAID="STAID", DATES="DATES", 
				Snames="Kjeldahl", 
				type="seasonal", Start="1974-10-01", End="1989-10-01")


###################################################
### code chunk number 5: Regional.Rnw:89-95
###################################################
# Which are OK?
estrend.st
# What seasonal definition?
estrend.ss
# What about censoring percentage?
estrend.cp


###################################################
### code chunk number 6: Regional.Rnw:105-111
###################################################
# Trend tests, accepting default seasons (6) and percent censoring (5)
SKTrends()
# repeat for the more highly censored stations
# The use of the log transform must be turned off to make it comparable
SKTrends(Stations="07297910", max.cens=10, use.logs=FALSE,
  report="kn_07297910")


###################################################
### code chunk number 7: Regional.Rnw:118-121
###################################################
# The log transform is turned on to make it comparable to the other stations
SKTrends(Stations="07297910", max.cens=10, use.logs=TRUE,
  report="kn_07297910_log")


###################################################
### code chunk number 8: Regional.Rnw:128-131
###################################################
# get the trends
kn.tnd <- getTrends()
print(kn.tnd)


###################################################
### code chunk number 9: Regional.Rnw:136-138
###################################################
# Any regional trend?
RSKTrends(kn.tnd$Station, "Kjeldahl")


###################################################
### code chunk number 10: Regional.Rnw:147-149
###################################################
# get the history
estrend.cl


