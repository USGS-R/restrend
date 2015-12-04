### R code from vignette source 'Nutrients.Rnw'

###################################################
### code chunk number 1: Nutrients.Rnw:35-41
###################################################
# Load restrend and other packages and the data
library(restrend)
library(smwrBase)
library(smwrQW)
data(EstrendSub)
head(EstrendSub)


###################################################
### code chunk number 2: Nutrients.Rnw:49-60
###################################################
# Convert to class qw
EstrendSub.qw <- convert2qw(EstrendSub)
# Create the subset
Nuts <- subset(EstrendSub.qw, select=c("STAID", "DATES",
  "PN.organic", "PAmmonia", "PKjeldahl", "PTotal.P"))
# Rename to remove leading P, not required--just pretty
constituents <- c("N.organic", "Ammonia", "Kjeldahl", "Total.P")
names(Nuts)[3:6] <- constituents
# The sampling for nutrients started later, so remove the samples that
# have no nutrient data.
Nuts <- dropAllMissing(Nuts, constituents)


###################################################
### code chunk number 3: Nutrients.Rnw:65-67
###################################################
# Create the report
sampReport(Nuts, DATES="DATES", STAID="STAID", file="NutrientSampling")


###################################################
### code chunk number 4: Nutrients.Rnw:91-98
###################################################
# subset to a few selected stations:
Nuts <- subset(Nuts, STAID %in% c("07227500", "07228000", "07336820", "07343200",
    "07346070"))
# Set up the project
setProj("nutrients", Nuts, STAID="STAID", DATES="DATES", 
				Snames=constituents, 
				type="monthly", Start="1969-10-01", End="1989-10-01")


###################################################
### code chunk number 5: Nutrients.Rnw:118-122
###################################################
# Which are OK?
estrend.st
# What seasonal definition?
estrend.ms


###################################################
### code chunk number 6: Nutrients.Rnw:132-134
###################################################
# Trend tests, accepting default seasons
SKTrends()


###################################################
### code chunk number 7: Nutrients.Rnw:142-145
###################################################
# get the trends
nutrients.tnd <- getTrends()
print(nutrients.tnd)


###################################################
### code chunk number 8: Nutrients.Rnw:153-155
###################################################
# get the history
estrend.cl


