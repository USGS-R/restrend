### R code from vignette source 'MajorIons.Rnw'

###################################################
### code chunk number 1: MajorIons.Rnw:32-37
###################################################
# Load the restrend and smwrBase packages and the data
library(restrend)
library(smwrBase)
data(EstrendSub)
head(EstrendSub)


###################################################
### code chunk number 2: MajorIons.Rnw:45-50
###################################################
# Compute FLOW, the coalesce function is in smwrBase
EstrendSub <- transform(EstrendSub, FLOW=coalesce(QI, QD))
# Create the subset
Majors <- subset(EstrendSub, select=c("STAID", "DATES", "FLOW", 
    "Calcium", "Chloride"))


###################################################
### code chunk number 3: MajorIons.Rnw:55-57
###################################################
# Create the report
sampReport(Majors, DATES="DATES", STAID="STAID", file="MajorIonSampling")


###################################################
### code chunk number 4: MajorIons.Rnw:71-75
###################################################
# Set up the project
setProj("majors", Majors, STAID="STAID", DATES="DATES", 
        Snames=c("Calcium", "Chloride"), FLOW="FLOW", 
        type="seasonal", Start="1968-10-01", End="1989-10-01")


###################################################
### code chunk number 5: MajorIons.Rnw:95-99
###################################################
# Which are OK?
estrend.st
# What seasonal definition?
estrend.ss


###################################################
### code chunk number 6: MajorIons.Rnw:110-112
###################################################
# Do the flow adjustment accepting all defaults
flowAdjust()


###################################################
### code chunk number 7: MajorIons.Rnw:118-120
###################################################
# Do the flow adjustment accepting all defaults
flowAdjust(Station="07331600", Snames=c("Calcium", "Chloride"), span=1)


###################################################
### code chunk number 8: MajorIons.Rnw:130-132
###################################################
# Trend tests, accepting default seasons
SKTrends()


###################################################
### code chunk number 9: MajorIons.Rnw:140-143
###################################################
# get the trends
majors.tnd <- getTrends()
print(majors.tnd)


###################################################
### code chunk number 10: MajorIons.Rnw:151-153
###################################################
# get the history
estrend.cl


