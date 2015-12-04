### R code from vignette source 'Metals.Rnw'

###################################################
### code chunk number 1: Metals.Rnw:33-39
###################################################
# Load the restrend and other packages and the data
library(restrend)
library(smwrBase)
library(smwrQW)
data(EstrendSub)
head(EstrendSub)


###################################################
### code chunk number 2: Metals.Rnw:49-61
###################################################
# Compute FLOW, the coalese function is in smwrBase
EstrendSub <- transform(EstrendSub, FLOW=coalesce(QI, QD))
# Convert, the default scheme is "booker"
# The convert2qw function is in smwrQW
EstrendSub.qw <- convert2qw(EstrendSub)
# Create the subset, the Pcolumn name is preserved
Metals <- subset(EstrendSub.qw, select=c("STAID", "DATES", "FLOW", 
                                         "PIron", "PCopper"))
# Rename metals to remove the leading P
names(Metals)[4:5] <- c("Iron", "Copper")
# Show the first few rows of the data
head(Metals)


###################################################
### code chunk number 3: Metals.Rnw:66-71
###################################################
# Subset the data and show first few lines
Metals <- subset(Metals, !(is.na(Iron) & is.na(Copper)))
head(Metals)
# Create the report
sampReport(Metals, DATES="DATES", STAID="STAID", file="MetalsSampling")


###################################################
### code chunk number 4: Metals.Rnw:85-89
###################################################
# Set up the project
setProj("metals", Metals, STAID="STAID", DATES="DATES", 
        Snames=c("Iron", "Copper"), FLOW="FLOW", 
        type="tobit", Start="1974-10-01", End="1989-10-01")


###################################################
### code chunk number 5: Metals.Rnw:107-109
###################################################
# Which are OK?
estrend.st


###################################################
### code chunk number 6: Metals.Rnw:119-121
###################################################
# Trend tests, accepting default
tobitTrends()


###################################################
### code chunk number 7: Metals.Rnw:129-131
###################################################
# Trend tests, accepting default
plotTT(Station="07297910", Sname="Iron", device="pdf")


###################################################
### code chunk number 8: Metals.Rnw:139-142
###################################################
# get the trends
metals.tnd <- getTrends()
print(metals.tnd)


###################################################
### code chunk number 9: Metals.Rnw:150-152
###################################################
# get the history
estrend.cl


