\name{legendPFDStatsGroupNames}
\alias{legendPFDStatsGroupNames}
\title{Create Group Names With Embedded Stats to Use in a Plot Legend}
\description{
   This function returns a vector of groupNames of the form:
   "Adults (25) 1.5/.6", which represents the group name, number of subjects in
   the group, the pfd mean / pfd standard deviation, where pfd = polyfunctional degree.  
}
\usage{
   legendPFDStatsGroupNames(pfdGroupStatsList,groupNames)   
}
\arguments{
   \item{pfdGroupStatsList}{a list of vectors containing the pfd group stats of 
     group size, pfd group mean, and pfd group standard deviation}
   \item{groupNames}{a vector of group names, such as c("Adult", "Neonate")}
}

\value{character vector; the length is equal to the number of groups.  Each element of the
vector contains the formatted stats for a group.  For a subset of the adultsNeonates data,
the 2-element vector result looks like this: "Adults (24), 1.47/0.10"   "Neonates (30), 1.41/0.13"}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{
## Load the PFD data and prepare it for a call to ternaryplot()

data(pfdDF)
pfdDataSubset = subset(pfdDF, stim=="LPS" & concGroup==3 & cell=="mDC")
ternaryData = makeTernaryData(pfdDataSubset, 1, 2, 3:4)
colnames(ternaryData) = c("PFD1", "PFD2", "PFD3-4")

## Make a ternary plot

library(vcd)
ternaryplot(ternaryData, cex=.5, col=as.numeric(pfdDataSubset$group)*2, main="Stimulation = LPS, 
   Concentration Group = 3, Cell = mDC")

## Compute PFD Stats to print in plot legend

adultPFDData = subset(pfdDataSubset, group=="adult", select=c(PFD1:PFD3))
neoPFDData = subset(pfdDataSubset, group=="neonate", select=c(PFD1:PFD3))
groupPFDDataList = list(adultPFDData, neoPFDData)
pfdGroupStatsList = computePFDGroupStatsList(groupPFDDataList, pfdValues=1:3, numDigitsMean=3, 
   numDigitsSD=2)
groupNames = c("Adults","Neonates")

## Create group names including the PFD Stats for the legend

legendNames = legendPFDStatsGroupNames(pfdGroupStatsList,groupNames)


## Add the legend to the ternary plot

grid_legend(0.8, 0.7, pch=c(20,20), col=c(2,4), legendNames, title = "Group (n), mean/sd:", 
   gp=gpar(cex=.8))
}
