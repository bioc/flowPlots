\name{computePFDGroupStatsList}
\alias{computePFDGroupStatsList}
\title{Compute Group Stats on PFD Data to Be Used In a Legend}
\description{
   This function can be used with ternaryplot() to add PFD group stats to, say, the legend.
   The stats computed are group size (N), pfd group mean, and pfd group standard deviation.
}
\usage{
   computePFDGroupStatsList(groupPFDDataList, pfdValues=1:3, numDigitsMean=3, numDigitsSD=2)    
}
\arguments{
   \item{groupPFDDataList}{one list item per group, each list item contains a matrix of PFD       
     percentages; the rows are subjects, and the columns are pfd categories.}
   \item{pfdValues}{vector of the PFD values that the columns in each matrix in the 
     groupPFDDataList represent; eg. 1:3 for (PFD1,PFD2,PFD3).}
   \item{numDigitsMean}{return a mean rounded to this number of digits}
   \item{numDigitsSD}{return a standard deviation rounded to this number of digits}
}

\value{a list; each list item contains the stats for a group as a 3 element character vector containing
the size of the group, the mean PFD, and the standard deviation of the PFD.}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{pfdData}} }

\examples{
## Load PFD data to plot
data(pfdDF)
pfdDataSubset = subset(pfdDF, stim=="LPS" & concGroup==3 & cell=="mDC")

## Prepare the PFD Data for a call to ternaryplot()
ternaryData = makeTernaryData(pfdDataSubset, 1, 2, 3:4)
colnames(ternaryData) = c("PFD1", "PFD2", "PFD3-4")

## Make a ternary plot
library(vcd)
ternaryplot(ternaryData, cex=.5, col=as.numeric(pfdDataSubset$group)*2, main="Stimulation = LPS, 
   Concentration Group = 3, Cell = mDC")

## Compute Group Stats to use in the legend of the ternary plot
adultPFDData = subset(pfdDataSubset, group=="adult", select=c(PFD1:PFD3))
neoPFDData = subset(pfdDataSubset, group=="neonate", select=c(PFD1:PFD3))
groupPFDDataList = list(adultPFDData, neoPFDData)

## Specifically, compute the PFD Group Stats List
pfdGroupStatsList = computePFDGroupStatsList(groupPFDDataList, pfdValues=1:3, numDigitsMean=3, 
   numDigitsSD=2)
groupNames = c("Adults","Neonates")

## Create group names for the legend based on the PFD Group Stats List
legendNames = legendPFDStatsGroupNames(pfdGroupStatsList,groupNames)
grid_legend(0.8, 0.7, pch=c(20,20), col=c(2,4), legendNames, title = "Group (n), mean/sd:", 
   gp=gpar(cex=.8))

}


