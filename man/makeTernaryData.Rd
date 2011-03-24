\name{makeTernaryData}
\alias{makeTernaryData}
\title{Prepare PFD Data for a Call to ternaryplot()}
\description{
   This function takes a dataframe of polyfunctional degree (pfd) data and 
   prepares a matrix with 3 columns (egs. PFD1,PFD2,PFD3, or PFD1-2, PFD3-4, PFD5-6) 
   to use as input to ternaryplot() in the vcd pkg.  If multiple columns are specified
   for one of the three columns, the percentages for the multiple columns will be summed.  
}
\usage{
   makeTernaryData(pfdData, columns1, columns2, columns3, 
      columnNames=c("PFD=1","PFD=2","PFD=3"))
}
\arguments{
   \item{pfdData}{data frame containing data for 1 or more groups}
   \item{columns1}{numeric; column(s) of pfd data to place in the first column of the matrix}
   \item{columns2}{numeric; column(s) of pfd data to place in the second column of the matrix}
   \item{columns3}{numeric; column(s) of pfd data to place in the third column of the matrix}
   \item{columnNames}{(optional) character vector of names for the three columns of the matrix}
}

\value{matrix of 3 columns; rows = subjects, cols = data for each point on the triangle}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{
## Load the PFD Data to plot in a ternary plot

data(pfdDF)
pfdDataSubset = subset(pfdDF, stim=="LPS" & concGroup==3 & cell=="mDC")

## Prepare the PFD Data for a call to ternaryplot()

ternaryData = makeTernaryData(pfdDataSubset, 1, 2, 3:4)
colnames(ternaryData) = c("PFD1", "PFD2", "PFD3-4")

## Make the ternary plot

library(vcd)
ternaryplot(ternaryData, cex=.5, col=as.numeric(pfdDataSubset$group)*2, main="Stimulation = LPS, 
   Concentration Group = 3, Cell = mDC")

## Prepare a legend with group stats

adultPFDData = subset(pfdDataSubset, group=="adult", select=c(PFD1:PFD3))
neoPFDData = subset(pfdDataSubset, group=="neonate", select=c(PFD1:PFD3))
groupPFDDataList = list(adultPFDData, neoPFDData)
pfdGroupStatsList = computePFDGroupStatsList(groupPFDDataList, pfdValues=1:3, numDigitsMean=3, 
   numDigitsSD=2)
groupNames = c("Adults","Neonates")
legendNames = legendPFDStatsGroupNames(pfdGroupStatsList,groupNames)
grid_legend(0.8, 0.7, pch=c(20,20), col=c(2,4), legendNames, title = "Group (n), mean/sd:", 
   gp=gpar(cex=.8))
}

