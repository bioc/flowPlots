\name{makeDataList}
\alias{makeDataList}
\title{Prepare Data in a Data Frame for a Call to GroupListBoxplot()}
\description{
This function makes a list of data, where each item of the list is a data frame and 
contains data for a group.  The list can be used as input for a call to GroupListBoxplot().     
}
\usage{
    makeDataList(theData, groupVariableName, columnsToKeep)
}
\arguments{
   \item{theData}{data.frame containing data for 1 or more groups}
   \item{groupVariableName}{character; name of the column in theData identifying group}
   \item{columnsToKeep}{numeric vector specifying the column numbers in theData to include 
     in the list of data.}
}

\value{a list; each item of the list is a data frame and contains data for a group. The rows
  of each data frame are for subject.  The cols of each data frame represent categories to be
  plotted on the x-axis.}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{GroupListBoxplot}} }

\examples{

# Load the data
data(marginalDF)
marginalDataSubset = subset(marginalDF, stim=="LPS" & concGroup==3 & cell=="mDC")

# Make a data list
dataList = makeDataList(marginalDataSubset, "group", 1:5) 

# Make a plot using the data list
GroupListBoxplot(dataList, xlabel="Cytokine", ylabel="Percent of All Cells", 
   xAxisLabels=c("TNFa","IL6","IL12","IFNa","AnyMarker"), 
   mainTitle="Stimulation = LPS and Concentration Group = 3 and Cell = mDC", 
   legendGroupNames=c("Adults","Neonates"), pointColor=c(2,4), testTitleCEX=.8, nCEX=.8, 
   pCEX=.8, legendColor=c(2,4), legendCEX=.7)

}
