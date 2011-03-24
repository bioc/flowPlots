\name{makeBarplotData}
\alias{makeBarplotData}
\title{Prepare Profile Data for a Call to Barplot()}
\description{
   This function takes a dataframe of profile data and prepares a matrix of data for 
   input to barplot()
}
\usage{
   makeBarplotData(profileData, profileColumns, groupVariableName)
}
\arguments{
   \item{profileData}{dataframe of profile data, such as cell %'s in different 
     categories of cells}
   \item{profileColumns}{the columns of profileData to include in the barplot}
   \item{groupVariableName}{the column in the dataframe containing the group info}
}

\value{
   a matrix whose rows represent different profile categories and whose
   columns represent different groups.  Each cell in the matrix contains the 
   mean value for the group for a given profile category.  For profile data
   with 16 cytokine combinations for each of 2 groups, the matrix returned with have 
   dimensions (16,2). 
}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{GroupListBoxplot}} }

\examples{
data(profileDF)
profileDataSubset = subset(profileDF, stim=="LPS" & concGroup==3 & cell=="mDC")
profileColumns = 1:16
barplotData = makeBarplotData(profileDataSubset, profileColumns, groupVariableName="group")
barplotDataWithLegend = cbind(barplotData, NA, NA)
barColors = gray(0:15/15)[16:1]
barplot(barplotDataWithLegend, col=barColors, main="Stimulation = LPS 
   Concentration Group = 3 
   Cell = mDC")
legendNames = rownames(barplotData)
legend(2.75, 100, legend=legendNames[16:1], col=barColors[16:1], cex=.8, pch=20)
}
