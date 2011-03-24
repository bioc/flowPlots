\name{computePFDPartsData-methods}
\docType{methods}
\alias{computePFDPartsData-methods}
\alias{computePFDPartsData}
\alias{computePFDPartsData,StackedData,character,character,character,character-method}
\title{Method computePFDPartsData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which computes the pfdParts data 
   which can be stored in the pfdPartsData slot of a StackedData object.  This method relies on the
   marker data slot being assigned in the StackedData object.  
}
\section{Methods}{
 \describe{
    \item{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}{Compute the pfd parts data.}
 }
}

\usage{
   computePFDPartsData(object, byVarNames, idVarName, percentVarName, groupVarName)    
}
\arguments{
   \item{object}{an object of the StackedData class}
   \item{byVarNames}{character; the names of the variables specifying the subsets of interest in the data}
   \item{idVarName}{character; the name of the id variable in the data}
   \item{percentVarName}{character; the name of the variable holding the percentages to be summed when computing
        the pfd summary data}
   \item{groupVarName}{character; the name of the variable specifying the group assignment in the data}
}

\value{
   a list. Each element of the list is a data frame containing the component percents for
   a given degree of polyfunctionality (PFD), except for the max PFD since there is only 
   one possible combination for the max PFD.  
}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{pfdPartsData}} }

\examples{
# Load stacked data
data(adultsNeonates)
# Create a stacked data object
stackedDataObject = new("StackedData", stackedData=adultsNeonates)

# Compute the marker data and set the marker data slot
markerNames = c("TNFa","IL6","IL12","IFNa")
markers = computeMarkers(markerNames,includeAllNegativeRow=TRUE)
markers(stackedDataObject) = markers

# Compute the pfd parts data and set the pfd parts data slot
byVarNames = c("stim", "concGroup", "cell")
pfdPartsData = computePFDPartsData(stackedDataObject, byVarNames, "id", "percentAll", "group")
pfdPartsData(stackedDataObject) = pfdPartsData

}