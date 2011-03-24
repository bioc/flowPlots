\name{computePFDData-methods}
\docType{methods}
\alias{computePFDData-methods}
\alias{computePFDData}
\alias{computePFDData,StackedData,character,character,character,character-method}
\title{Method computePFDData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which computes the polyfunctional degree (pfd) data 
   (usually recorded as percentages of reactive cells) which can be stored in the pfdData slot of a StackedData object.  
   PFD=1 refers to cells which are producing only one marker.  PFD=2 refers to cells which are producing exactly two 
   markers.  Similarly, up to PFD=n, where n is the number of markers. This method relies on the marker data slot 
   being assigned in the StackedData object.
}
\section{Methods}{
 \describe{
    \item{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}{Computes the pfd data.}
 }
}

\usage{
   computePFDData(object, byVarNames, idVarName, percentVarName, groupVarName)    
}
\arguments{
   \item{object}{an object of the StackedData class}
   \item{byVarNames}{character; the names of the variables specifying the subsets of interest in the data}
   \item{idVarName}{character; the name of the id variable in the data}
   \item{percentVarName}{character; the name of the variable holding the percentages to be summed when computing
        the pfd summary data}
   \item{groupVarName}{character; the name of the variable specifying the group assignment in the data}
}

\value{data frame of pfd data}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{pfdData}} }

\examples{
# Load stacked data
data(adultsNeonates)
# Create a stacked data object
stackedDataObject = new("StackedData", stackedData=adultsNeonates)

# Compute the marker data and set the marker data slot
markerNames = c("TNFa","IL6","IL12","IFNa")
markers = computeMarkers(markerNames,includeAllNegativeRow=TRUE)
markers(stackedDataObject) = markers

# Compute the pfd data and set the pfd data slot
byVarNames = c("stim", "concGroup", "cell")  
pfdData = computePFDData(stackedDataObject, byVarNames, "id", "percentAll", "group")
pfdData(stackedDataObject) = pfdData
 
}
