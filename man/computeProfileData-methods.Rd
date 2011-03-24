\name{computeProfileData-methods}
\docType{methods}
\alias{computeProfileData-methods}
\alias{computeProfileData}
\alias{computeProfileData,StackedData,character,character,character,character-method}
\title{Method computeProfileData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which computes the profile data 
   which can be stored in the profileData slot of a StackedData object.  This method relies on the
   marker data slot being assigned in the StackedData object.  No computation is required.  Instead, the
   marker combination percentages data is re-organized to be 'horizontal' rather than 'stacked' (vertical). 
   This makes it ready for plotting via the GroupListBoxplot() function.    
}
\section{Methods}{
 \describe{
   \item{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}{Compute the profile data.}
 }
}

\usage{
   computeProfileData(object, byVarNames, idVarName, percentVarName, groupVarName)    
}
\arguments{
   \item{object}{an object of the StackedData class}
   \item{byVarNames}{character; the names of the variables specifying the subsets of interest in the data}
   \item{idVarName}{character; the name of the id variable in the data}
   \item{percentVarName}{character; the name of the variable holding the percentages to be summed when computing
        the pfd summary data}
   \item{groupVarName}{character; the name of the variable specifying the group assignment in the data}
}

\value{data frame of profile data}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{profileData}}  }

\examples{
# Load stacked data
data(adultsNeonates)
# Create a stacked data object
stackedDataObject = new("StackedData", stackedData=adultsNeonates)

# Compute the marker data and set the marker data slot
markerNames = c("TNFa","IL6","IL12","IFNa")
markers = computeMarkers(markerNames,includeAllNegativeRow=TRUE)
markers(stackedDataObject) = markers

# Compute the profile percent data and set the profile percent data slot
byVarNames = c("stim", "concGroup", "cell")
profilePercent = computeProfileData(stackedDataObject, byVarNames, "id", "percentAll", "group")
profileData(stackedDataObject) = profilePercent

}