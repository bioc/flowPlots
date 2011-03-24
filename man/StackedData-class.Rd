\name{StackedData-class}
\Rdversion{1.1}
\docType{class}
\alias{StackedData-class}
\alias{StackedData}
\title{Class "StackedData"}
\description{
   Stacked data refers to gated data originating from an ICS Flow Cytometry experiment where the marker
   combinations for a subset of interest, say a given cell type, stimulus, and concentration, are "stacked". 
   A common type of marker is a cytokine.  A subset of stacked data could look like this:

\preformatted{
   id   group stim concGroup cell percentAll count totalCount percentReactive            cytCombo

a2004 	adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6+IL12+IFNa+

a2004 	adult  LPS         3  mDC       0.43     3        700        0.940625 TNFa+IL6+IL12+IFNa-

a2004 	adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6+IL12-IFNa+

a2004 	adult  LPS         3  mDC      21.86   153        700       47.818750 TNFa+IL6+IL12-IFNa-

a2004 	adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6-IL12+IFNa+

a2004 	adult  LPS         3  mDC       0.29     2        700        0.634375 TNFa+IL6-IL12+IFNa-

a2004 	adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6-IL12-IFNa+

a2004 	adult  LPS         3  mDC      19.71   138        700       43.115625 TNFa+IL6-IL12-IFNa-
}

}

\section{Objects from the Class}{
Objects can be created by calls of the form: 
	
   \code{stackedDataObject = new("StackedData", stackedData=NA, profileData=NA, marginalData=NA, 
       pfdData=NA, pfdPartsData=NA, markers=NA )}
}
\section{Slots}{
  \describe{
    \item{\code{stackedData}:}{\code{"data.frame"} of stacked data}
    \item{\code{profileData}:}{\code{"data.frame"} of profile data}
    \item{\code{marginalData}:}{\code{"data.frame"} of marginal data}
    \item{\code{pfdData}:}{\code{"data.frame"} of pfd data}
    \item{\code{pfdPartsData}:}{\code{"list"} of \code{"data.frame"} 's of pfd parts data}
    \item{\code{markers}:}{\code{"matrix"} of marker data}
  }
}

\section{Methods}{
  \describe{
    \item{computeMarginalData}{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}
    \item{computeMarkers}{\code{signature(markerNames = "character", includeAllNegativeRow = "logical")}}
    \item{computePFDData}{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}
    \item{computePFDPartsData}{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}
    \item{computeProfileData}{\code{signature(object = "StackedData", byVarNames = "character", idVarName = "character", percentVarName = "character", groupVarName = "character")}}
    # getters
    \item{marginalData}{\code{signature(object = "StackedData")}}
    \item{markers}{\code{signature(object = "StackedData")}}
    \item{pfdData}{\code{signature(object = "StackedData")}}
    \item{pfdPartsData}{\code{signature(object = "StackedData")}}
    \item{profileData}{\code{signature(object = "StackedData")}}
    \item{stackedData}{\code{signature(object = "StackedData")}}
    # setters
    \item{marginalData}{\code{signature(object = "StackedData", value = "data.frame")}}
    \item{markers}{\code{signature(object = "StackedData", value = "matrix")}}
    \item{pfdData}{\code{signature(object = "StackedData", value = "data.frame")}}
    \item{pfdPartsData}{\code{signature(object = "StackedData", value = "list")}}
    \item{profileData}{\code{signature(object = "StackedData", value = "data.frame")}}
    \item{stackedData}{\code{signature(object = "StackedData", value = "data.frame")}}
  }
}


\details{
The marker combinations in the stacked data should be ordered within each subset of 
interest (for example: subjectID, celltype, concentration, and stimulation) to match
the marker matrix.  If the data are not in this order, it should be sorted into this
order before using the computeProfileData, computeMarginalData, computePFDData, 
computePFDPartsData methods.  The computeMarkers method can be used to compute the
marker matrix used by these methods.  If the matrix computed does not match the order
of your data, then you can supply your own marker matrix, assign it to the marker data
slot of a StackedData object, and then use the 'compute' methods to compute the other
types of data.
}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{markers}}, \code{\link{profileData}}, \code{\link{marginalData}}, \code{\link{pfdData}}, 
   \code{\link{pfdPartsData}},  \code{\link{adultsNeonates}} }

\examples{
# View the Data Slots in StackedData
showClass("StackedData")

# Load stacked data
data(adultsNeonates)
# Create a stacked data object
stackedDataObject = new("StackedData", stackedData=adultsNeonates)

# Compute the marker data and set the marker data slot
# The marker matrix computed here matches the order of the marker combinations
# in the adultsNeonates stacked data within each subset of interest
markerNames = c("TNFa","IL6","IL12","IFNa")
markers = computeMarkers(markerNames,includeAllNegativeRow=TRUE)
markers(stackedDataObject) = markers

# Compute the profile data and set the profile data slot
byVarNames = c("stim", "concGroup", "cell")
profileData = computeProfileData(stackedDataObject, byVarNames, "id", "percentAll", "group")
profileData(stackedDataObject) = profileData

# Compute the marginal data and set the marginal data slot
byVarNames = c("stim", "concGroup", "cell")
marginalData = computeMarginalData(stackedDataObject, byVarNames, "id", "percentAll", "group")
marginalData(stackedDataObject) = marginalData

# Compute the pfd data and set the pfd data slot
byVarNames = c("stim", "concGroup", "cell")  
pfdData = computePFDData(stackedDataObject, byVarNames, "id", "percentAll", "group")
pfdData(stackedDataObject) = pfdData

# Compute the pfd parts data and set the pfd parts data slot
byVarNames = c("stim", "concGroup", "cell")
pfdPartsData = computePFDPartsData(stackedDataObject, byVarNames, "id", "percentAll", "group")
pfdPartsData(stackedDataObject) = pfdPartsData

# Get the data from the stacked data object

markers = markers(stackedDataObject)
profileData = profileData(stackedDataObject)
marginalData = marginalData(stackedDataObject)
pfdData = pfdData(stackedDataObject)
pfdPartsData = pfdPartsData(stackedDataObject)

}


