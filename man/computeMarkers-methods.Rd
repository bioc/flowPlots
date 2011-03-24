\name{computeMarkers-methods}
\docType{methods}
\alias{computeMarkers-methods}
\alias{computeMarkers}
\alias{computeMarkers,character,logical-method}
\title{Method computeMarkers from Class "StackedData"}
\description{
   This function is a method of the StackedData class which computes the markers which 
   can be stored in the markers data slot of a StackedData object.  The marker matrix should 
   match the order of the rows in the stacked data file.  The stacked data should be 
   sorted so that the order of each 'stack' in the file is the same.  This method can 
   be used to compute a marker matrix.  If that matrix does not match the order of the 'stack', 
   then the user can generate the marker matrix separately and assign it to the marker the data slot
   in a StackedData object. 
}
\section{Methods}{
 \describe{
     \item{\code{signature(markerNames = "character", includeAllNegativeRow = "logical")}}{Compute the markers. Include the all negative row if the data includes the all negative case, such as: TNFa-IFNg-IL2-, in this case with 3 markers.}
 }
}

\usage{
   computeMarkers(markerNames, includeAllNegativeRow)    
}
\arguments{
   \item{markerNames}{character; vector of the names of the markers}
   \item{includeAllNegativeRow}{logical; TRUE, if the stacked data contains the all-negative row of markers; 
          for example, TNFa-IL6-IL12-IFNa- }
}

\value{matrix of 0's and 1's; rows represent marker combinations, cols represent markers.}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{markers}} }

\examples{
# Compute the marker data
markerNames = c("TNFa","IL6","IL12","IFNa")
markers = computeMarkers(markerNames,includeAllNegativeRow=TRUE)

## If you're using a StackedData object to compute summary data

# Create a stacked data object
stackedDataObject = new("StackedData")

# Assign the markers to the marker data slot
markers(stackedDataObject) = markers
}
