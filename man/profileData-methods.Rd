\name{profileData-methods}
\docType{methods}
\alias{profileData-methods}
\alias{profileData}
\alias{profileData,StackedData-method}
\alias{profileData<--methods}
\alias{profileData<-}
\alias{profileData<-,StackedData-method}
\title{Method profileData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which retrieves the profileData
   from a StackedData object or which assigns the profileData data slot of a StackedData object.  
}
\section{Methods}{

 \describe{
    \item{\code{signature(object = "StackedData")}}{Get the profileData from the object.}
    \item{\code{signature(object = "StackedData", value = "data.frame")}}{Set the value of the profileData slot in the object.}
 }
}

 \usage{
   # Get the profileData from a StackedData object.
   profileData(object)    

   # Set the profileData slot of a StackedData object.
   profileData(object) = value
 }

 \arguments{
   \item{object}{an object of the StackedData class}
 }

 \value{
   # Get the profileData from a StackedData object.
   data frame of profile data.
 }

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{profileData}} }

\examples{

# Load the profile data, create a stackedData object and set the profile data slot
data(profileDF)
# Create a stacked data object
stackedDataObject = new("StackedData")
# Set the profile data slot
profileData(stackedDataObject) = profileDF

# Get the profile data from the stacked data object
profileData = profileData(stackedDataObject)

}