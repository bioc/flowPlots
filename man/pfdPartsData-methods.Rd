\name{pfdPartsData-methods}
\docType{methods}
\alias{pfdPartsData-methods}
\alias{pfdPartsData}
\alias{pfdPartsData,StackedData-method}
\alias{pfdPartsData<--methods}
\alias{pfdPartsData<-}
\alias{pfdPartsData<-,StackedData-method}
\title{Method pfdPartsData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which retrieves the pfdPartsData
   from a StackedData object or which assigns the pfdPartsData data slot of a StackedData object.    
}
\section{Methods}{ 
 \describe{
   \item{\code{signature(object = "StackedData")}}{Get the pfdPartsData from the object.}
   \item{\code{signature(object = "StackedData", value = "list")}}{Set the value of the pfdPartsData slot in the object.}
 }
}

 \usage{
   # Get the pfdPartsData from a Stacked Data object.
   pfdPartsData(object)    
 
   # Set the pfdPartsData of a Stacked Data object.
   # pfdPartsData(object) = value   
   \S4method{pfdPartsData}{StackedData}(object) <- value
 }

 \arguments{
   \item{object}{an object of the StackedData class}
   \item{value}{a replacement value}
 }

 \value{
   a list. Each element of the list is a data frame containing the component percents for
   a given degree of polyfunctionality (PFD), except for the max PFD since there is only 
   one possible combination for the max PFD.  
 }


\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{pfdPartsData}} }

\examples{

# Load the pfdParts data and set the pfdParts data slot
data(pfdPartsList)

# Create a stacked data object
stackedDataObject = new("StackedData")

# Set the pfdParts data slot
pfdPartsData(stackedDataObject) = pfdPartsList

# Get the pfd parts data
pfdPartsData = pfdPartsData(stackedDataObject)

}