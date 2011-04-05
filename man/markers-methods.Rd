\name{markers-methods}
\docType{methods}
\alias{markers-methods}
\alias{markers}
\alias{markers,StackedData-method}
\alias{markers<--methods}
\alias{markers<-}
\alias{markers<-,StackedData-method}
\title{Method markers from Class "StackedData"}
\description{
   This function is a method of the StackedData class which retrieves the markers
   from a StackedData object or which assigns the markers data slot of a StackedData object. 
}
\section{Methods}{ 
 \describe{
    \item{\code{signature(object = "StackedData")}}{Get the marker matrix from the object.}
    \item{\code{signature(object = "StackedData", value = "matrix")}}{Set the value of the markers data slot in the object.}
 }
}
 \usage{
   # Get the marker matrix from a StackedData object.
   markers(object)    

   # Set the marker matrix slot of a StackedData object.
   \S4method{markers}{StackedData}(object) <- value
 }
 \arguments{
   \item{object}{an object of the StackedData class}
   \item{value}{a replacement value}
 }
 \value{
   matrix of markers data.
 }

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{markers}} }

\examples{

# Load the marker data and set the marker data slot
data(markerMatrix)

# Create a stacked data object
stackedDataObject = new("StackedData")

# Set the marker data slot
markers(stackedDataObject) = markerMatrix

# Get the marker data from the stacked data object
markers = markers(stackedDataObject)

}
