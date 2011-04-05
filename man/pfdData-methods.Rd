\name{pfdData-methods}
\docType{methods}
\alias{pfdData-methods}
\alias{pfdData}
\alias{pfdData,StackedData-method}
\alias{pfdData<--methods}
\alias{pfdData<-}
\alias{pfdData<-,StackedData-method}
\title{Method pfdData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which retrieves the pfdData
   from a StackedData object or which assigns the pfdData data slot of a StackedData object.  
}
\section{Methods}{ 

 \describe{
   \item{\code{signature(object = "StackedData")}}{Get the pfdData from the object. }
   \item{\code{signature(object = "StackedData", value = "data.frame")}}{Set the value of the pfdData slot in the object.}
  }
}

 \usage{
    # Get the pfdData from a StackedData object.
    pfdData(object)   

    # Set the pfdData slot of a StackedData object.
    # pfdData(object) = value
    \S4method{pfdData}{StackedData}(object) <- value
 }

 \arguments{
   \item{object}{an object of the StackedData class}
   \item{value}{a replacement value}
 }

 \value{
   data frame of pfd data.
 }


\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{pfdData}} }

\examples{
# Load the pfd data and set the pfd data slot
data(pfdDF)

# Create a stacked data object
stackedDataObject = new("StackedData")

# Set the pfd data slot
pfdData(stackedDataObject) = pfdDF

# Get the pfd data from the stacked data object
pfdData = pfdData(stackedDataObject)

}
