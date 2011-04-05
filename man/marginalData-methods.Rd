\name{marginalData-methods}
\docType{methods}
\alias{marginalData-methods}
\alias{marginalData<--methods}
\alias{marginalData}
\alias{marginalData<-}
\alias{marginalData,StackedData-method}
\alias{marginalData<-,StackedData-method}
\title{Method marginalData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which can retrieve the marginalData
   from a StackedData object or which can assign the marginalData data slot of a StackedData object.
}
\section{Methods}{

 \describe{
    \item{\code{signature(object="StackedData")}}{Get the value of the marginalData slot in the stackedDataObject.} 
    \item{\code{signature(object = "StackedData", value = "data.frame")}}{Set the value of the marginalData slot in the stackedDataObject.}
 }

}

 \usage{
  # Get the marginal data from a StackedData object
   marginalData(object)     

  # Set the marginal data slot of a StackedData object
  # marginalData(object) <- value
  \S4method{marginalData}{StackedData}(object) <- value 

 }

 \arguments{
   \item{object}{an object of the StackedData class}
   \item{value}{a replacement value} 
 }

 \value{
   data frame of marginal data.
 }


\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{marginalData}}}

\examples{
# Load the marginal data and set the marginal data slot
data(marginalDF)
# Create a stacked data object
stackedDataObject = new("StackedData")
# Set the marginal data slot
marginalData(stackedDataObject) = marginalDF

# Get the marginal data from the stacked data object
marginalData = marginalData(stackedDataObject)

}