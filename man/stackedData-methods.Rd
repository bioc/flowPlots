\name{stackedData-methods}
\docType{methods}
\alias{stackedData-methods}
\alias{stackedData}
\alias{stackedData,StackedData-method}
\alias{stackedData<--methods}
\alias{stackedData<-}
\alias{stackedData<-,StackedData-method}
\title{Method stackedData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which retrieves the stackedData
   from a StackedData object or which assigns the stackedData data slot of a StackedData object. 
}
\section{Methods}{
 \describe{
     \item{\code{signature(object = "StackedData")}}{Get the stackedData from the object.}
     \item{\code{signature(object = "StackedData", value = "data.frame")}}{Set the value of the stackedData slot in the object.}
 }
}

\usage{
   # Get the stacked data from the StackedData object   
   stackedData(object)    

   # Set the stacked data slot of a StackedData object
   # stackedData(object) = value
   \S4method{stackedData}{StackedData}(object) <- value
}
\arguments{
   \item{object}{an object of the StackedData class}
   \item{value}{a replacement value}
}

\value{
   data frame of stacked data.
}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{

## Set the stacked data slot WHILE creating a new stacked data object

# Load stacked data
data(adultsNeonates)
# Create a stacked data object
stackedDataObject = new("StackedData", stackedData=adultsNeonates)

## Set the stacked data slot AFTER creating a new stacked data object

# Load stacked data
data(adultsNeonates)
# Create a stacked data object
stackedDataObject = new("StackedData")
# Set the stacked data slot
stackedData(stackedDataObject) = adultsNeonates

## Set the stacked data slot after creating a new stacked data object

stackedDataFrame = readStackedData(fileName=system.file("extdata","adultsNeonates.csv", package="flowPlots"))    
stackedDataObject = new("StackedData")
stackedData(stackedDataObject) = stackedDataFrame

}
