\name{readStackedData-methods}
\docType{methods}
\alias{readStackedData}
\alias{readStackedData,character-method}
\title{Method readStackedData from Class "StackedData"}
\description{
   This function is a method of the StackedData class which reads a csv file of stacked data
   into a data frame, which can be stored in the stackedData data slot of a StackedData object. 
   This function is a wrapper for the R base function, read.csv(), for users not so familiar with R.
}
\section{Methods}{
 \describe{
    \item{\code{signature(fileName = "character")}}{Wrapper function for read.csv() to read a file of stacked data.}
 }
}

\usage{
   readStackedData(fileName)    
}
\arguments{
   \item{fileName}{character; the full name of the file containing stacked data in csv format.}
}

\value{
   data frame of stacked data.
}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}}, \code{\link{read.csv}} }

\examples{
 stackedDataFrame = readStackedData(fileName=system.file("extdata","adultsNeonates.csv", package="flowPlots"))   
}