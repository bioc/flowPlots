\name{pfdDF}
\alias{pfdDF}
\title{An example of pfdData}
\description{
   The polyfunctional degree (pfd) data computed from the adultsNeonates data.   
}
\usage{
   data(pfdDF) 
}

\value{data frame; The data frame contains n columns of PFD percentages of reactive cells, where n is the number of 
markers in the data set, as well as 'demographic' data for each row of percentages.  PFD=1 refers to
the percentages of reactive cells which are producing only one marker, irregardless of which marker.  PFD=2 refers 
to cells which are producing exactly two markers.  Similarly, up to PFD=n markers.}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{
# Load pfd data
data(pfdDF)
}
