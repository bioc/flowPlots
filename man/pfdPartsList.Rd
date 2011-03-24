\name{pfdPartsList}
\alias{pfdPartsList}
\title{An example of pfdPartsData}
\description{
   The pfd parts data computed from the adultsNeonates data.  
}
\usage{
   data(pfdPartsList) 
}

\value{a list; pdfPartsData is a list of data frames.
   The first data frame holds the compositional percentages for PFD=1; that is, of the cells producing only 
   one marker, the percentage of cells which express marker1, the percentage of cells which produce marker2, etc.
   For example, if there are 4 markers, the data frame for PFD=1 will have a percentage column for each marker, and 
   the sum of those 4 cols will equal 100%.  The second data frame holds the compositional percentages for PFD=2.  
   The second data frame has n-choose-2 columns of percentages, where n is the number of markers and 2 is the 
   polyfunctional degree; i.e. the 2 in PFD=2.  And, so on.  The length of the list is (maxPFD-1), since 
   there is only 1 way to achieve maxPFD; that is, all markers are positive when PFD=maxPFD. Each data frame 
   also contains 'demographic' data describing each row of percentages.
}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{
# Load pfd parts data
data(pfdPartsList)
}
