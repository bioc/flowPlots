\name{marginalDF}
\alias{marginalDF}
\title{An example of marginalData}
\description{
   The marginal data computed from the adultsNeonates data.
}
\usage{
   data(marginalDF) 
}

\value{data frame; There are (n+1) columns in the data frame for the marginal percentages; one
col for each marker and 1 col for "anyMarker", which is the sum of the individual marker cols.  The
data frame also has cols of 'demographic' data describing each row of percentages.}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{
# Load marginal data
data(marginalDF)
}
