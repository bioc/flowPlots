\name{profileDF}
\alias{profileDF}
\title{An example of profileData}
\description{
   The profile data computed from the adultsNeonates data.
}
\usage{
   data(profileDF) 
}

\value{a data frame containing cols for the percentages for each of the marker combinations (the 
profile percentages), as well as 'demographic' data describing each row of percentages.  For example, 
if there are 4 markers in the data set, there will be 2^4 columns for the profile percentages if all
possible combinations are included.  If the all-negative combination is excluded, there will be (2^4-1)
cols of profile percentages.}

\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}

\seealso{\code{\link{StackedData}} }

\examples{
# Load profile data
data(profileDF)
}
