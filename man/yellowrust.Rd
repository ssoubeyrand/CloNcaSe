\name{yellowrust}
\alias{yellowrust}
\docType{data}
\title{
Yellow rust data set
}
\description{
Counts of strains of yellow rust disease collected at two different sampling dates.
}
\usage{data("yellowrust")}
\format{
  The format is:
 num [1:436, 1:2] 18 0 1 7 6 0 3 3 1 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "counts1" "counts2"
}
\details{
\code{yellowrust} is a 2-column matrix whose 1st column gives the counts of strains of yellow rust disease collected at 1st sampling date and 2nd column gives the counts of strains of yellow rust disease collected at 2nd sampling date.
}
\source{
	Ali S, Soubeyrand S, Gladieux P, Giraud T, Leconte M, Gautier A, Mboup M, de Vallavieille-Pope C and Enjalbert J (2015). CloNcaSe: Estimation of sex frequency and effective population size by clonemate re-sampling in partially clonal organisms. Research Report.
}
\references{
	Ali S, Soubeyrand S, Gladieux P, Giraud T, Leconte M, Gautier A, Mboup M, de Vallavieille-Pope C and Enjalbert J (2015). CloNcaSe: Estimation of sex frequency and effective population size by clonemate re-sampling in partially clonal organisms. Research Report.
}
\examples{
data(yellowrust)
yellowrust
}
\keyword{datasets}
