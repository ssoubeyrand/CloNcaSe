\name{simul}
\alias{simul}

\title{
Function simulating data sets
}
\description{
This function simulates the evolution and the sampling of a partially clonal population.
}
\usage{
simul(start.counts, param, nb1ClonalCycles, nb2ClonalCycles, nb1Sample, nb2Sample, 
	nbPrelCycles)
}

\arguments{
  \item{start.counts}{
Vectors of non-negative integers providing counts of strains used to initialize the simulations.
}
  \item{param}{
2-component numeric vector with the values of Nc and s used for the simulations.
}
  \item{nb1ClonalCycles}{
Non-negative integer indicating the number of clonal cycles separating the 1st sample and the mixed (i.e. clonal-sexual) cycle.
}
  \item{nb2ClonalCycles}{
Non-negative integer indicating the number of clonal cycles separating the mixed (i.e. clonal-sexual) cycle and the 2nd sample.
}
  \item{nb1Sample}{
Non-negative integer indicating the size of the 1st sample.
}
  \item{nb2Sample}{
Non-negative integer indicating the size of the 2nd sample.
}
  \item{nbPrelCycles}{
Number of preliminary series of cycles (one series consists of \code{nb1ClonalCycles} clonal cycles, 1 mixed cycle and \code{nb2ClonalCycles} clonal cycles) performed before the collection of the 1st sample.
}
}

\value{
A list consisted of:
\item{population }{Counts of strains in the population at the 1st and 2nd sampling dates.}
\item{samples }{Counts of sampled strains in the population at the 1st and 2nd sampling dates.}
\item{statistics }{Numeric vector providing the following statistics: the empirical probabilities of sampling two clonemates at the 1st and 2nd sampling dates (1st and 2nd components of  the vector), and the empirical probability that two random isolates, the 1st one being collected at the 1st sampling date and the 2nd one being collected at the 2nd sampling date, are clonemates (5th component of the vector). The 3rd and 4th components of the vector can be ignored.}
}
\references{
	Ali S, Soubeyrand S, Gladieux P, Giraud T, Leconte M, Gautier A, Mboup M, de Vallavieille-Pope C and Enjalbert J (2015). CloNcaSe: Estimation of sex frequency and effective population size by clonemate re-sampling in partially clonal organisms. Research Report.
}
\author{
Samuel Soubeyrand
}


\seealso{
\code{\link{ncase4simul}}.
}
\examples{
## count data used to initialize the simulation
counts1=c(20,10,5,1,1,1,1,1)

## simulated data set
simul(start.counts=counts1,param=c(1000,0.8),nb1ClonalCycles=17,nb2ClonalCycles=7,
	nb1Sample=100,nb2Sample=100,nbPrelCycles=5)
}

\keyword{ ~misc }
