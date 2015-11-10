\name{ncase4simul}
\alias{ncase4simul}

\title{
Function estimating effective size (Nc) and sex rate (s) for simulated data
}
\description{
This function simulates the evolution and the sampling of a partially clonal population, and computes raw point estimates of Nc and s for each simulated data set. 
}
\usage{
ncase4simul(param, nb1ClonalCycles, nb2ClonalCycles, nb1Sample, nb2Sample, 
            start.counts, nbPrelCycles, nbSimul, start.param, maxit = 500, trace = TRUE)
}

\arguments{
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
  \item{start.counts}{
Vectors of non-negative integers providing counts of strains used to initialize the simulations.
}
  \item{nbPrelCycles}{
Number of preliminary series of cycles (one series consists of \code{nb1ClonalCycles} clonal cycles, 1 mixed cycle and \code{nb2ClonalCycles} clonal cycles) performed before the collection of the 1st sample.
}
  \item{nbSimul}{
Non-negative integer giving the number of simulated data sets that are generated and for which Nc and s are estimated.
}
  \item{start.param}{
2-column matrix whose rows give initial values for Nc and s that are used in the optimization of the contrast.
}
  \item{maxit}{
Maximum number (non-negative integer) of iterations in the optimization of the contrast. Default to 500.
}
  \item{trace}{
Logical indicating if tracing information are displayed (\code{TRUE}) or not displayed (\code{FALSE}). Default to \code{TRUE}.
}
}
\details{
The optimization method that is used is the Nelder and Mead (1965) method implemented in the \code{optim} function (\code{stats} package).
}
\value{
A 4-column matrix with true Nc (1st col.), true s (2nd col.), the raw estimate of Nc (3rd col.) and the raw estimate of s (4th col.). The raw estimates are the estimates of Nc and s without bias correction. The matrix contains as many lines as the number of simulations requested by the argument \code{nbSimul}.
}
\references{
Ali S, Soubeyrand S, Gladieux P, Giraud T, Leconte M, Gautier A, Mboup M, de Vallavieille-Pope C and Enjalbert J (2015). CloNcaSe: Estimation of sex frequency and effective population size by clonemate re-sampling in partially clonal organisms. Research Report.

Nelder J A and Mead R (1965). A simplex algorithm for function minimization. Computer Journal 7: 308-313.
}
\author{
Samuel Soubeyrand
}

\seealso{
\code{\link{simul}}, \code{\link{ncase}}.
}
\examples{
## count data used to initialize the simulations
counts1=c(20,10,5,1,1,1,1,1)

## simulation and raw estimation
ncase4simul(param=c(1000,0.5),nb1ClonalCycles=17,nb2ClonalCycles=7,
            nb1Sample=200,nb2Sample=200,start.counts=counts1,nbPrelCycles=5,
            nbSimul=100,start.param=expand.grid(c(100,1000),c(0.2,0.8)))
}
\keyword{ ~misc }

