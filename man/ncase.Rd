\name{ncase}
\alias{ncase}

\title{
Function for estimating effective size (Nc) and sex rate (s)
}
\description{
This function computes point estimates of Nc and s for a partially clonal population, confidence intervals (using parametric bootstrap) and the coverage rates of the intervals. The parametric bootstrap is also used to propose a bias correction of the point estimates and their confidence intervals. 
}
\usage{
ncase(counts, nb1ClonalCycles, nb2ClonalCycles, start.param, maxit = 500, 
	nbBootstrap = 0, conf.level = 0.95, nbCoverage = 0, trace = TRUE)
}

\arguments{
  \item{counts}{
2-column matrix whose 1st (resp. 2nd) column gives the counts of strains in sample 1 (resp. 2).
}
  \item{nb1ClonalCycles}{
Non-negative integer indicating the number of clonal cycles separating the 1st sample and the mixed (i.e. clonal-sexual) cycle.
}
  \item{nb2ClonalCycles}{
Non-negative integer indicating the number of clonal cycles separating the mixed (i.e. clonal-sexual) cycle and the 2nd sample.
}
  \item{start.param}{
2-column matrix whose rows give initial values for the effective size (Nc) and the sex rate (s) that are used in the optimization of the contrast.
}
  \item{maxit}{
Maximum number (non-negative integer) of iterations in the optimization of the contrast. Default to 500.
}
  \item{nbBootstrap}{
Number (non-negative integer) of bootstrap replications to assess confidence intervals (based on a normal approximation) and to carry out a bias correction. Default to 0 (no bootstrap procedure).
}
  \item{conf.level}{
Confidence level (between 0 and 1) of the marginal confidence intervals for Nc and s. Default to 0.95. Used if \code{nbBootstrap}>0.
}
  \item{nbCoverage}{
Number of (non-negative integer) iterations to assess the coverage rates of the marginal confidence intervals. Default to 0 (no coverage assessment). Used if \code{nbBootstrap}>0.
}
  \item{trace}{
Logical indicating if tracing information are displayed (\code{TRUE}) or not displayed (\code{FALSE}). Default to \code{TRUE}. Tracing information consist of a message every 100 iterations in the bootstrap procedure and in the assessment of the coverage rates of the confidence intervals.
}
}
\details{
The optimization method that is used is the Nelder and Mead (1965) method implemented in the \code{optim} function (\code{stats} package).
}
\value{
A list whose items are:
\item{bootstrap.estim}{(if \code{nbBootstrap>0}) Estimates of Nc (1st col.) and s (2nd col.) obtained for bootstrap replications.}
\item{cover}{(if \code{nbCoverage>0}) Coverage rates of the marginal intervals for Nc and s when no bias correction is applied (1st col.) and when bias correction is applied (2nd col.).}
\item{optim.output}{A list consisted of (i) the output of the \code{optim} function (\code{optim} is applied to the logarithm of Nc and the logit of s) and (ii) the optimal values of Nc and s (transformed back to the original units of Nc and s).}
\item{estimates}{A numeric vector providing the estimates of Nc and s (if \code{nbBootstrap=0}) or a matrix providing the estimates of Nc and s, their bias-corrected versions, and the associated confidence intervals (if \code{nbBootstrap>0}).}
}
\references{
	Ali S, Soubeyrand S, Gladieux P, Giraud T, Leconte M, Gautier A, Mboup M, de Vallavieille-Pope C and Enjalbert J (2015). CloNcaSe: Estimation of sex frequency and effective population size by clonemate re-sampling in partially clonal organisms. Research Report.
	}
\author{
Samuel Soubeyrand
}

\seealso{
\code{\link{simul}}, \code{\link{ncase4simul}}.
}
\examples{
## Yellow rust data
data(yellowrust)

## Point estimates
ncase(yellowrust,nb1ClonalCycles=17,nb2ClonalCycles=7,
	start=expand.grid(c(100,1000),c(0.2,0.8)))

## Point estimates, Parametric bootstrap, Bias correction 
## and confidence intervals
\donttest{
ncase(yellowrust,nb1ClonalCycles=17,nb2ClonalCycles=7,
	start=expand.grid(c(100,1000),c(0.2,0.8)),
	nbBootstrap=5000,conf.level=0.95)
}

## Point estimates, Parametric bootstrap, Bias correction
## and confidence intervals and their coverage
\donttest{
ncase(yellowrust,nb1ClonalCycles=17,nb2ClonalCycles=7,
	start=expand.grid(c(100,1000),c(0.2,0.8)),
	nbBootstrap=5000,conf.level=0.95,nbCoverage=100)
}
}
\keyword{ ~misc }

