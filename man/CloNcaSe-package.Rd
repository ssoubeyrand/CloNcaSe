\name{CloNcaSe-package}
\alias{CloNcaSe-package}
\alias{CloNcaSe}
\docType{package}
\title{Estimation of effective size and sex rate of a partially clonal population}
\description{
\code{CloNcaSe}: Nc (effective size) and s (sex rate) estimator. Code to estimate effective size (Nc) and sex rate (s) of a partially clonal population sampled at two different times. The samples are supposed to be separated by several clonal cycles and one mixed (i.e. sexual-clonal) cycle.



}
\details{
\tabular{ll}{
Package: \tab CloNcaSe\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-07-07\cr
License: \tab GPL-3\cr
}

The main function of this package is the function  \code{\link{ncase}}.

}
\author{
Samuel Soubeyrand

Maintainer: Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
}
\references{
Ali S, Soubeyrand S, Gladieux P, Giraud T, Leconte M, Gautier A, Mboup M, de Vallavieille-Pope C and Enjalbert J (2015). CloNcaSe: Estimation of sex frequency and effective population size by clonemate re-sampling in partially clonal organisms. Research Report.
}

\keyword{ package }
\seealso{
\code{\link{ncase}}, \code{\link{simul}}, \code{\link{ncase4simul}}.
}
\examples{
## Yellow rust data
data(yellowrust)

## Point estimates
ncase(yellowrust,nb1ClonalCycles=17,nb2ClonalCycles=7,
	start=expand.grid(c(100,1000),c(0.2,0.8)))
}
