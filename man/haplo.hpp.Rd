\name{haplo.hpp}
\alias{haplo.hpp}
\title{Obtain Summaries of Phase Ambiguity}
\description{
  
  'haplo.hpp' obtains summary statistics of phase ambiguity. The
  proportion of subjects whose highest posterior diplotype probability
  is greater than or equal to a specified probability is reported.

}
\usage{

haplo.hpp(model, prob=0.95)

}
\arguments{

  \item{model}{a fitted model of class '"haplo.ccs"'.}

  \item{prob}{the probability to which to compare the highest posterior
    diplotype probability for each subject. Note that the default
    probability is 0.95. Note also that either a single probability or a
    list of probabilities may be specified.}

}
\value{
  
  The proportion of subjects whose highest posterior diplotype
  probability is greater than or equal to the specified probability or
  probabilities.

}
\author{

  Benjamin French and Thomas Lumley, University of Washington

}
\references{



}
\seealso{

  \code{\link{haplo}},
  \code{\link{haplo.ccs}}

}
\examples{

data(renin)

haplo.hpp(model=haplo.ccs(case ~ haplo(geno)), prob=c(0.90, 0.95))

}
\keyword{models}
\keyword{regression}
