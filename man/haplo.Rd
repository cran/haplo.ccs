\name{haplo}
\alias{haplo}
\title{Specify an Allele Matrix and Inheritance Mode}
\description{
  
  'haplo' specifies an allele matrix and an inheritance mode for the
  'haplo.ccs' model formula.

}
\usage{

haplo(..., mode)

}
\arguments{

  \item{\dots}{a matrix of alleles or list of columns of an allele matrix.
    Each locus on the chromosome has a pair of adjacent columns in the
    allele matrix, so that the number of columns of the allele matrix is
    twice the number of loci on the chromosome. The order of the columns
    corresponds to the order of the loci on the chromosome, and each row
    corresponds to the alleles for each subject.  The alleles should be
    numerically coded, i.e., 1, 2, 3, or 4 for 'A', 'C', 'G', or 'T'.}

  \item{mode}{the inheritance mode, either 'additive', 'dominant', or
    'recessive'. Note that although 'haplo' does not specify a default
    inheritance mode, 'additive' is the default inheritance mode for
    'haplo.ccs'.}
  
}
\value{
  
  A matrix of alleles with the inheritance mode assigned as an
  attribute.

}
\author{

  Benjamin French and Thomas Lumley

  Department of Biostatistics

  University of Washington

}
\seealso{\code{\link{haplo.ccs}}}
\examples{

data(Renin)

## Specify an allele matrix with additive inheritance in a model fit by
## 'haplo.ccs' to estimate haplotype effects.

m1 <- haplo.ccs(case ~ haplo(geno[,1:12], mode = "additive"),
                control = haplo.em.control(min.posterior=1e-4), referent = "223144")

}
\keyword{models}
\keyword{regression}
