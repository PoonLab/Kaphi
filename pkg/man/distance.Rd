\name{distance}
\alias{distance}

\title{
    Compute Distance Between Two Trees
}

\description{
    \code{distance} calculates the kernel distance between two
    trees, t1 and t2, using the kernel parameters set in config.
    
    This function uses the kernel function (\code{tree.kernel})
    as a similarity measure for the shape of phylogenetic trees.
}

\usage{
    distance(t1, t2, config)
}

\arguments{
    \item{t1}{
        An object of class 'phylo'.
}
    \item{t2}{
        An object of class 'phylo'.
}
    \item{config}{
        An object of class 'smc.config'.
}
}

\details{}

\value{
    A numeric vector containing the mesaure of distance between
    the two trees.
}

\references{
    Del Moral, P., Doucet, A., & Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}(5), 1009--1020.
    
    Poon, A. F., Walker, L. W., Murray, H., McCloskey, R. M., Harrigan, P. R., & Liang, R. H. (2013) Mapping the shapes of phylogenetic trees from human and zoonotic RNA viruses. \emph{PLoS one}, \bold{8}(11), e78122.
}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{tree.kernel}} for the kernel distance calculation.
    
    \code{\link{load.config}} for the configuration of kernel parameters.
    
    \code{\link{initialize.smc}} for use of this function to populate 'dists'.
}

\examples{
## Set trees to compare
t1 <- read.tree(text="(A:0.1,B:0.2):0;")
t2 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")

## Set kernel configurations
config <- list(
     decay.factor=1.0,
     rbf.variance=2.0,
     sst.control=1,
     norm.mode='NONE'
)
class(config) <- 'smc.config'

## Compute tree self kernels
pt1 <- parse.input.tree(t1, config)
pt2 <- parse.input.tree(t2, config)

result <- distance(pt1, pt2, config)
}

\keyword{Distance}
\keyword{Kernel}