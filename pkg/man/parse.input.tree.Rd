\name{parse.input.tree}
\alias{parse.input.tree}

\title{
    Parse Input Tree for \code{init.workspace}
}

\description{
    \code{parse.input.tree} validates the config argument,
    checks the class of obs.tree, and preprocesses obs.tree
    for use in init.workspace.
}

\usage{
    parse.input.tree(obs.tree, config)
}

\arguments{
    \item{obs.tree}{
        An object of class "phylo".
    }
    \item{config}{
        An object of class "smc.config".
    }
}

\details{
    The smc.config object is validated by checking the presence
    of required config arguments. The user is warned if any
    required arguments are missing from config.
    
    If obs.tree is not of class "phylo" the function will attempt
    to coerce the object into a "phylo" object. If this fails, the
    user is warned.
    
    obs.tree is preprocessed by ladderizing, rescaling, and 
    computing the self-kernel for the tree.
}

\value{
    An object of class "phylo". The obs.tree after preprocessing.
}

\references{}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{init.workspace}} for use of output.
    \code{\link{tree.kernel}} for \code{preprocess.tree} function.
}

\examples{
## Set observed tree:
t1 <- read.tree(text="(A:0.1,B:0.2):0;")

## Set SMC configurations:
config <- list(
     params=c('N'),
     nparticle=10,
     nsample=3,
     decay.factor=1.0,
     rbf.variance=2.0,
     sst.control=1,
     norm.mode='NONE'
)
class(config) <- 'smc.config'

result <- parse.input.tree(t1, config)
}

\keyword{workspace}
\keyword{obs.tree}