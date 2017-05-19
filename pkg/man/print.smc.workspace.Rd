\name{print.smc.workspace}
\alias{print.smc.workspace}

\title{
    Print Initialized SMC Workspace
}
\description{
    \code{print.smc.workspace} prints a summary of the data
    contained in an "smc.workspace" object. This includes 
    information about the observed tree and the SMC 
    configurations.
}
\usage{
    print.smc.workspace(workspace, ...)
}

\arguments{
    \item{workspace}{
        An object of class 'smc.workspace'.
}
    \item{\dots}{
      Opptional arguments.
    }
}

\details{}

\value{
    A summary of workspace printed to the screen that contains
    information on the observed tree and the configuration used.
}

\references{}

\author{
    Art Poon and Mathias Renaud
}

\note{}

\seealso{
    \code{\link{init.workspace}} for the building of "workspace".
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

## Initialize workspace
ws <- init.workspace(t1, config, regex=NA)

print.smc.workspace(ws)
}

\keyword{SMC}
\keyword{smc.workspace}
\keyword{workspace}