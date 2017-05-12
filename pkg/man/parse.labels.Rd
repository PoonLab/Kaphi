\name{parse.labels}
\alias{parse.labels}

\title{
    Parse labels for \code{init.workspace}
}

\description{
    \code{parse.labels} compares a list of tip labels from a
    phylogenetic tree to a regular expression to parse labels
    to be used by \code{init.workspace}.
}

\usage{
    parse.labels(labels, regex)
}

\arguments{
    \item{labels}{
        A character vector containing the tip labels of obs.tree.
    }
    \item{regex}{
        A character string containing a \link{regular expression}.
    }
}

\details{}

\value{
    A character vector containing parsed labels.
}

\references{}

\author{
    Art Poon and Mathias Renaud
}

\note{
    Perl-compatible regular expressions should be used.
}

\seealso{
    \code{\link{init.workspace}} for use of output.
}

\examples{
## Set labels:
labels <- c('1_A', '2_B', '3_A', '4_')

## Set regular expression:    
regex <- '_([A-Z]*)$'

result <- parse.labels(labels, regex)
}

\keyword{labels}
\keyword{workspace}