\name{rescale.tree}
\alias{rescale.tree}
\alias{RESCALE.TREE}
\title{
    Rescales the edge length of the tree relative to the specified mode
}
\description{
    Tree edge length is rescaled as a ratio of the current tree edge length over the specified mode of tree lengths. Mode can be specified as the mean, median, max, or none.
}
\usage{
    rescale.tree(tree, mode)
}
\arguments{
    \item{tree}{
        An R object of type tree.  
    }
    \item{mode}{
        A \code{character} variable mode to be specified for scaling the tree edge lengths.  
        Written as \sQuote{mean}, \sQuote{median}, \sQuote{max}, or \sQuote{none}
    }
}  
\details{
    An explicit description of the input format for the `tree` object is shown under \sQuote{Examples} below.
}  
\value{
    Returns an object of tyoe tree with modified edge length characteristic
}  
\author{
    author(s) go here
}  
\seealso{
    \code{\link{tree.kernel}
    }
}  
\examples{
}      
\keyword{file}
\note{}

