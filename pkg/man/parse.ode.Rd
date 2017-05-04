\name{parse.ode}
\alias{parse.ode}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Parse string matrices into unevaluated R expressions}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parse.ode(births, deaths, ndd, migrations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{births}{
%%     ~~Describe \code{births} here~~
}
  \item{deaths}{
%%     ~~Describe \code{deaths} here~~
}
  \item{ndd}{
%%     ~~Describe \code{ndd} here~~
}
  \item{migrations}{
%%     ~~Describe \code{migrations} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (births, deaths, ndd, migrations) 
{
    result <- list(demeNames = rownames(births), nonDemeNames = names(ndd), 
        pbirths = parse(text = t(births)), pmigrations = parse(text = t(migrations)), 
        pdeaths = parse(text = deaths), pndd = ifelse(!is.na(ndd) & 
            length(ndd) > 0, parse(text = ndd), NA))
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
