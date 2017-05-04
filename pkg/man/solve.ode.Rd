\name{solve.ode}
\alias{solve.ode}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Numerical solution of ODE system}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
solve.ode(expr, t0, t1, x0, parms, time.pts = 2000, integrationMethod = "rk4")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
%%     ~~Describe \code{expr} here~~
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{t1}{
%%     ~~Describe \code{t1} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{parms}{
%%     ~~Describe \code{parms} here~~
}
  \item{time.pts}{
%%     ~~Describe \code{time.pts} here~~
}
  \item{integrationMethod}{
%%     ~~Describe \code{integrationMethod} here~~
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
function (expr, t0, t1, x0, parms, time.pts = 2000, integrationMethod = "rk4") 
{
    demeNames <- expr$demeNames
    nonDemeNames <- expr$nonDemeNames
    m <- length(demeNames)
    mm <- length(nonDemeNames)
    if (length(x0) != m + mm) {
        stop("Initial conditions incorrect dimension", x0, m, 
            mm)
    }
    if (any(!is.element(c(demeNames, nonDemeNames), names(x0)))) {
        stop("Initial conditions vector is incomplete", names(x0), 
            demeNames, nonDemeNames)
    }
    x0 <- x0[c(demeNames, nonDemeNames)]
    if (time.pts < 1) 
        stop("time.pts argument must be >= 1")
    if (t0 > t1) 
        stop("Start time t0 must be less than end time t1")
    times <- get.times(t0, t1, time.pts)
    .birth.matrix <- function(x, t) {
        with(as.list(x), t(matrix(sapply(1:m^2, function(k) eval(expr$pbirths[k])), 
            nrow = m, ncol = m)))
    }
    .migration.matrix <- function(x, t) {
        with(as.list(x), t(matrix(sapply(1:m^2, function(k) eval(expr$pmigrations[k])), 
            nrow = m, ncol = m)))
    }
    tBirths <- function(x, t) {
        colSums(.birth.matrix(x, t))
    }
    tMigrationsIn <- function(x, t) {
        colSums(.migration.matrix(x, t))
    }
    tMigrationsOut <- function(x, t) {
        rowSums(.migration.matrix(x, t))
    }
    tDeaths <- function(x, t) {
        with(as.list(x, t), sapply(1:m, function(k) eval(expr$pdeaths[k])))
    }
    dNonDeme <- function(x, t) {
        with(as.list(x, t), sapply(1:mm, function(k) eval(expr$pndd[k])))
    }
    dx <- function(t, y, parms, ...) {
        dxdeme <- setNames(tBirths(y, t) + tMigrationsIn(y, t) - 
            tMigrationsOut(y, t) - tDeaths(y, t), demeNames)
        dxnondeme <- ifelse(mm > 0, setNames(dNonDeme(y, t), 
            nonDemeNames), NULL)
        list(c(dxdeme, dxnondeme))
    }
    sol <- ode(y = x0, times = times, func = dx, parms = parms, 
        method = integrationMethod)
    index <- nrow(sol):1
    result <- list(times = rev(times), Y = lapply(index, function(i) sol[i, 
        demeNames]), F = lapply(index, function(i) .birth.matrix(sol[i, 
        ], sol[i, 1])), G = lapply(index, function(i) .migration.matrix(sol[i, 
        ], sol[i, 1])), sol = sol)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
