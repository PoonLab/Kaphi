\name{simulate.ode.tree}
\alias{simulate.ode.tree}
\title{
Simulate Tree with Ordinary Differential Equations
}
\description{
Simulates a tree with a given ODE solution and returns an \code{ape} phylo object.
}
\usage{
simulate.ode.tree(sol, sample.times, sample.states, integration.method='rk4', simulate.migrations=FALSE)
}
\arguments{
  \item{sol}{
    Given ODE solution, either as the solution in a returned list from the \code{solve.ode} function or a direct return value from ode() from the \code{deSolve} package.
  }
  \item{sample.times}{
    An n-vector of sample collection times, where n is sample size.
  }
  \item{sample.states}{
    An n times m matrix of sample deme states where m represents the number of demes.
  }
  \item{integration.method}{
    The integrator to use, either a \bold{function} that performs integration, or a \bold{list} of class \code{\link{rkMethod}}, or a \bold{string} ("lsoda", "lsode", "lsodes","lsodar","vode", "daspk", "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" ,"iteration"). Options "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" are the backward differentiation formula, the BDF with diagonal representation of the Jacobian, the (explicit) Adams and the implicit Adams method, and the implicit Adams method with diagonal representation of the Jacobian respectively (see details). The default integrator used is \code{\link{rk4}}. See \code{ode} in the \code{deSolve} package for more information.
  }
  \item{simulate.migrations}{
    Default is FALSE.
  }
}
\details{
See the selected integrator for the additional options.

The default integrator used is \code{rk4}.

The option method = "bdf" provdes a handle to the backward differentiation formula (it is equal to using method = "lsode"). It is best suited to solve stiff (systems of) equations.

The option method = "bdf_d" selects the backward differentiation formula that uses Jacobi-Newton iteration (neglecting the off-diagonal elements of the Jacobian (it is equal to using method = "lsode", mf = 23). It is best suited to solve stiff (systems of) equations.

method = "adams" triggers the Adams method that uses functional iteration (no Jacobian used); (equal to method = "lsode", mf = 10. It is often the best choice for solving non-stiff (systems of) equations. Note: when functional iteration is used, the method is often said to be explicit, although it is in fact implicit.

method = "impAdams" selects the implicit Adams method that uses Newton- Raphson iteration (equal to method = "lsode", mf = 12.

method = "impAdams_d" selects the implicit Adams method that uses Jacobi- Newton iteration, i.e. neglecting all off-diagonal elements (equal to method = "lsode", mf = 13.

For very stiff systems, method = "daspk" may outperform method = "bdf".
}
\value{
Returns a tree (\code{ape} phylo object).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\seealso{
\code{\link{solve.ode}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
