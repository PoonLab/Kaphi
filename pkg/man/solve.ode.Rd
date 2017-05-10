\name{solve.ode}
\alias{solve.ode}
\alias{.birth.matrix}
\alias{.migration.matrix}
\alias{tBirths}
\alias{tMigrationsIn}
\alias{tMigrationsOut}
\alias{tDeaths}
\alias{dNonDeme}
\alias{dx}

\title{
    Numerical Solution of ODE System
}

\description{
    \code{solve.ode()} solves a system of ordinary differential
    equations using the parsed expressions from the model, 
    as given by \code{parse.ode()}.
    
    Derived from rcolgem:\code{make.fgy()} (\link{https://github.com/rforge/colgem}).
}

\usage{
    solve.ode(expr, t0, t1, x0, parms, time.pts = 2000, integrationMethod = "rk4")
}

\arguments{
    \item{expr}{
        parsed expressions from \code{parse.ode()}.
    }
    \item{t0}{
        initial time for ODE solution.
    }
    \item{t1}{
        end time for ODE solution.
    }
    \item{x0}{
        initial state vector
    }
    \item{parms}{
        model parameter settings
    }
    \item{time.pts}{
        number of time points in range [t0, t1]
    }
    \item{integrationMethod}{
        the integrator to use; passed to \code{ode()}
    }
}

\details{}

\value{
    A list containing (K = time.pts):
    \item{times}{
        (K) time points of numerical solutions in
        reverse order.
    }
    \item{Y}{
        K (m) vectors of deme population sizes.
    }
    \item{F}{
        K (m x m) matrices of birth rates.
    }
    \item{G}{
        K (m x m) matrices of migration rates.
    }
    \item{sol}{
        return value of \code{ode()}.
    }
}

\references{
    Volz, E. M. (2012) Complex population dynamics and the coalescent under neutrality. \emph{Genetics}, \bold{190}(1), 187--201.
}

\author{
    Art Poon, Mathias Renaud
}

\note{
    A deme refers to an isolated subpopulation. For example, in 
    a standard SIS model, two demes are susceptible ("S") 
    and infected ("I").
}

\seealso{
    \code{\link{parse.ode}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
## rcolgem (https://github.com/rforge/colgem) definition 
## for a standard SIS model with recovery (single-deme model)

## Model parameters
parms <- list(beta=0.1, mu=0.01)

## "Infected" deme name
demes <- ("I")

births <- matrix(c('parms$beta*S*I/(S+I)'), nrow=1, ncol=1)
rownames(births)=colnames(births) <- demes

# migrations undefined for single-deme models
migrations <- matrix(c('0'), nrow=1, ncol=1)
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('parms$mu*I')
names(deaths) <- demes

# non-deme dynamics
ndd <- c(S='(parms$mu)*I - parms$beta*S*I/(S+I)')

expr <- parse.ode(births, deaths, ndd, migrations)
x0 <- c(S=999, I=1)
sol <- solve.ode(expr, t0=0, t1=100, x0=x0, parms=parms, time.pts=100, integrationMethod='rk4')
}

\keyword{math}