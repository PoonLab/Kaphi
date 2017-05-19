\name{parse.ode}
\alias{parse.ode}

\title{
    Parse String Matrices into Unevaluated R Expressions
}

\description{
    \code{parse.ode} returns a list of deme labels and parsed 
    but unevaluated R expressions from the chosen model to be 
    used by \code{solve.ode}.
}

\usage{
    parse.ode(births, deaths, ndd, migrations)
}

\arguments{
    \item{births}{
        a character matrix containing the expression describing births
    }
    \item{deaths}{
        a character containing the expression describing deaths
    }
    \item{ndd}{
        a character containing the expression describing non-deme dynamics
    }
    \item{migrations}{
        a character matrix containing the expression describing migrations
    }
}

\details{}

\value{
    A list containing:
    \item{demeNames}{
        a vector containing the name of each deme
    }
    \item{nonDemeNames}{
        a vector containing the name of each non-deme
    }
    \item{pbirths}{
    the parsed, unevaluated expression describing births
    }
    \item{pmigrations}{
    the parsed, unevalutated expression describing migrations
    }
    \item{pdeaths}{
    the parsed, unevaluated expression describing deaths
    }
    \item{pndd}{
    the parsed, unevaluated expression describing non-deme dynamics
    }
}

\references{}

\author{
    Art Poon, Mathias Renaud
}

\note{
    A deme refers to an isolated subpopulation. For example, in 
    a standard SIS model, two demes are susceptible ("S") 
    and infected ("I").
}

\seealso{
\code{\link{solve.ode}}
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
}

\keyword{}