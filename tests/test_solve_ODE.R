# This file is part of Kaphi.

# Kaphi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Kaphi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.

require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/SI-matrices.R')
source('tests/fixtures/kingman.R')

test.parse.ode <- function() {
    parms <- list(beta=0.1, mu=0.01)
    S<-99
    I<-1
    result <- parse.ode(births, deaths, ndd, migrations)

    checkEquals('I', result$demeNames)
    checkEquals('S', result$nonDemeNames)

    checkEquals(c(0.1*99*1/100), eval(result$pbirths))
    checkEquals(c(0), eval(result$pmigrations))
    checkEquals(c((0.01)*1), eval(result$pdeaths))
    checkEquals(c(0.01-0.1*99/100), eval(result$pndd))
}

test.get.times <- function() {
    result <- .get.times(0, 5, 9)
    expected <- c(0, 0.625, 1.25, 1.875, 2.5, 3.125, 3.75, 4.375, 5)
    checkEquals(expected, result)
}

test.solve.ode <- function() {
    expr <- parse.ode(births, deaths, ndd, migrations)
    parms <- list(beta=0.1, mu=0.01)
    x0 <- c(S=999, I=1)
    result <- solve.ode(expr, t0=0, t1=500, x0=x0, parms=parms, time.pts=100, integrationMethod='rk4')

    # check time axis
    checkEquals(100, length(result$times))
    checkEquals(0, result$times[100])
    checkEquals(500, result$times[1])
    checkEquals(500/99, result$times[99])

    checkEquals(100, length(result$Y))

    # check equilibrium solution for SIS model
    checkEquals((1-0.01/0.1)*1000., result$Y[[1]], checkNames=FALSE)

    # pure birth SI model without dampening by S is exponential growth
    result <- unlist(k.sol$Y)
    expected <- exp(k.parms$beta*unlist(k.sol$times))
    # rk4 tends to underestimate I(t) -- improves at higher resolutions
    checkEquals(expected, result, tolerance=1e-6, checkNames=FALSE)
}


