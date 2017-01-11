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

test.parse.ode <- function() {
    parms <- list(beta=0.1, mu=0.01, gamma=0.02)
    S<-99
    I<-1
    result <- parse.ode(births, deaths, ndd, migrations)

    checkEquals('I', result$demeNames)
    checkEquals('S', result$nonDemeNames)

    checkEquals(c(0.1*99*1/(99+1)), eval(result$pbirths))
    checkEquals(c(0), eval(result$pmigrations))
    checkEquals(c((0.01+0.02)*1), eval(result$pdeaths))
    checkEquals(c(0.03-0.1*99/100), eval(result$pndd))
}

test.solve.ode <- function() {
    expr <- parse.ode(births, deaths, ndd, migrations)
    parms <- list(beta=0.1, mu=0.01, gamma=0.02)
    x0 <- c(S=999, I=1)
    result <- solve.ode(expr, t0=0, t1=100, x0=x0, parms=parms, time.pts=100, integrationMethod='rk4')

    # check time axis
    checkEquals(100, length(result$times))
    checkEquals(0, result$times[100])
    checkEquals(100, result$times[1])
    checkEquals(100/99, result$times[99])

    checkEquals(100, length(result$Y))
    checkTrue(all(sapply(result$Y, ncol)==1))

    checkTrue(all(sapply(result$F, ncol)==1))
    # TODO: is analytical solution feasible for this problem?
}
