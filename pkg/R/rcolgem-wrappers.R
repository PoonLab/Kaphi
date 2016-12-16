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

rcolgem.SI <- function(theta, nsim, n.tips, labels=NA, seed=NA, fgyResolution=500, integrationMethod='adams') {
	"
	Use rcolgem to simulate coalescent trees under susceptible-infected (SI)
	model.
	@x
	@param fgyResolution : time resolution of ODE solution
	@param integrationMethod : method for numerical solution of ODE
	@param t.end : boundary condition for ODE solution, time scale of simulation
	@param N : total population size (constant)
	@param beta : transmission rate
	@param gamma : additional mortality from infection
	@param mu : baseline mortality rate
	"
	t0=0
	
	## initial population frequencies
	S=N-1  # number of susceptibles
	I=1  # epidemic originates from a single infection
	x0 <- c(I=I, S=S)  # state vector
	if (any(x0<0)) {
		stop('Population sizes cannot be less than 0.')
	}
	
	# for parsing R expressions representing ODE system
	params <- list(beta=beta, gamma=gamma, mu=mu)  
	if (any(parms<0)) {
		stop ("No negative values permitted for model rate parameters.")
	}
	
	## define ODE system
	
	# demes are subpopulations from which we can sample virus
	demes <- c('I')
	
	# birth is the rate of lineage splitting - in this case, infection of a susceptible
	births <- rbind(c('parms$beta*S*I / (S+I)'))
	rownames(births)=colnames(births) <- demes
	
	# migration is state transition of a lineage without splitting
	migrations <- rbind(c('0'))
	rownames(migrations)=colnames(migrations) <- demes
	
	nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$mu*S + (parms$mu+parms$gamma)*I', 
	'-S*(parms$beta*I) / (S+I)')
	names(nonDemeDynamics) <- 'S'
	
	
}