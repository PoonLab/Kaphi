library(deSolve)

## SIR model w/out vital dynamics
SIR.nondynamic <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / (S +I)
    dI <- beta * S * I / (S + I) - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# testing SIR.nondynamic
init <- c(S=(1 - 1e-6), I=1e-6, 0.0)
parameters <- c(beta = 1.4247, gamma = 0.14286)
times <- seq(0, 70, by = 1)
out <- as.data.frame(ode(y=init, times = times, func = SIR.nondynamic, parms = parameters))
out$time <- NULL

matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)


######################################################################################################################
## SIR model w/ vital dynamics and constant population
SIR.dynamic <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- delta - mu * S - beta * S * I
    dI <- beta * S * I - (gamma + mu) * I
    dR <- gamma * I - mu * R
    
    return(list(c(dS, dI, dR)))
  })
}

# testing SIR.dynamic
init <- c(S=(1 - 1e-6), I=1e-6, R=0.0)
parameters <- c(beta = 1.4247, gamma = 0.14286, delta = 0.08, mu = 0.06)
times <- seq(0, 70, by = 1)
out <- as.data.frame(ode(y=init, times = times, func = SIR.dynamic, parms = parameters))
out$time <- NULL

matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)


######################################################################################################################
## SIS model
SIS <- function(time,state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I / (S + I) + gamma * I
    dI <- (beta * (S + I) - gamma) * I - beta * I^2
    
    return(list(c(dS, dI)))
  })
}

# testing SIS
init <- c(S=(1 - 1e-6), I=1e-6)
parameters <- c(beta = 1.4247, gamma = 0.14286, delta = 0.15)
times <- seq(0, 70, by = 1)
out <- as.data.frame(ode(y=init, times = times, func = SIS, parms = parameters))
out$time <- NULL

matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovered-Susceptibles", main = "SIS Model", lwd = 1, lty = 1, bty = "l", col = 2:3)
legend(40, 0.7, c("Susceptibles", "Infecteds"), pch = 1, col = 2:3)


######################################################################################################################
## SEIR model, assuming presence of vital dynamics w/ birth rate equal to the death rate
SEIR <- function(time, state, parameters){
  with(as.list(c(state, parameters)), {
    dS <- mu * (S + I) - mu * S - beta * S * I / (S + I)
    dE <- beta * S * I / (S + I) - (mu + alpha) * E
    dI <- alpha * E - (gamma + mu) * I 
    dR <- gamma * I - mu * R
    
    return(list(c(dS, dE, dI, dR)))
  })
}

# testing SEIR
init <- c(S=(1 - 1e-6), E=0, I=1e-6, R=0.0)
parameters <- c(beta = 1.4247, gamma = 0.14286, delta = 0.15, mu = 0.15, alpha = 5)
times <- seq(0, 70, by = 1)
out <- as.data.frame(ode(y=init, times = times, func = SEIR, parms = parameters))
out$time <- NULL

matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SEIR Model", lwd = 1, lty = 1, bty = "l", col = 2:5)
legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds", "Exposed"), pch = 1, col = 2:5)