### Application of the SAEM_KM function to the SIR model
#####################################################################

rm(list=objects())

library(deSolve)      # load this library to use the R ode function
library(GillespieSSA) # load this library to use the R ssa function
library(doParallel)   # load this library to parallelize the code

source('KalmanFunctions.R') 
source('ModelFunctions.R')

s0 = 0.9 # initial proportion of susceptible individuals 
i0 = 0.1 # initial proportion of infectious individuals
x0 = c(s0, i0)

lambda = 0.6   # transmission rate
gamma  = 0.4 # recovery rate

p = 0.8 # reporting rate

N  = 10000 # population size
tf = 40    # final time     

U  = 30 # number of epidemics
dt = 2  # time interval between observations

Gamma_1 = log(1 + 0.2 ^ 2)
beta_1  = log(lambda) - 1 / 2 * Gamma_1 

Gamma_2 = 0.8
beta_2  = 1.45

set.seed(50)

# lambda and p are random parameters, named lambda_u and p_u respectively
psi1 = rnorm(U, beta_1, sqrt(Gamma_1)) # h(beta_1, xi_{1,u}) = exp(beta_1 + xi_{1,u}) := lambda_u with targeted mean and variance: 1 and 0.2 ^ 2
psi2 = rnorm(U, beta_2, sqrt(Gamma_2)) # h(beta_2, xi_{2,u}) = 1 / (1 + exp(- beta_2 + xi_{2,u})) := p_u with targeted mean and variance: 0.5 and 0.16 ^ 2

### a- Simulation of the SIR jump process with the GillespieSSA package features

nu = matrix(c(- 1, 0, 1, - 1, 0, 1), nrow = 3, byrow = TRUE) # SIR jumps
a  = c("lambda * S * I / N", "gamma * I")                    # transition rates

Y = matrix(NA, nrow = U, ncol = tf / dt + 1)

for(u in 1 : U){
  out = ssa(
    x0 = c(S = s0 * N, I = i0 * N, R = 0),
    a = a,
    nu = nu,
    parms = c(lambda = exp(psi1[u]), gamma = gamma),
    tf = tf,
    method = ssa.d(),
    verbose = FALSE,
    consoleInterval = 1
  ) 
  
  ### b- Simulation of the observations
  
  # Generation of the observations in regurlarly spaced times 
  
  data = rep(NA, tf / dt + 1)
  S    = rep(NA, tf / dt + 1)
  
  for (k in 1 : length(data)){
    ind  = max(which(out$data[, 1] < (k * dt)))
    S[k] = out$data[ind, 2]   # Markov jump process of the susceptible individuals
    
    if(k == 1){
      data[k] = rbinom(1, s0 * N - S[k], 1 / (1 + exp(- psi2[u])))  # observations
    }else{
      data[k] = rbinom(1, S[k - 1] - S[k], 1 / (1 + exp(- psi2[u])))  # observations 
    }
  }
  
  Y[u, ] = - data / N # the data are normalized by the population size N 
}

### c- SAEM algorithm combined with Kalman filtering techniques (SAEM_KM function)
### lambda and p are (unknown) random effects;  gamma is fixed and known
### intial starting points of the epidemics are fixed and known

# Initialization of the parameters values (fixed effects and variances) for using the SAEM_KM function

theta.init = c(log(0.8), 1.6, 0.4, 1.2) # starting points for the model parameters c(beta_1, beta_2, sqrt(Gamma_1), sqrt(Gamma_2))

# Tuning parameters of the SAEM algorithm (see Appendix C of the manuscript for more details)

fixed_effect = c(FALSE, FALSE) # (lambda, p) is a vector of random parameters
M_max        = 2000            # maximal number of iterations of the SAEM algorithm
M_0          = 1500            # number of burn-in iterations of the SAEM algorithm
nu_0         = 0.6             # real value in [0.5, 1] implied in the expression of the step-size of the SAEM algorithm after a number M_0 of iterations
tau_0        = 0.98            # real value in ]0, 1[ implied in the simulated annealing version of SAEM
K_0          = 0.87            # real value in ]0, 1[ if there exist fixed parameters

res_SAEM = SAEM_KM(theta.init, gamma, dt, N, Y, x0, fixed_effect, M_max, M_0, nu_0, tau_0, K_0,
                   ODE, b, gradb, Sigma, B, Vobs) 

print(res_SAEM$theta[res_SAEM$iter_max, ]) # estimation of the parameters in the following order: beta_1, beta_2, sqrt(Gamma_1), sqrt(Gamma_2)

# Convergence plots of the SAEM algorithm for each model parameter
par(mfrow = c(2, 2))
plot(x = 1 : res_SAEM$iter_max, y = res_SAEM$theta[1 : res_SAEM$iter_max, 1], col = "blue", type = "l", xlab = "Iteration", ylab = expression(beta[1]), ylim = c(- 0.7, 0))
abline(a = mean(psi1), b = 0, col = "red", lty = 2)

plot(x = 1 : res_SAEM$iter_max, y = res_SAEM$theta[1 : res_SAEM$iter_max, 2], col = "blue", type = "l", xlab = "Iteration", ylab = expression(beta[2]), ylim = c(1, 2))
abline(a = mean(psi2), b = 0, col = "red", lty = 2)

plot(x = 1 : res_SAEM$iter_max, y = res_SAEM$theta[1 : res_SAEM$iter_max, 3], col = "blue", type = "l", xlab = "Iteration", ylab = expression(sqrt(Gamma[1])), ylim = c(0, 0.5))
abline(a = sd(psi1), b = 0, col = "red", lty = 2)

plot(x = 1 : res_SAEM$iter_max, y = res_SAEM$theta[1 : res_SAEM$iter_max, 4], col = "blue", type = "l", xlab = "Iteration", ylab = expression(sqrt(Gamma[2])), ylim = c(0.6, 1.5))
abline(a = sd(psi2), b = 0, col = "red", lty = 2)
