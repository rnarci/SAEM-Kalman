### Functions implementing the SAEM algorithm combined with Kalman filtering techniques
########################################################################################

### The following functions are generalisable to any compartmental model, the latter being specified by the functions described in ModelFunctions.R

# Kalman: function implementing the Kalman filter (cf Section 3.1.1 in https://doi.org/10.1016/j.csda.2021.107319 )
Kalman = function(eta,   # numerical vector of parameters ruling the dynamics of the epidemic model
                  p,     # real value of the reporting rate in [0,1]
                  dt,    # real value of the sampling interval
                  N,     # population size
                  Y,     # vector of observations (values in Y must be normalized with respect to the population size)
                  x0,    # initial state vector (values in x0 must be normalized with respect to the population size)
                  ODE,   # name of the generic function used to compute the ordinary differential equations of the deterministic epidemic model
                  b,     # name of the generic function used to compute the drift of the diffusion process
                  gradb, # name of the generic function used to compute the gradient of b
                  Sigma, # name of the generic function used to compute the diffusion matrix of the diffusion process
                  B,     # name of the generic function used to compute the projection operator linking the observations to the states of the epidemic model
                  Vobs   # name of the generic function used to compute the variance of the observations
)
{

  n         = length(Y)  # number of observations
  nb_states = length(x0) # number of states in the model
  
  # Key quantities computed by the Kalman filter

  epsilon    = list() # innovation
  Gamma      = list() # innovation covariance
  H          = list() # Kalman Gain
  X_pred     = list() # predicted mean state estimation
  Sigma_pred = list() # predicted error covariance
  
  # Deterministic model
  
  x = abs(ode(func = ODE, y = x0,
              times = seq(0, (n + 1) * dt, dt),
              parms = eta))
  
  x = x[-1, -1] 
  
  h = 0.5 
  if(dt > h){
    x_full = abs(ode(func = ODE, y = x0,
                     times = c(0, seq(dt, (n + 1) * dt, h)),
                     parms = eta)) 
    x_full = x_full[-1, -1] 
    
    Phi_th = list() # list containing each resolvent matrix computed at times dt,dt + h,...,(n + 1) * dt
    for(i in 1 : length(seq(dt, (n + 1) * dt, h))){ 
      gradb_km1_full = gradb(x_full[i, ], eta)
      Phi_th[[i]]    = diag(1, nb_states) + h * gradb_km1_full
    }
  }
  
  # Initialization 
  
  X_pred[[1]]     = matrix(x0, nrow = nb_states, ncol = 1) 
  Sigma_pred[[1]] = diag(1 / N, nb_states) 
  
  # Recursions on the Kalman Filter
  
  for (k in 1 : n)
  {
    
    x_km1     = x[k, ]            # ODE solution at step k - 1
    gradb_km1 = gradb(x_km1, eta) # gradient of the drift function b at step k - 1 
    x_k       = x[k + 1, ]        # ODE solution at step k
    b_km1     = b(x_km1, eta)     # drift function at step k - 1
    Sigma_km1 = Sigma(x_km1, eta) # diffusion matrix at step k - 1
    
    if(dt > h){ 
      A_km1 = res_matrix_h(Phi_th, dt, k, h) 
    }else{
      A_km1 = diag(1, nb_states) + dt * gradb_km1 
    }
    F_k   =  x_k - A_km1 %*% x_km1 
    T_km1 = 1 / N * dt * Sigma_km1 
    
    # Filter equations
    
    epsilon             = Y[k] - B(p) %*% X_pred[[k]]                                    # innovation
    Gamma               = B(p) %*% Sigma_pred[[k]]  %*% t(B(p)) + 1 / N * Vobs(x_km1, p) # innovation covariance
    H                   = A_km1 %*% Sigma_pred[[k]] %*% t(B(p)) %*% solve(Gamma)         # Kalman Gain
    X_pred[[k + 1]]     = F_k + A_km1 %*% X_pred[[k]] + H %*% epsilon                    # predicted mean state estimation
    Sigma_pred[[k + 1]] = (A_km1 - H %*% B(p)) %*% Sigma_pred[[k]] %*% t(A_km1) + T_km1  # predicted error covariance
    
  }
  
  output = list(X_pred = X_pred, Sigma_pred = Sigma_pred, x = x)
  
}

# Likelihood: function to compute the log-likelihood
Likelihood = function(eta,   # numerical vector of parameters ruling the dynamics of the epidemic model
                      p,     # real value of the reporting rate in [0,1]
                      dt,    # real value of the sampling interval
                      N,     # population size
                      Y,     # vector of observations (values in Y must be normalized with respect to the population size)
                      x0,    # initial state vector (values in x0 must be normalized with respect to the population size)
                      ODE,   # name of the function used to compute the ordinary differential equations of the deterministic epidemic model
                      b,     # name of the function used to compute the drift of the diffusion process
                      gradb, # name of the function used to compute the gradient of b
                      Sigma, # name of the function used to compute the diffusion matrix of the diffusion process
                      B,     # name of the function used to compute the projection operator linking the observations to the states of the epidemic  model
                      Vobs   # name of the function used to compute the variance of the observations
)
{
  n = length(Y)
  
  res = Kalman(eta, p, dt, N, Y, x0,
               ODE, b, gradb, Sigma, B, Vobs)
  
  X_pred     = res$X_pred
  Sigma_pred = res$Sigma_pred
  x          = res$x
  
  # Recursive computation of the log-likelihood 
  
  LL = 0 
  
  for(k in 1 : n) 
  {
    s2 = abs(B(p) %*% Sigma_pred[[k]] %*% t(B(p))) + abs(1 / N * Vobs(x[k, ], p)) # variance of Y_k given observations y_{k-1},...,y_0
    m  = B(p) %*% X_pred[[k]]                                           # mean of Y_k given observations y_{k-1},...,y_0 
    LL = LL + log(2 * pi) + 0.5 * log(det(s2)) + 0.5 * t(Y[k] - m) %*% solve(s2) %*% (Y[k] - m)
  }
  return(- LL)
}

# res_matrix_h: function to compute the resolvent matrix for high sampling interval dt  
res_matrix_h = function(Phi_th, # list containing each resolvent matrix computed at times dt,dt + h,...,(n + 1) * dt
                        dt,     # real value of the sampling interval
                        k,      # integer value specifying the time the resolvent matrix is computed (from (k - 1) * dt to k * dt)
                        h = 0.1 # step of integration used for approximation of the resolvent matrix (default to 0.5)
)
{
  A_km1 = Phi_th[[floor(dt / h) * (k - 1) + 1]] 
  
  for(i in 2 : (floor(dt / h))){ 
    A_km1 = A_km1 %*% Phi_th[[floor(dt / h) * (k - 1) + i]]
  }
  return(A_km1)
}

# SAEM_KM: functions implementing the SAEM algorithm (cf Appendix C for the notations of the tuning parameters)
SAEM_KM = function(theta.init,   # numerical vector of starting points of the fixed effects beta and the diagonal of the covariance matrix Gamma
                   gamma,        # real value of 1 / infectious period, with the infectious period is known
                   dt,           # real value of the sampling interval
                   N,            # population size
                   Y,            # vector of observations (values in Y must be normalized with respect to the population size)
                   x0,           # initial state vector (values in x0 must be normalized with respect to the population size)
                   fixed_effect, # logical vector of size the number of parameters, with FALSE = random parameter and TRUE = fixed parameter
                   M_max,        # maximal number of iterations of the SAEM algorithm
                   M_0,          # number of burn-in iterations of the SAEM algorithm
                   nu_0,         # real value in [0.5, 1] implied in the expression of the step-size of the SAEM algorithm after a number M_0 of iterations
                   tau_0,        # real value in ]0, 1[ implied in the simulated annealing version of SAEM
                   K_0,          # real value in ]0, 1[ if there exist fixed parameters
                   ODE,          # name of the generic function used to compute the ordinary differential equations of the deterministic epidemic model
                   b,            # name of the generic function used to compute the drift of the diffusion process
                   gradb,        # name of the generic function used to compute the gradient of b
                   Sigma,        # name of the generic function used to compute the diffusion matrix of the diffusion process
                   B,            # name of the generic function used to compute the projection operator linking the observations to the states of the epidemic  model
                   Vobs          # name of the generic function used to compute the variance of the observations
)
{
  
  no_cores = detectCores() - 1 # number of cores for using a parallelized code
  
  U         = dim(Y)[1]            # number of epidemics
  nb_params = length(theta.init)/2 # number of random parameters
  alpha     = 1                    # initial step size
  
  s1 = t(t(rep(0, nb_params)))
  s2 = t(t(rep(0, nb_params)))
  
  theta     = matrix(NA, nrow = M_max, ncol = length(theta.init))
  theta[1, ] = theta.init
  
  psi = matrix(NA, nrow = nb_params, ncol = U) # initial vector of random parameters
  
  for(a in 1 : nb_params){
    psi[a, ] = rnorm(U, mean = theta[1, a], sd = theta[1, nb_params + a])
  }
  
  stop_success = 0
  
  cl = makeCluster(no_cores, type="FORK")
  
  for(k in 2 : M_max){
    c_psi = matrix(NA, nrow = nb_params, ncol = U) # initial vector of candidate random parameters
    
    for(a in 1 : nb_params){
      c_psi[a, ] = rnorm(U, mean = theta[k - 1, a], sd = theta[k - 1, nb_params + a]) 
    }
    clusterExport(cl, "c_psi", envir = environment())
    clusterExport(cl, "k", envir = environment())
    
    res = parSapply(cl, 1 : U, function(u){ # Metropolis-Hastings algorithm (a single iteration is used)
      for(a in 1:nb_params){
        if(a == 1){
          proba_1 = Likelihood(eta = c(exp(c_psi[1, u]), gamma), p = 1 / (1 + exp(- psi[2, u])), x0 = x0, dt = dt, N = N, Y = Y[u, ],
                               ODE = ODE_SIR, b = b_SIR, gradb = gradb_SIR, Sigma = Sigma_SIR, 
                               B = B_SIR, Vobs = Vobs_SIR)
          proba_2 = Likelihood(eta = c(exp(psi[1, u]), gamma), p = 1 / (1 + exp(- psi[2, u])), x0 = x0, dt = dt, N = N, Y = Y[u, ],
                               ODE = ODE_SIR, b = b_SIR, gradb = gradb_SIR, Sigma = Sigma_SIR, 
                               B = B_SIR, Vobs = Vobs_SIR) 
        }else{
          proba_1 = Likelihood(eta = c(exp(psi[1, u]), gamma), p = 1 / (1 + exp(- c_psi[2, u])), x0 = x0, dt = dt, N = N, Y = Y[u, ],
                               ODE = ODE_SIR, b = b_SIR, gradb = gradb_SIR, Sigma = Sigma_SIR, 
                               B = B_SIR, Vobs = Vobs_SIR)
          proba_2 = Likelihood(eta = c(exp(psi[1, u]), gamma), p = 1 / (1 + exp(- psi[2, u])), x0 = x0, dt = dt, N = N, Y = Y[u, ],
                               ODE = ODE_SIR, b = b_SIR, gradb = gradb_SIR, Sigma = Sigma_SIR, 
                               B = B_SIR, Vobs = Vobs_SIR)
        }
        
        log_rate_MH = min(0, proba_1 - proba_2) # log(rate of acceptation of the Metropolis-Hastings algorithm)
        
        if(exp(log_rate_MH) > runif(1, 0, 1)){
          psi[a, u] = c_psi[a, u]
        }
      }
      
      c(psi[, u])
    })
    
    psi = res[1 : nb_params, ]
    clusterExport(cl, "psi", envir = environment())
    
    # Computation of the sufficient statistics
    stat1 = t(t((rep(NA, nb_params))))
    stat2 = t(t((rep(NA, nb_params))))
    for(a in 1 : nb_params){
      stat1[a, 1] = sum(psi[a, ])
      stat2[a, 1] = sum(psi[a, ] ^ 2)
    }
    
    s1 = s1 + alpha * (stat1 - s1)
    s2 = s2 + alpha * (stat2 - s2)
    
    # Update of the model parameters using the simulated annealing version of SAEM
    theta[k, 1 : nb_params] = s1[1 : nb_params, 1] / U
    for (a in 1 : nb_params){
      theta[k, nb_params + a] = (k <= M_0) * max(tau_0 * theta[k - 1, nb_params + a], sqrt(s2[a, 1] / U - (s1[a, 1] / U) ^ 2)) + (k > M_0) * sqrt(s2[a, 1] / U - (s1[a, 1] / U) ^ 2)
    }
    
    if(k <= M_0){
      alpha = 1
    }else{
      alpha = 1 / (k - M_0) ^ nu_0
      if(length(which(fixed_effect == TRUE)) != 0){ # update for the fixed parameters
        theta[k, nb_params + which(fixed_effect == TRUE)] = K_0 * theta[k - 1, nb_params + which(fixed_effect == TRUE)]
      }
      
      stop_params = c()
      for(a in 1 : (2 * nb_params)){
        stop_params[a] = abs(theta[k, a] - theta[k - 1, a]) / abs(theta[k, a]) # stopping criterion of the SAEM algorithm for each model parameter
      }
      
      if(length(which(fixed_effect == TRUE)) != 0){
        stop_params = stop_params[ - (nb_params + which(fixed_effect == TRUE))]
      }
      
      if(max(stop_params)<0.0005){
        stop_success = stop_success + 1
      }else{
        stop_success = 0
      }
    }
    
    clusterExport(cl, "theta", envir = environment())
    print(k)
    
    if(stop_success == 20){
      iter_max = k
      break
    }else{
      iter_max = M_max
    }
  }
  
  stopCluster(cl)
  
  res_SAEM = list(theta = theta, iter_max = iter_max)
  
  return(res_SAEM)
}

