### SIR model with lambda the transmission rate and gamma the recovery rate

# Deterministic system 

ODE_SIR = function(times, # finite time interval for the resolution of the ordinary differential equations
                   y,     # vector of the initial (state) values for the ODE system
                   parms  # numerical vector of parameters ruling the dynamics of the epidemic model
)
{ 
  
  lambda = parms[1]
  gamma  = parms[2]
  
  s = y[1]
  i = y[2]
  
  ds = - lambda * s * i
  di = lambda * s * i - gamma * i
  
  list (c(ds, di))
}

# Drift function

b_SIR = function(x,  # solution of the ordinary differential equations of the deterministic epidemic model at a given time point
                 eta # numerical vector of parameters ruling the dynamics of the epidemic model
)
{ 
  
  lambda = eta[1]
  gamma  = eta[2]
  
  s = x[1]
  i = x[2]
  
  b = matrix(0, nrow = 2, ncol = 1)
  
  b[1, 1] = - lambda * s * i
  b[2, 1] = lambda * s * i - gamma * i
  
  return(b)
}

# Gradient of the drift function

gradb_SIR = function(x,  # solution of the ordinary differential equations of the deterministic epidemic model at a given time point
                     eta # numerical vector of parameters ruling the dynamics of the epidemic model
)
{ 
  lambda = eta[1]
  gamma  = eta[2]
  
  s = x[1]
  i = x[2]
  
  gradb = matrix(0, 2, 2)
  
  gradb[1, 1] = - lambda * i
  gradb[1, 2] = - lambda * s
  gradb[2, 1] = lambda * i
  gradb[2, 2] = lambda * s - gamma
  
  return(gradb)
}

# Diffusion matrix

Sigma_SIR = function(x,  # solution of the ordinary differential equations of the deterministic epidemic model at a given time point
                     eta # numerical vector of parameters ruling the dynamics of the epidemic model
)
{ 
  
  lambda = eta[1]
  gamma  = eta[2]
  
  s = x[1]
  i = x[2]
  
  S = matrix(0, 2, 2)
  
  S[1, 1] = lambda * s * i
  S[1, 2] = - lambda * s * i 
  S[2, 1] = - lambda * s * i
  S[2, 2] = lambda * s * i + gamma * i

  return(S)
}

# Projection operator linking the observations to the states of the epidemic model

B_SIR = function(p # real value of the reporting rate in [0,1]
)
{ 
  return(matrix(c(p, 0), nrow = 1, ncol = 2))
}

# Variance of the observations

Vobs_SIR = function(x, # solution of the ordinary differential equations of the deterministic epidemic model at a given time point
                    p  # real value of the reporting rate in [0,1]
)
{ 
  
  s = x[1]
  i = x[2]
  
  return(p * (1 - p) * s)
}


