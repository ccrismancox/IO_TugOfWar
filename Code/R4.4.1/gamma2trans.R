gamma2trans <- function(gamma, sigma, mars.states, d=.1, discretize=TRUE, F=pnorm, sweep=TRUE, bound=0.05){
  #' Function to input 1st stage parameters and return a transition matrix
  #' Assumes the transition model is of the form 
  #' $$ s^{t+1} = \gamma_0 + a_\text{Hamas} \gamma_1 + a_\text{Fatah}\gamma_2 + 
  #'               \text{L.states} \gamma_3 + 
  #'               a_\text{Hamas}\cdot\text{L.states}\gamma_4+
  #'               a_\text{Fatah}\cdot\text{L.states}\gamma_5+
  #'               \epsilon$$
  #' Inputs: 
  #'
  #'     - Gamma: 1st stage regression coefficients.  Should be of the order
  #'
  #'         1. Intercept
  #'         2. Hamas attacks
  #'         3. Fatah attacks
  #'         4. lag.state parameter 
  #'         5. Hattacks:lag.states
  #'         6. Fattacks:lag.states
  #'         
  #' If you're not including interactions, set these parameters to zero.
  #'     
  #'     - sigma: Standard error of regresssion
  #'     - mars.states: States from the measurement model
  #'     - d: the space between states
  #'
  # functions i need
  
  
  interiorTransition <- function(to, mu, sigma, d){
    hi <- F((to + d/2 - mu ) / sigma)
    lo <- F((to - d/2 - mu ) / sigma)
    return(hi-lo)
  }
  
  if (discretize){
    bounds <- quantile( mars.states, c(bound,1-bound))
    states <- seq(from=bounds[1], to=bounds[2], by=d)
  } else{
    states <- mars.states
  }
  
  G <- expand.grid(states, c(0,1), c(0,1))
  names(G) = c("state", "aH", "aF")
  G <- G[order(G$state, G$aH),]
  rownames(G) <- c()
  nG <- dim(G)[1]
  # MU <- cbind(1,  G$aH, G$aF,G$state, G$state*G$aH, G$state*G$aF) %*% c(mod0$coef[1:3],1, mod0$coef[4:5])
  # sigma <- summary(mod0)$sigma
  MU <- cbind(1,  G$aH, G$aF,G$state, G$state*G$aH, G$state*G$aF) %*% gamma
  Trans <- t(mapply(function(mu){interiorTransition(states, mu, sigma, d)}, MU))
  Trans[,1] <- F((states[1] + d/2 - MU ) / sigma) 
  Trans[,length(states)] <- 1- F((states[length(states)]  - d/2 - MU) / sigma)
  
  # did it work?
  min(rowSums(Trans))
  max(rowSums(Trans))
  max(abs(rowSums(Trans) - 1))
  
  if(sweep){
    # should we do this?
    Trans[Trans < 1e-8] <- 0
    Trans <- sweep(Trans, 1, rowSums(Trans), "/")
  }
  return(Trans)
}
