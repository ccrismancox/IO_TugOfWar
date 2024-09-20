liml.fit <- function( X1, X2, Z, y, C=0, vcov=c("robust", "bekker", "CSE", "iid", "none"), fixk=FALSE){
  #' X1: an N by K* matrix of EXOGENOUS variables (include constant)
  #' X2: an N by K** matrix of ENDOGENOUS variables
  #' Z: an N by L matrix of INSTRUMENTS (including self instruments i.e., all of X1)
  #' y: length N vector of outcomes
  #' C: adjustment to k estimator, 0 is liml 1 is Fuller's adjustment
  #' vcov: Do you want robust standard errors, homoskedastic, or none
  
  vcov <- match.arg(vcov)
  if(var(X1[,1])==0 || tolower(colnames(X1)[1]) %in% c("constant", "intercept", "(intercept)", "(constant)")){
    X <- cbind(X1[,1,drop=FALSE], X2, X1[,-1, drop=FALSE])  
  }else{
    X <- cbind(X2, X1)  
  }
  xy <- cbind(y,X2) #endogneous variables
  N <- length(y)
  K <- ncol(Z)
  G <- ncol(X)
  MZ <- -Z %*%solve(crossprod(Z)) %*% t(Z)
  PZ <- -1*MZ
  diag(MZ) <- diag(MZ) + 1
  xy_MZ_xy <- t(xy) %*% MZ %*% xy
  

  
  MX1 <- -X1 %*%solve(crossprod(X1)) %*% t(X1)
  diag(MX1) <- diag(MX1) + 1
  xy_MX1_xy <- t(xy) %*% MX1 %*% xy
  xy_MZ_xy.inv <- solve(xy_MZ_xy)
  xyMZxy.5 <- t(chol(xy_MZ_xy.inv))
   
  if(!fixk){
    k <- min(eigen(t(xyMZxy.5) %*% xy_MX1_xy %*% xyMZxy.5)$values) - C/(N-K)
  }else{
    k <- fixk
    }
  X_kM_X_inv <- solve(crossprod(X) - k*t(X) %*% MZ %*% X)
  out <- y - k* (MZ %*% y)
  liml.est <- drop(X_kM_X_inv%*%crossprod(X, out))
  liml.est <- drop(liml.est)
  names(liml.est) <- colnames(X)
  

  if(vcov=="none"){
    output <- list(coef.table = cbind(liml.est),
                   VCOV = NULL,
                   k=k)
    return(output)
  }
  if(vcov=="robust"){ #Robust to heteroskedasticity 
      # taken from Hansen
      e.hat <- drop(y-X %*% liml.est)
      
      Xtilde <- cbind(X1, X2-k*MZ%*%X2)
      V2 <- solve(t(Xtilde) %*% X) %*% crossprod(Xtilde*e.hat) %*% solve(t(X) %*% Xtilde)
      se.liml <- sqrt(diag(V2))
      z.tests <- liml.est/se.liml
      p.val <- 2*pnorm(abs(z.tests), lower=F)
      output <- list(coef.table = cbind(liml.est, se.liml,z.tests,p.val),
                     VCOV = V2,
                     k=k)
  }
  if(vcov=="bekker"| vcov=="CSE"){
      # Robust to heteroskedasticity 
      # Taken from Hansen, Hausman, and Newey (2008)
    # For Hausman and Newey's notation
    # a.tilde <- min(eigen(solve(crossprod(X.bar)) %*% ( t(X.bar) %*% P %*% X.bar))$values)
    # a.hat <- a.tilde
    # liml2 <-solve(t(X) %*% PZ %*% X - a.hat *crossprod(X)) %*% (t(X) %*% PZ %*% y -a.hat *t(X)%*%y)
    # a.hat <- (a.tilde-((1-a.tilde)*C)/nrow(X))/  ((1-((1-a.tilde)*C)/nrow(X)))
    # fuller <-  solve(t(X) %*% PZ %*% X - a.hat *crossprod(X)) %*% (t(X) %*% PZ %*% y -a.hat *t(X)%*%y)
    
    u <- y - X  %*% liml.est
    sigma2 <- drop(crossprod(u))/(N-G)
    ahat <- drop(t(u) %*% PZ %*%u / drop(crossprod(u)))
    GAMMA <- PZ%*% X
    Xtilde <- X - u %*%(t(u) %*%X)/ drop(crossprod(u))
    V <- MZ %*%Xtilde
    kappa <- sum(diag(PZ)^2/K)
    tau <- K/N
    
    H <- t(X) %*% PZ %*% X - (ahat) *crossprod(X)
    SIGMAb <- sigma2 *( (1-ahat)^2 *(t(Xtilde) %*% PZ %*% Xtilde) + ahat^2 *( t(Xtilde) %*% MZ %*% Xtilde))
    vcovBekker <- solve(H) %*% SIGMAb %*% solve(H)
    if(vcov=="bekker"){
      V.liml <- vcovBekker
      se.liml <- sqrt(diag(V.liml))
      z.tests <- liml.est/se.liml
      p.val <- 2*pnorm(abs(z.tests), lower=F)
      output <- list(coef.table = cbind(liml.est, se.liml,z.tests,p.val),
                     VCOV = V.liml,
                     k=k)
    }else{
      A <- colSums(diag(PZ-tau)*GAMMA) %*% t(colSums(drop(u^2) * V/N))
      B <- K*(kappa-tau)* ( t((drop(u^2)-sigma2)*V ) %*% V)/(N*(1-2*tau+kappa*tau))
      VCOV.CSE <- solve(H) %*% (SIGMAb + A + t(A) +B) %*% solve(H)
      V.liml <- VCOV.CSE
      se.liml <- sqrt(diag(V.liml))
      z.tests <- liml.est/se.liml
      p.val <- 2*pnorm(abs(z.tests), lower=F)
      output <- list(coef.table = cbind(liml.est, se.liml,z.tests,p.val),
                     VCOV = V.liml,
                     k=k)
    }
  }
  if(vcov=="iid"){
      e.hat <- drop(y-X %*% liml.est)
      sigma <- drop(sqrt(crossprod(e.hat)/(nrow(X)-G)))
      V.liml <- sigma^2 * X_kM_X_inv
      se.liml <- sqrt(diag(V.liml))
      z.tests <- liml.est/se.liml
      p.val <- 2*pnorm(abs(z.tests), lower=F)
      output <- list(coef.table = cbind(liml.est, se.liml,z.tests,p.val),
                     VCOV = V.liml,
                     k=k)
      
  
  }
  
  return(output)
}
