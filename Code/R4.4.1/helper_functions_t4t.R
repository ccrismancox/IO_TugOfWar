PSI <-  function(V, theta, Trans, G){
  nS <- length(V)/4
  # Values, vi[j,k] := actor i, action j, state k
  # note: row 2 is payoff from attack.  Row 1 is payoff from do nothing
  vH <- matrix(V[1:(2*nS)],  ncol=nS) 
  vF <- matrix(V[(2*nS+1):(4*nS)],  ncol=nS)
  
  # # normalize
  # note: given our parameterization, we know the answer to the max step
  # mvH should equal vH[1,]
  mvH <- colMaxs(vH)
  mvF <- colMaxs(vF)
  evHn <- exp(sweep(vH, 2, mvH, "-"))
  evFn <- exp(sweep(vF, 2, mvF, "-"))
  
  # Condition choice probabilties
  PH <- matrix(sweep(evHn, 2, colSums(evHn), "/"), ncol=1)
  PF <- matrix(sweep(evFn, 2, colSums(evFn), "/"), ncol=1)
  
  # U_i(a|s)
  GF <- G$state >=3
  GH <- G$state == 2 | G$state==4
  
  uas <- cbind(theta$tauH*GH*G$aH + theta$kappaH*G$aH,
               theta$tauF*GF*G$aF + theta$kappaF*G$aF)
  colnames(uas) <- c("uH", "uF")
  UHsa <- uas[,'uH'] + theta$delta * (Trans %*% (log(colSums(evHn)) + mvH))
  UFsa <- uas[,'uF'] + theta$delta * (Trans %*% (log(colSums(evFn)) + mvF))
  
  # Multiply by p_{-i}( a_{-1}|s)
  # verify this with PF and PH equal to 1:2*S
  VH <- matrix(matrix(PF, ncol=nS) %x% matrix(c(1,1), nrow=1), ncol=1)*UHsa
  VF <- matrix(matrix(PH, ncol=nS) %x% matrix(c(1,1), ncol=1), ncol=1)*UFsa
  
  # Sum over a_{-i}
  VH <- rowSums(matrix(VH, ncol=2,byrow=T))
  VF <- matrix(apply(array(VF, c(2,2,nS)), 3, rowSums), ncol=1)
  return(c(VH,VF))
}

PsiDer <- function(V, theta, Trans, G){
  # Jacbian of F(v) = Psi(v;theta), with respect to v ****NOT Psi(v)-v****
  #
  #nAct1 <- max(given$stateActions[,1])
  #nAct2 <- max(given$stateActions[,2])
  #nActs <- nAct1*nAct2
  #nS <- length(given$statesInData)/given$mun
  nAct1 <- nAct2 <- 2
  nActs <- nAct1*nAct2
  nS <- length(V)/4
  
  vH <- matrix(V[1:(2*nS)],  ncol=nS)
  vF <- matrix(V[(2*nS+1):(4*nS)],  ncol=nS)
  
  
  
  ##normalization
  mvH <- colMaxs(vH)
  mvF <- colMaxs(vF)
  evHn <- exp(sweep(vH, 2, mvH, "-"))
  evFn <- exp(sweep(vF, 2, mvF, "-"))
  
  
  
  GH <- mvH + log(colSums(evHn))
  GF <- mvF + log(colSums(evFn))
  
  #probs
  pH <- sweep(evHn, 2, colSums(evHn), FUN="/")
  pF <- sweep(evFn, 2, colSums(evFn), FUN="/")
  
  
  # U_i(a|s)
  G$GF <- G$state >=3
  G$GH <- G$state == 2 | G$state==4
  
  uas <- cbind(theta$tauH*G$GH*G$aH + theta$kappaH*G$aH,
               theta$tauF*G$GF*G$aF + theta$kappaF*G$aF)
  colnames(uas) <- c("uH", "uF")
  UHsa <- uas[,'uH'] + theta$delta * (Trans %*% (log(colSums(evHn)) + mvH))
  UFsa <- uas[,'uF'] + theta$delta * (Trans %*% (log(colSums(evFn)) + mvF))
  
  
  #Set up Jacobian of Phi
  #JPhi <- Matrix(0, nrow=nS*(nAct1+nAct2), ncol=nS*(nAct1+nAct2))
  pF.1 <- t(pF) %x% rep(1, nAct1)
  pF.2 <- matrix(t(t(pF) %x% rep(1, nAct2)), ncol=1)
  P <- new("dgCMatrix",
           Dim=c(nrow(pF.1), nrow(pF.1)*2L),
           i = rep(0:(nrow(pF.1)-1), each=2),
           p= 0:(nrow(pF.1)*2L),
           x= c(t(pF.1)))
  Q <- Trans
  idx <- as.numeric(rep(1,2) %x% matrix(rep(1:(2*nS)) ,2))
  
  ##dV1/dV1
  JPhi0 <- (theta$delta*P%*% Q %x% t(rep(1,nAct1))) %*% Diagonal(nS*nAct1,as.numeric(pH))
  #JPhi[1:(nS*nAct1),1:(nS*nAct1)] <- JPhi0
  
  # max(abs(JPhi[1:(nS*nAct1),1:(nS*nAct1)] - J[1:(nS*nAct1),1:(nS*nAct1)] ))
  #dV1/dV2
  jsp <- (1:(nS*nAct2))%x%rep(1,nAct1)
  isp <- as.numeric(matrix(1:(nS*nAct1), nrow=nAct1,ncol=nS) %x% t(rep(1,nAct2)))
  hsp <- pF.2* (UHsa - P %*% UHsa %x% rep(1,2))
  JPhi1 <- sparseMatrix(jsp,isp,x=drop(hsp),dims=c(nS*nAct1,nS*nAct2))
  #JPhi[1:(nS*nAct1),(nS*nAct1+1):(nS*nAct1+nS*nAct2)] <- JPhi1
  # max(abs(JPhi[1:(nS*nAct1),(nS*nAct1+1):(nS*nAct1+nS*nAct2)]  - J[1:(nS*nAct1),(nS*nAct1+1):(nS*nAct1+nS*nAct2)]  ))
  
  #dV2/dV2
  jsp <- rep(1:(nS*nActs), nAct2)
  isp <- matrix( matrix(1:(nS*nAct2),nS,nAct2, byrow=TRUE), ncol=1) %x% rep(1, nActs)
  hsp <- Diagonal(nAct2) %x% pH
  dim(hsp) <- c(nS*nAct2*nActs,1)
  P <- sparseMatrix(isp,jsp,x=drop(hsp),dims=c(nS*nAct2,nS*nActs))
  JPhi2 <-  ((theta$delta*P %*% Q[order(idx),] %x% t(rep(1,nAct2))) %*% Diagonal(nS*nAct2, matrix(pF,nS*nAct2,1)))
  #JPhi[(nS*nAct1+1):(nS*nAct1+nS*nAct2),(nS*nAct1+1):(nS*nAct1+nS*nAct2)] <- JPhi2
  
  # max(abs(JPhi[(nS*nAct1+1):(nS*nAct1+nS*nAct2),(nS*nAct1+1):(nS*nAct1+nS*nAct2)] -
  # J[(nS*nAct1+1):(nS*nAct1+nS*nAct2),(nS*nAct1+1):(nS*nAct1+nS*nAct2)] ))
  
  ## dV2/dV1
  jsp <- as.numeric( (matrix(1:(nS*nAct1),nAct1,nS) %x% t(rep(1,nAct2))))
  isp <-  rep(1:(nS*nAct2), each=nAct1)
  x2 <- UFsa[order(idx)] - P %*% UFsa[order(idx)]%x% rep(1,2)
  pH.F <- matrix(t(t(pH) %x% rep(1, 2)), ncol=1)
  hsp <- pH.F*x2
  JPhi3 <- sparseMatrix(isp,jsp, x=drop(hsp), dims= c(nS*nAct2,nS*nAct1))
  
  ## JPhi[(nS*nAct1+1):(nS*nAct1+nS*nAct2),1:(nS*nAct1)] <- Jphi3
  # max(abs(JPhi[(nS*nAct1+1):(nS*nAct1+nS*nAct2),1:(nS*nAct1)]- J[(nS*nAct1+1):(nS*nAct1+nS*nAct2),1:(nS*nAct1)] ))
  JPhi <-  rbind(cbind(JPhi0, JPhi1), cbind(JPhi3, JPhi2))
  
  
  return(as.matrix(JPhi))
} 
