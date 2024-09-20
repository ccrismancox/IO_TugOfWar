#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Numerical example searching for equilibria
#' ---

#### DO NOT RUN UNLESS NEEDED

rm(list=ls())

######################################
# packages i need
library("foreach")
library("doParallel")
library("doRNG")
library("rootSolve")
library("Matrix")
source("helperFunctions.r")
source("gamma2trans.R")
###############################################
# set up model
states <- seq(from = -50, to = 50, by = 1) 
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()
nG <- dim(G)[1]

thetaStar = list(betaH = -1/500, # hamas state payoff 
                 betaF = 1/500, # fatah state payoff
                 kappaH = c(-2, 0), # hamas cost of attack
                 kappaF = c(-2, 0), # fatah cost of attack
                 delta=0.999)  

gammaStar <- c(0, -1, 1, 1, 0, 0)
sigmaStar <- 2

# transition matrix
Trans <- gamma2trans(gammaStar, sigmaStar, states, d=1, discretize=F)

# per period utilities
uas <- cbind(thetaStar$betaH*G$state + (thetaStar$kappaH[1] + thetaStar$kappaH[2]*G$state)*G$aH,
             thetaStar$betaF*G$state + (thetaStar$kappaF[1] + thetaStar$kappaF[2]*G$state)*G$aF)
colnames(uas) <- c("uH", "uF")

###############################################
# the symmetric eq
vstart <- c(rowMeans(matrix(uas[,'uH'], ncol=2, byrow=T)),
            c(apply(array(uas[,'uF'], c(2,2,length(states))), 3, rowMeans)))

thetaStar$delta <- 0.90
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, vstart, 
                   maxiter = 100, jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.99
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 100, jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.9945
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 200,  jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.99675
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 200,  jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.999
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 200,  jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})



###############################################
# other eq

ng <- 100000
np <- length(EQsym$root)
set.seed(1)

cl <- makeCluster(min(detectCores() - 1, 18))
registerDoParallel(cl)


sols <- foreach(x=1:ng, .combine = 'cbind', .packages = 'rootSolve', .inorder=F) %dorng% {
  sv <- runif(np, min = 0, max=200)
  use <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, sv, 
                   maxiter = 250, rtol=1e-10, atol=1e-12, ctol=1e-12,
                   jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})
  if (max(abs(use$f.root)) < 1e-10){
    out <- AttackProbs(use$root)
  } else {
    out <- rep(NA, np/2)
  }
  out
}
stopCluster(cl)
sols <- sols[, !is.na(sols[1,])]
solsU <- unique(round(sols,digits=6), MARGIN=2)
dim(solsU)
