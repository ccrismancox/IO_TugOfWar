#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Counterfactual on kappa
#' ---x
rm(list=ls())

###############################################
# SET absolute change

mag0 <- .13
mag1 <- .13


######################################
# packages i need
library(rootSolve)
library(Matrix)
library(xtable)
library(matrixStats)
library(data.table)
library(ggplot2)
library(sandwich)
library(lmtest)
######################################
source("gamma2trans.R")
source("helperFunctions.r")
###############################################
# results
load("Results/firstStageOutput.Rdata")
gammaStar <- coefficients(mod0)
sigmaStar <- summary(mod0)$sigma
rm( mod0, states, Trans)


load("Results/mainModel.rdata")
states <- model.main[["states"]]
Trans <- model.main[["Trans"]]
thetaEst <- model.main[["regtable"]]$V1
delta <-  as.numeric(model.main[["params"]]$delta)
vEst <- model.main[["v"]]
Obs <-  model.main[["regData"]]

percentages <- c((thetaEst[3] - mag0)/thetaEst[3] -1 ,
                 (thetaEst[4] - mag0)/thetaEst[4] -1 ,
                 1- (thetaEst[3] + mag1)/thetaEst[3]  ,
                 1- (thetaEst[4] + mag1)/thetaEst[4] )*100
                
## mentioned in paper
cat("This magnitude represents a \n",
    round(percentages[1]), "% increase in Hamas' costs\n",
    round(percentages[2]), "% increase in Fatah's costs\n",
    round(percentages[3]),"% decrease in Hamas' costs\n",
    round(percentages[4]), "% decrease in Fatah's costs\n")

mainData <- model.main$regData
mainData$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
mainData$states.int <- sapply(mainData$states.discrete, function(x){which(x==states)})
Tperiod <- dim(mainData)[1]

###############################################
# create states
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()
nG <- dim(G)[1]

# does everything work?
thetaStar = list(betaH = thetaEst[1], # hamas state payoff
                 betaF = thetaEst[2], # fatah state payoff
                 kappaH = c(thetaEst[3], 0), # hamas cost of attack
                 kappaF = c(thetaEst[4], 0), # fatah cost of attack
                 delta=delta)


max(abs(PSI(vEst, thetaStar, Trans, G) - vEst)) < 1e-8
max(abs(Trans - gamma2trans(gammaStar, sigmaStar, states, discretize=F, d=.05, bound=0.025))) < 1e-10
max(abs(PSI(vEst, thetaStar, gamma2trans(gammaStar, sigmaStar, states, discretize=F, d=.05, bound=0.025), G) - vEst)) < 1e-8

EQCP <- AttackProbs(vEst)
EQIV <- invarDist(vEst, Trans)
vH <- log(colSums(exp(matrix(vEst[1:(2*length(states))],  ncol=length(states)))))
vF <- log(colSums(exp(matrix(vEst[(2*length(states)+1):(4*length(states))],  ncol=length(states)))))

EQall <- data.frame(base = c(EQCP$prAH, EQCP$prAF),
                    state = states,
                    actor = rep(c("Hamas", "Fatah"), each = length(states)))


###############################################
# Counterfact: Hamas costs change

steps0 <- seq(from=thetaEst[3], to=thetaEst[3] - mag0, length.out=10)
vCF0 <- matrix(NA, nrow=length(vEst), ncol=length(steps0))
vCF0[,1] <- vEst
theta <- thetaStar


for (i in 2:length(steps0)){
  JV <- PsiDer(vCF0[,i-1], theta, Trans, G)-Diagonal(length(vCF0[,i-1]))
  JKH <- Matrix(c(rep(c(0,1), length(states)), rep(0, length(states)*2)), ncol=1)
  TJ <-  solve(-JV, JKH)
  predicted <- vCF0[,i-1] + as.numeric((steps0[i] - steps0[i-1]) * TJ)

  theta$kappaH[1] <- steps0[i]
  EQ <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, predicted,
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))})
 
  vCF0[,i] <- EQ$root
}

steps1 <- seq(from=thetaEst[3], to=thetaEst[3] + mag1, length.out=10)
vCF1 <- matrix(NA, nrow=length(vEst), ncol=length(steps1))
vCF1[,1] <- vEst
theta <- thetaStar

for (i in 2:length(steps1)){
  JV <- PsiDer(vCF1[,i-1], theta, Trans, G)-Diagonal(length(vCF1[,i-1]))
  JKH <- matrix(c(rep(c(0,1), length(states)), rep(0, length(states)*2)), ncol=1)
  TJ <-  solve(-JV, JKH)
  predicted <- vCF1[,i-1] + as.numeric((steps1[i] - steps1[i-1]) * TJ)

  theta$kappaH[1] <- steps1[i]
  EQ <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, predicted,
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))})
 

  vCF1[,i] <- EQ$root
}

CFCP0 <- AttackProbs(vCF0[,length(steps0)])
CFCP1 <- AttackProbs(vCF1[,length(steps1)])

EQall$H1 <- c(CFCP1$prAH, CFCP1$prAF)
EQall$H0 <- c(CFCP0$prAH, CFCP0$prAF)



###############################################
# Counterfact: Fatah has different costs

steps0 <- seq(from=thetaEst[4], to= thetaEst[4] - mag0, length.out=10)
vCF0 <- matrix(NA, nrow=length(vEst), ncol=length(steps0))
vCF0[,1] <- vEst
theta <- thetaStar

for (i in 2:length(steps0)){
  JV <- PsiDer(vCF0[,i-1], theta, Trans, G)-Diagonal(length(vCF0[,i-1]))
  JKF <- matrix(c(rep(0, length(states)*2), rep(c(0,1), length(states))), ncol=1)
  JVinv <-  solve(-JV, JKF)
  predicted <- vCF0[,i-1] + as.numeric((steps0[i] - steps0[i-1]) * (JVinv))

  theta$kappaF[1] <- steps0[i]
  EQ <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, predicted,
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))})
 
  vCF0[,i] <- EQ$root
}

steps1 <- seq(from=thetaEst[4], to= thetaEst[4] + mag1, length.out=10)
vCF1 <- matrix(NA, nrow=length(vEst), ncol=length(steps1))
vCF1[,1] <- vEst
theta <- thetaStar

for (i in 2:length(steps1)){
  JV <- PsiDer(vCF1[,i-1], theta, Trans, G)-diag(length(vCF1[,i-1]))
  JKF <-  matrix(c(rep(0, length(states)*2), rep(c(0,1), length(states))), ncol=1)
  JVinv <-  solve(-JV, JKF)


  predicted <- vCF1[,i-1] +
    as.numeric((steps1[i] - steps1[i-1]) * (JVinv))

  theta$kappaF[1] <- steps1[i]
  EQ <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, predicted,
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))})
 

  vCF1[,i] <- EQ$root
}

CFCP0 <- AttackProbs(vCF0[,length(steps0)])
CFCP1 <- AttackProbs(vCF1[,length(steps1)])

EQall$F1 <- c(CFCP1$prAH, CFCP1$prAF)
EQall$F0 <- c(CFCP0$prAH, CFCP0$prAF)


###############################################
# Counterfacult: Both

steps0F <- seq(from=thetaEst[4], to= thetaEst[4] - mag0, length.out=30)
steps0H <- seq(from=thetaEst[3], to=thetaEst[3] - mag0, length.out=30)

vCF0 <- matrix(NA, nrow=length(vEst), ncol=length(steps0F))
vCF0[,1] <- vEst
theta <- thetaStar

for (i in 2:length(steps0F)){
  JV <- PsiDer(vCF0[,i-1], theta, Trans, G)-diag(length(vCF0[,i-1]))
  JKH <- matrix(c(rep(c(0,1), length(states)), rep(0, length(states)*2)), ncol=1)
  JKF <- matrix(c(rep(0, length(states)*2), rep(c(0,1), length(states))), ncol=1)

  predicted <- vCF0[,i-1] +
    as.numeric((steps0H[i] - steps0H[i-1]) * (solve(-JV, JKH))) +
      as.numeric((steps0F[i] - steps0F[i-1]) * (solve(-JV, JKF)))

  theta$kappaF[1] <- steps0F[i]
  theta$kappaH[1] <- steps0H[i]
  EQ <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, predicted,
                 jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))})
  vCF0[,i] <- EQ$root
}

steps1F <- seq(from=thetaEst[4], to= thetaEst[4] + mag1, length.out=20)
steps1H <- seq(from=thetaEst[3], to= thetaEst[3] + mag1, length.out=20)
vCF1 <- matrix(NA, nrow=length(vEst), ncol=length(steps1F))
vCF1[,1] <- vEst
theta <- thetaStar

for (i in 2:length(steps1F)){
  JV <- PsiDer(vCF1[,i-1], theta, Trans, G)-diag(length(vCF1[,i-1]))
  JKF <-  matrix(c(rep(0, length(states)*2), rep(c(0,1), length(states))), ncol=1)
  JKH <- matrix(c(rep(c(0,1), length(states)), rep(0, length(states)*2)), ncol=1)
  
  predicted <- vCF1[,i-1] +
    as.numeric((steps1H[i] - steps1H[i-1]) * (solve(-JV, JKH))) +
    as.numeric((steps1F[i] - steps1F[i-1]) * (solve(-JV, JKH)))

  theta$kappaF[1] <- steps1F[i]
  theta$kappaH[1] <- steps1H[i]

  EQ <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, predicted,
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){as.matrix(PsiDer(V, theta, Trans, G)-Diagonal(length(V)))})
 

  vCF1[,i] <- EQ$root
}

CFCP0 <- AttackProbs(vCF0[,length(steps0)])
CFCP1 <- AttackProbs(vCF1[,length(steps1)])

EQall$Both1 <- c(CFCP1$prAH, CFCP1$prAF)
EQall$Both0 <- c(CFCP0$prAH, CFCP0$prAF)






###############################################
# Summarize changes

l1 <- c(mean(EQall$base[EQall$actor=="Hamas"]),
    mean(EQall$base[EQall$actor=="Fatah"]),
    mean(EQall$base[EQall$actor=="Hamas"] + EQall$base[EQall$actor=="Fatah"] - (EQall$base[EQall$actor=="Hamas"]*EQall$base[EQall$actor=="Fatah"])))

l2 <- c(mean(EQall$H0[EQall$actor=="Hamas"]),
        mean(EQall$H0[EQall$actor=="Fatah"]),
        mean(EQall$H0[EQall$actor=="Hamas"] + EQall$H0[EQall$actor=="Fatah"] - (EQall$H0[EQall$actor=="Hamas"]*EQall$H0[EQall$actor=="Fatah"])))

l3 <- c(mean(EQall$F0[EQall$actor=="Hamas"]),
        mean(EQall$F0[EQall$actor=="Fatah"]),
        mean(EQall$F0[EQall$actor=="Hamas"] + EQall$F0[EQall$actor=="Fatah"] - (EQall$F0[EQall$actor=="Hamas"]*EQall$F0[EQall$actor=="Fatah"])))

l4 <- c(mean(EQall$Both0[EQall$actor=="Hamas"]),
        mean(EQall$Both0[EQall$actor=="Fatah"]),
        mean(EQall$Both0[EQall$actor=="Hamas"] + EQall$Both0[EQall$actor=="Fatah"] - (EQall$Both0[EQall$actor=="Hamas"]*EQall$Both0[EQall$actor=="Fatah"])))



EQall <- data.table(EQall)
Either <- EQall[actor=="Hamas", !c("state", "actor")] + EQall[actor=="Fatah", !c("state", "actor")]-
    EQall[actor=="Hamas", !c("state", "actor")]* EQall[actor=="Fatah", !c("state", "actor")]


Either[, `:=`(actor= "Either", state=unique(EQall$state))]
EQall <- merge(EQall, Either, all=TRUE)
X0 <- EQall[, .(Baseline=mean(base), 
         Costly.Hamas = mean(H0), Costly.Fatah=mean(F0), Costly.both=mean(Both0),
         Easier.Hamas = mean(H1), Easier.Fatah = mean(F1),Easier.both= mean(Both1)), 
      by=actor]
X <- t(X0[,!c("actor")])
colnames(X) <- X0$actor
X <- num2str(X[,c("Hamas", "Fatah", "Either")])
colnames(X) <- c("Pr(Hamas attacks)", "Pr(Fatah attacks)", "Pr(Either attack)")
X <- cbind(c("Baseline",
             "Increase costs for", "", "",
             "Decrease costs for", "", ""),
           c("", rep(c("Hamas", "Fatah", "Both"), 2)),
           X)
X <- rbind(colnames(X), X)
rownames(X) <- NULL

cat(kable(X, format="pipe",
          caption="Average attack probabilities as kappa changes"),
    file="../../Output/Tables/Table5.txt",
    sep="\n")
