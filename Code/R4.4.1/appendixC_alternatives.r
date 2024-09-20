#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Measurement model robustness
#' ---
#' 
#' Clear workspace and load packages:
rm(list=ls())

library(MARSS)
library(data.table)
library(ggplot2)
library(knitr)
library(zoo) #for working with year-months
load("../../Data/actionsSetup.Rdata")


jmccWB <- subset(read.csv("../../Data/jmcc_WB.csv"), year <= 2018)
jmccGZ <- subset(read.csv("../../Data/jmcc_GAZA.csv"), year <= 2018)

cpsrWB <- subset(read.csv("../../Data/cpsr_WB.csv"), year <= 2018)
cpsrGZ <- subset(read.csv("../../Data/cpsr_GAZA.csv"), year <= 2018)


d <- .05
lower <- 0.025
#' Check the factors that measure s

use <- subset(dat, select=c(trustHamas,trustFatah,
                            supportHamas,supportFatah,
                            legisHamas, legisFatah))

setnames(jmccWB,
         c("trust_1", "trust_2", "trust_5", "legis_1", "legis_2", "legis_4"),
         c("trustFatah", "trustHamas", "trustNone",  "legisFatah", "legisHamas", "legisNone"))

setnames(jmccGZ,
         c("trust_1", "trust_2", "trust_5", "legis_1", "legis_2", "legis_4"),
         c("trustFatah", "trustHamas", "trustNone",  "legisFatah", "legisHamas", "legisNone"))

use.gaza <- cbind.data.frame(with(jmccGZ, cbind(trustHamas,trustFatah)),
                             with(cpsrGZ, cbind(supportHamas,supportFatah)),
                             with(jmccGZ, cbind(legisHamas,legisFatah)))
use.WB <- cbind.data.frame(with(jmccWB, cbind(trustHamas,trustFatah)),
                           with(cpsrWB, cbind(supportHamas,supportFatah)),
                           with(jmccWB, cbind(legisHamas,legisFatah)))

X1 <- data.frame(Constant=1,
                 HA=as.numeric(dat$lag.Hattacks>0), 
                 FA = as.numeric(dat$lag.Fattacks>0))


## use B is the correlation AR(1) term
## use R is the variance  of xi
## use C is alpha 
specifications <- data.frame(useB=c("unconstrained", "identity", "unconstrained", "unconstrained","unconstrained"),
                             useR=c("identity", "identity", "diagonal and equal", "diagonal and unequal","identity"),
                             useC=c("unconstrained","unconstrained","unconstrained","unconstrained","zero"))
output <- matrix(0, nrow=300, ncol=8)
for(i in 1:8){
  j <- ifelse(i >= 6, 1, i)
  if(i >= 6){
    if(i==6){USE <- use[,1:4]} #remove plan to vote
    if(i==7){USE <- use.gaza} # use Gaza only
    if(i==8){USE <- use.WB} #use WB only
  } else{
    USE <- use
  }
  mod <-  list(Z=matrix(names(USE)),
               A="zero", 
               R=specifications$useR[j],
               B=specifications$useB[j],
               U="zero",
               Q="identity", 
               x0="zero",
               C=specifications$useC[j],
               c = t(X1),
               V0=matrix(5,1,1),
               tinitx=1)
  cntl.list = list(maxit=10000, minit=500, abstol = 1e-5, conv.test.slope.tol = 0.05)
  #' Fit the measurement model
  
  
  
  mars <- MARSS(t(scale(USE)), model=mod, control=cntl.list, silent=T)
  
  
  regData <- data.table(date=dat$date, 
                        states = mars$states[1,],
                        lag.Hattacks = mod$c["HA",],
                        lag.Fattacks = mod$c["FA",],
                        Hkills = dat$Hkills, Fkills=dat$Fkills)
  #' MARSS latent variable is not sign identified so let's make sure that it moves 
  #' in the way it should. Observation 153 is the 2006
  if(sign(regData$states[153])==1){
    d <-d
    
    bounds <- quantile( -mars$states[1,], c(lower, 1-lower))
    states <- seq(from=bounds[1], to=bounds[2], by=d)
    
    G <- expand.grid(states, c(0,1), c(0,1))
    names(G) = c("state", "aH", "aF")
    G <- G[order(G$state, G$aH),]
    rownames(G) <- c()
    nG <- dim(G)[1]
    regData$states <- -1* regData$states
    mars$par$Z <- -1*mars$par$Z
    mars$par$U <- -1*mars$par$U
  }else{
    d <- d
    
    bounds <- quantile( mars$states[1,], c(lower, 1-lower))
    states <- seq(from=bounds[1], to=bounds[2], by=d)
    
    G <- expand.grid(states, c(0,1), c(0,1))
    names(G) = c("state", "aH", "aF")
    G <- G[order(G$state, G$aH),]
    rownames(G) <- c()
    nG <- dim(G)[1]
  }
  output[,i] <- regData$states
}
X <- cor(output)
rownames(X) <- colnames(X) <- paste("Model", 1:8)


cat(kable(X, digits=2,
          caption="Correlations across measurement model specifications"),
    file="../../Output/Tables/TableC4.txt",
    sep="\n")
          
