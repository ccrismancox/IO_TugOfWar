#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Second stage robustness to time windows
#' output: pdf_document
#' ---

library(data.table)
library(stringr)
library(knitr)
library(numDeriv)
library(Matrix)
library(zoo)
rm(list=ls())
load("Results/firstStageOutput.Rdata")
load("Results/mainStateActions.Rdata")
load("Results/mainModel.rdata")
load("Results/startingvalues.rdata")
start0 <- c(startmod$regtable2$V1, startmod$v)
load("Results/TimePeriods.rdata")
source("helperFunctions.r")

system("mkdir -p ipoptTEMP")

# Common outputs
py.params <- data.table(delta=0.999,
                        nkappa=2, nbeta=2,
                        twostep=FALSE)
write.csv(Trans, file=paste("ipoptTEMP/trans.csv", sep=""),row.names = F)
write.csv(states, file=paste("ipoptTEMP/statespace.csv", sep=""),row.names = F)
write.csv(py.params, file=paste("ipoptTEMP/params.csv", sep=""),row.names = F)
write.csv(c(mod0$coef, summary(mod0)$sigma), file="ipoptTEMP/gamma.csv",row.names = F)
write.csv(bootOut, file="ipoptTEMP/V1.csv",row.names = F)

model.main$Dates <- mainData$Date



######limit the time#######
## 1.  Reconciliation deal April 2014
regData <- mainData[Date < "2014/04/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)]
Date <- regData$Date
regData$Date <- NULL

start <- try(with(shortT.main_2014, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}
write.csv(start, file=paste("ipoptTEMP/start.csv", sep=""),row.names = F)
write.csv(regData, file=paste("ipoptTEMP/regData.csv", sep=""),row.names = F)

cat("IPOPT running, please wait"\n)
system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
cat("IPOPT complete\n")

regtable <- read.csv(paste("ipoptTEMP/regtable.csv", sep=""),header = F)
conv <- read.csv(paste("ipoptTEMP/convergence.csv", sep=""),header = F)
v <- read.csv(paste("ipoptTEMP/v.csv", sep=""),header = F)[,1]
shortT.main_2014 <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                         conv=conv, v=v, regtable=regtable,
                         Dates=Date)

save(list=c("model.main",  "shortT.main_2014"), file="Results/TimePeriods.rdata")




## 2. Reconciliation deal signed in May 2011
regData <- mainData[Date < "2011/05/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)] 
Date <- regData$Date
regData$Date <- NULL

start <- try(with(shortT.main_2011, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}
write.csv(start, file=paste("ipoptTEMP/start.csv", sep=""),row.names = F)
write.csv(regData, file=paste("ipoptTEMP/regData.csv", sep=""),row.names = F)

cat("IPOPT running, please wait"\n)
system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
cat("IPOPT complete\n")


regtable <- read.csv(paste("ipoptTEMP/regtable.csv", sep=""),header = F)
conv <- read.csv(paste("ipoptTEMP/convergence.csv", sep=""),header = F)
v <- read.csv(paste("ipoptTEMP/v.csv", sep=""),header = F)[,1]
shortT.main_2011 <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                                 conv=conv, v=v, regtable=regtable, 
                         Dates=Date)

save(list=c("model.main", "shortT.main_2011",  "shortT.main_2014"), file="Results/TimePeriods.rdata")




## 3.  Last Fatah attack Mar. 2009
regData <- mainData[Date < "2009/04/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)]
Date <- regData$Date
regData$Date <- NULL

start <- try(with(shortT.main_2009, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}
write.csv(start, file=paste("ipoptTEMP/start.csv", sep=""),row.names = F)
write.csv(regData, file=paste("ipoptTEMP/regData.csv", sep=""),row.names = F)

cat("IPOPT running, please wait"\n)
system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
cat("IPOPT complete\n")

regtable <- read.csv(paste("ipoptTEMP/regtable.csv", sep=""),header = F)
conv <- read.csv(paste("ipoptTEMP/convergence.csv", sep=""),header = F)
v <- read.csv(paste("ipoptTEMP/v.csv", sep=""),header = F)[,1]
shortT.main_2009 <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                         conv=conv, v=v, regtable=regtable, 
                         Dates=Date)

save(list=c("model.main", "shortT.main_2011", "shortT.main_2014", "shortT.main_2009"), 
     file="Results/TimePeriods.rdata")



## 4. elections in Jan. 2006
regData <- mainData[Date < "2006/01/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)]
Date <- regData$Date
regData$Date <- NULL

start <- try(with(shortT.main_2006, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}
write.csv(start, file=paste("ipoptTEMP/start.csv", sep=""),row.names = F)
write.csv(regData, file=paste("ipoptTEMP/regData.csv", sep=""),row.names = F)

cat("IPOPT running, please wait"\n)
system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
cat("IPOPT complete\n")


regtable <- read.csv(paste("ipoptTEMP/regtable.csv", sep=""),header = F)
conv <- read.csv(paste("ipoptTEMP/convergence.csv", sep=""),header = F)
v <- read.csv(paste("ipoptTEMP/v.csv", sep=""),header = F)[,1]
shortT.main_2006 <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                         conv=conv, v=v, regtable=regtable, 
                         Dates=Date)


save(list=c("model.main", "shortT.main_2011", "shortT.main_2014","shortT.main_2009", "shortT.main_2006"),
     file="Results/TimePeriods.rdata")



## 5. 2nd intifada starts in Sept 2000
regData <- mainData[Date < "2000/9/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)] #2011 was OK, try 2014?
Date <- regData$Date
regData$Date <- NULL
start <- try(with(shortT.main_2000, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}

### CMLE failed try two step
vhat <- start[-c(1:4)]



# create states
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()
nG <- dim(G)[1]

given <- list(G=G, P=Trans)

Y <- regData[,.(states.discrete, Hattacks, Fattacks)]
Y$states.int <- sapply(Y$states.discrete, \(x){which(states==x)})

statesInData <- table(factor(Y$states.discrete, levels=states))

stateActions <- expand.grid(states,0:1)
stateActions <- stateActions[order(stateActions[,1]),]
stateActionsInData <- matrix(0, nrow=dim(stateActions)[1]*2, ncol=1)
### ^ *2 for number of players
addon <- nrow(stateActions)
for (i in 1:dim(stateActions)[1]){
  stateActionsInData[i] <- sum(Y$states.discrete==stateActions[i,1] & Y$Hattacks==stateActions[i,2])
  stateActionsInData[i+addon] <- sum(Y$states.discrete==stateActions[i,1] & Y$Fattacks==stateActions[i,2])
}

given$statesInData <- statesInData
given$stateActions <- stateActions
given$stateActionsInData <- stateActionsInData
HM <- function(theta, vhat, Y, given){
    thetaStar = list(betaH = theta[1], # hamas state payoff
                     betaF = theta[2], # fatah state payoff
                     kappaH = c(theta[3], 0), # hamas cost of attack
                     kappaF = c(theta[4], 0), # fatah cost of attack
                     delta=.999)
 
  Vout <- PSI(V=vhat, theta=thetaStar, Trans=given$P, G=given$G)
  Pout <- AttackProbs(Vout)

  S = ncol(given$P)

  VH = Vout[1:(2*S)]
  VF = Vout[(2*S+1):(4*S)]
  vHm <- matrix(VH, ncol=S)
  vFm <- matrix(VF, ncol=S)

  # # normalize
  mvH <- apply(vHm, 2, max)
  mvF <- apply(vFm, 2, max)
  evHn <- exp(sweep(vHm, 2, mvH, "-"))
  evFn <- exp(sweep(vFm, 2, mvF, "-"))


  cSAD1 = matrix(given$stateActionsInData[1:(S*2)], ncol=1)
  cSAD2 = matrix(given$stateActionsInData[(S*2+1):(S*4)], ncol=1)

  LL = VH %*% cSAD1 -  sum((log(colSums(evHn))+mvH) * given$statesInData)
  LL = LL + VF %*% cSAD2  - sum((log(colSums(evFn))+mvF) * given$statesInData)
  return(-LL)

}



given$P <- Trans
tol <- 100
while(tol > 1e-5){
  fun <- function(x){HM(theta=x, vhat=vhat, given=given, Y)}
  twostep <- optim(par=start[1:4], fn = fun, method="BFGS")
  theta <- twostep$par
  thetaStar = list(betaH = theta[1], # hamas state payoff
                   betaF = theta[2], # fatah state payoff
                   kappaH = c(theta[3], 0), # hamas cost of attack
                   kappaF = c(theta[4], 0), # fatah cost of attack
                   delta=.999)
  vnew <- PSI(vhat, theta = thetaStar,G = given$G, Trans = given$P)
  tol <- max(abs(c(twostep$par,vnew)-start))
  start <- c(twostep$par,vnew)
  vhat <- vnew
  cat("tol: ", tol, "\n")
}

conv <- data.frame(V1=c(twostep$convergence, twostep$value))
v <- vhat
regtable <- data.frame(V1 = twostep$par)



HM_byobs <- function(theta, vhat, given,Y){
  thetaStar = list(betaH = theta[1], # hamas state payoff
                   betaF = theta[2], # fatah state payoff
                   kappaH = c(theta[3], 0), # hamas cost of attack
                   kappaF = c(theta[4], 0), # fatah cost of attack
                   delta=.999)
  

  Vout <- PSI(V=vhat, theta=thetaStar, Trans=given$P, G=given$G)
  Pout <- AttackProbs(Vout)


  S = ncol(given$P)

  VH = Vout[1:(2*S)]
  VF = Vout[(2*S+1):(4*S)]
  vHm <- matrix(VH, ncol=S)
  vFm <- matrix(VF, ncol=S)

    Vidx <-expand.grid(actions=c(0,1),given$states)[,c(2,1)]
  ll <- rep(0, nrow(Y))
  for(i in 1:nrow(Y)){
    col.idx <- which(given$states==Y$states.discrete[i])
    Hrow.idx <- Y$Hattacks[i]+1
    Frow.idx <- Y$Fattacks[i]+1
    lPH <- -log(1+exp(vHm[-Hrow.idx, col.idx]-vHm[Hrow.idx, col.idx]))
    lPF <- -log(1+exp(vFm[-Frow.idx, col.idx]-vFm[Frow.idx, col.idx]))
    ll[i] <- lPH+lPF
  }


  return(-ll)

}

given$states <- states
thetaStar = list(betaH = twostep$par[1], # hamas state payoff
                 betaF = twostep$par[2], # fatah state payoff
                 kappaH = c(twostep$par[3], 0), # hamas cost of attack
                 kappaF = c(twostep$par[4], 0), # fatah cost of attack
                 delta=.999)
psi_v <- PsiDer(V=vhat, theta=thetaStar, Trans=given$P, G=given$G)
psiDertheta <- jacobian(func=function(theta){
  thetaStar = list(betaH = theta[1], # hamas state payoff
                   betaF = theta[2], # fatah state payoff
                   kappaH = c(theta[3], 0), # hamas cost of attack
                   kappaF = c(theta[4], 0), # fatah cost of attack
                   delta=.999)
  return(PSI(V=vhat, theta=thetaStar, Trans=given$P, G=given$G))
  },
  twostep$par
  )

Jtheta <- jacobian(twostep$par, func=HM_byobs,vhat = vhat,given = given,Y=Y)
Jv <- jacobian(vhat, func=HM_byobs,theta=twostep$par, given = given,Y=Y)
S <- nrow(Trans)
OMEGA_theta <-  crossprod(Jtheta)
OMEGA_v <- crossprod(Jtheta,Jv)
left <- solve(OMEGA_theta + OMEGA_v %*% solve(diag(S)-t(psi_v)) %*%psiDertheta)
right <- solve(OMEGA_theta + t(psiDertheta) %*% solve(diag(S)-(psi_v)) %*%t(OMEGA_v))
VCOV_2000 <- left %*% OMEGA_theta %*% right
se <- sqrt(diag(VCOV_2000))

regtable$V2 <- se
regtable$V3 <- regtable$V1/regtable$V2
regtable$V4 <- 2*pnorm(abs(regtable$V3), lower=F)
shortT.main_2000 <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                         conv=conv, v=v, regtable=regtable,
                         Dates=Date)


save(list=c("model.main", "shortT.main_2011",
            "shortT.main_2014", "shortT.main_2006", "shortT.main_2009",
            "shortT.main_2000"), file="Results/TimePeriods.rdata")

## 6.  Start the same time as Bloom
regData <- mainData[Date >= "1997/04/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)] 
Date <- regData$Date
regData$Date <- NULL

start <- try(with(shortT.main_bloom, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}
write.csv(start, file=paste("ipoptTEMP/start.csv", sep=""),row.names = F)
write.csv(regData, file=paste("ipoptTEMP/regData.csv", sep=""),row.names = F)

cat("IPOPT running, please wait"\n)
system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
cat("IPOPT complete\n")


regtable <- read.csv(paste("ipoptTEMP/regtable.csv", sep=""),header = F)
conv <- read.csv(paste("ipoptTEMP/convergence.csv", sep=""),header = F)
v <- read.csv(paste("ipoptTEMP/v.csv", sep=""),header = F)[,1]
shortT.main_bloom <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                         conv=conv, v=v, regtable=regtable, 
                         Dates=Date)

save(list=c("model.main", "shortT.main_2011", "shortT.main_2014",
            "shortT.main_2009", "shortT.main_2006",
            "shortT.main_2000",
            "shortT.main_bloom"),
     file="Results/TimePeriods.rdata")


## 7. start from 2001
regData <- mainData[Date >= "2001/01/01",list(Date,states, Hattacks, Fattacks, lag.states, states.discrete)] 
Date <- regData$Date
regData$Date <- NULL
start <- try(with(shortT.main_from2001, c(regtable[,1], v)), silent = TRUE)
if("try-error" %in% class(start)){
  start <- start0
}
write.csv(start, file=paste("ipoptTEMP/start.csv", sep=""),row.names = F)
write.csv(regData, file=paste("ipoptTEMP/regData.csv", sep=""),row.names = F)

cat("IPOPT running, please wait"\n)
system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
cat("IPOPT complete\n")


regtable <- read.csv(paste("ipoptTEMP/regtable.csv", sep=""),header = F)
conv <- read.csv(paste("ipoptTEMP/convergence.csv", sep=""),header = F)
v <- read.csv(paste("ipoptTEMP/v.csv", sep=""),header = F)[,1]
shortT.main_from2001 <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                         conv=conv, v=v, regtable=regtable, 
                         Dates=Date)

save(list=c("model.main", "shortT.main_2011", "shortT.main_2014",
            "shortT.main_2009", "shortT.main_2006",
            "shortT.main_2000",
            "shortT.main_bloom",
            "shortT.main_from2001"),
     file="Results/TimePeriods.rdata")
system(paste("rm ipoptTEMP/*.csv", sep=""))

system("rm ipoptTEMP -r")


############## Make a table ############
model.main$Dates <- seq.Date(as.Date("1994-01-01"), as.Date("2018-12-01"), "months")
models <- list(`Full Sample`=model.main,
               `2014 Agreement`=shortT.main_2014,
               `2011 Agreement`=shortT.main_2011,
               `2009 Last Fatah attack`=shortT.main_2009,
               `2006 Elections`=shortT.main_2006,
               `Second intifada`=shortT.main_2000,
               `Bloom's (2004) start year`=shortT.main_bloom,
               `Start 2001`=shortT.main_from2001)

tableE1 <- do.call(cbind,lapply(models,
       function(x){
         c(num2str(matrix(t(x$regtable[,1:2]),ncol=1)), 
           nrow(x$regData), 
                  paste(as.yearmon(c(min(x$Dates), 
                                      max(x$Dates)))),
           num2str(-1*x$conv$V1[2]))
       }))
tableE1[c(2,4,6,8), ] <- paste0("(",tableE1[c(2,4,6,8), ],")")
tableE1 <- cbind(" " =c("beta_H", "",
                   "beta_F", "",
                   "kappa_H", "",
                   "kappa_F", "",
                   "T", 
                   "Start date",
                   "End date",
                   "LL"),
                 tableE1)
tableE1[12,7]  <- paste0(tableE1[12,7], "^dagger")

cat(kable(tableE1, format="pipe",
          caption="Robustness to different time periods"),
    file="../../Output/Tables/tableH1.txt", sep="\n")
cat("Note: ^dagger CMLE did not converge, estimates from nested-pseudo-likelihood (NPL) estimator.  BH standard errors in parentheses.\n", file="../../Output/Tables/tableH1.txt", append=TRUE)


