#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Sensitivity of the second stage to changes in the first 
#' ---

library(foreach)
library(doParallel)
library(doRNG)
library(data.table)
library(MASS)
rm(list=ls())
source("helperFunctions.r")
source("gamma2trans.R")
source("firststageboot.r")

load("Results/mainStateActions.Rdata")
load("Results/firstStageOutput.Rdata")
load("Results/mainModel.rdata")
load("Results/startingvalues.rdata")

states <- model.main$states
firstStageData <- regData[, list(states, lag.Hattacks, lag.Fattacks,lag.states)]

start <- c(startmod$regtable$V1, startmod$v)
B <- 500 

## parametric bootstrap
workers <- floor(detectCores() *7/8)
cl <- makeCluster(workers)
registerDoParallel(cl)

set.seed(1234)
bootOut <- foreach(b= 1:B, .combine = 'rbind', 
                   .packages = c("data.table"), 
                   .errorhandling = "remove") %dorng%{
                     
                     
  state <- c(firstStageData$states[1], rep(0, 299)) 

  for(t in 2:300){
    state[t] <- rnorm(1,mean=predict(mod0, newdata=data.frame(lag.states=state[t-1],
                                                              lag.Hattacks=firstStageData$lag.Hattacks[t],
                                                              lag.Fattacks=firstStageData$lag.Fattacks[t])),
                      sd=summary(mod0)$sigma)
  }

   
  regData.boot <- data.table(states=state,
                             lag.Hattacks=firstStageData$lag.Hattacks,
                             lag.Fattacks=firstStageData$lag.Fattacks)
  regData.boot[,lag.states := shift(state)]
  mod0a <- lm( states ~ lag.Hattacks + lag.Fattacks + lag.states + 
                 lag.Hattacks:lag.states + lag.Fattacks:lag.states, data=regData.boot)
  
  mainData.boot <- data.table(states=mainData$states,
                              Hattacks=mainData$Hattacks,
                              Fattacks=mainData$Fattacks,
                              lag.states=mainData$lag.states)
  states.boot <- model.main$states #use the existing discretization as our state space
  
  Trans.boot <- gamma2trans(mod0a$coef, summary(mod0a)$sigma,mars.states = states.boot, discretize = F)
 
  mainData.boot$states.discrete <- sapply(mainData.boot$states, 
                                          function(x){return(states.boot[which.min((states.boot-x)^2)])})
  py.params <- data.table(delta=.999,
                          nkappa=2, nbeta=2)

  system(paste0("mkdir -p ipoptTEMP", b))
  write.csv(mainData.boot, file=paste0("ipoptTEMP",b, "/regData.csv"),row.names = F)
  write.csv(Trans.boot, file=paste0("ipoptTEMP",b,"/trans.csv"),row.names = F)
  write.csv(states.boot, file=paste0("ipoptTEMP",b,"/statespace.csv"),row.names = F)
  write.csv(py.params, file=paste0("ipoptTEMP",b,"/params.csv"),row.names = F)
  write.csv(start, file=paste0("ipoptTEMP",b,"/start.csv"),row.names = F)
  
  system(paste("python ../Python3/fitSensitivity.py --k",b), ignore.stdout = T)
  
  regtable <- read.csv(paste0("ipoptTEMP",b,"/regtable.csv"),header = F)
  conv <- read.csv(paste0("ipoptTEMP",b,"/convergence.csv"),header = F)
  v <- read.csv(paste0("ipoptTEMP",b,"/v.csv"),header = F)[,1]
  system(paste0("rm ipoptTEMP",b,"/ -r"))  

  c(regtable$V1, conv$V1)
  
}
stopCluster(cl)
colnames(bootOut) <- c("beta_Hamas", "beta_Fatah", "kappa_Hamas", "kappa_Fatah",
                       "Convergence","logLik")

save(bootOut,file="Results/sensitivityOutput2ndStage.rdata")



rm(list=ls())
load("Results/sensitivityOutput2ndStage.rdata")
converged <- bootOut[bootOut[,"Convergence"]==0,]
converged <- data.frame(converged)

pdf("../../Output/Figures/figureF1.pdf", height=6, width=11)
par(mfrow=c(2,2), cex=.8, cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
truehist(converged$beta_Hamas, main=expression("Histogram of"~beta[H]~"estimates"),
         xlab=expression(beta[H]~"estimates"), ylab="Density")
legend("topright", cex=1.3,adj=c(-.025,.5),bty="n",
       legend = c(paste(formatC(mean(converged$beta_Hamas<0)*100,1,format="f"),
                        c("% estimates < 0"),
                        sep = "")
                  )
       )

truehist(converged$beta_Fatah,  main=expression("Histogram of"~beta[F]~"estimates"),
         xlab=expression(beta[F]~"estimates"), ylab="Density")
legend("topleft", cex=1.3,adj=c(.15,0.5),bty="n",
       legend = c(paste(formatC(mean(converged$beta_Fatah>0)*100,1,format="f"),
                        c("% estimates > 0"),
                        sep = "")
       )
)

truehist(converged$kappa_Hamas,
         main=expression("Histogram of"~kappa[H]~"estimates"),
         xlab=expression(kappa[H]~"estimates"), ylab="Density")
legend("topright", cex=1.3,adj=c(-.025,.5), bty="n",
       legend = c(paste(formatC(mean(converged$kappa_Hamas<0)*100,1,format="f"),
                        c("% estimates < 0"),
                        sep = "")
       )
)

truehist(converged$kappa_Fatah,
         main=expression("Histogram of"~kappa[F]~"estimates"),
         xlab=expression(kappa[F]~"estimates"), ylab="Density")
legend("topleft", cex=1.3, adj=c(.15,0.5),bty="n",
       legend = c(paste(formatC(mean(converged$kappa_Fatah<0)*100,1,format="f"),
                        c("% estimates < 0"),
                        sep = "")
       )
)

dev.off()




