#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Fit no competition models
#' ---
#' 
rm(list=ls())
library(data.table)
library(knitr)
load("Results/mainStateActions.Rdata")
load("Results/firstStageOutput.Rdata")
load("Results/mainModel.rdata")
load("Results/noCompetition.Rdata")
source("gamma2trans.R")
source("helperFunctions.r")



gamma0 <- mod0$coef
sigma0 <- summary(mod0)$sigma

##### No one is competitive #####
modRW <- lm(states~lag.states, data=regData)
gamma <- c(modRW$coef[1],0,0, modRW$coef[2],0,0)
sigma <- summary(modRW)$sigma

Trans <- gamma2trans(gamma, sigma, 
                     d = 0.05,
                     bound = 0.025,
                     mars.states = mainData$states)
mainData$states.discrete <- sapply(mainData$states, 
                                   function(x){return(states[which.min((states-x)^2)])})

py.params <- data.table(delta=.999,
                        nkappa=2, nbeta=0,
                        agent=0) #agent code; 0 means no competition
pyData <- mainData[,list(states, Hattacks, Fattacks, lag.states, states.discrete)]
start0 <- c(noCompModel$regtable$V1, noCompModel$v)


system("mkdir -p ipoptTEMP")
write.csv(pyData, file="ipoptTEMP/regData.csv",row.names = F)
write.csv(Trans, file="ipoptTEMP/trans.csv",row.names = F)
write.csv(states, file="ipoptTEMP/statespace.csv",row.names = F)
write.csv(py.params, file="ipoptTEMP/params.csv",row.names = F)
write.csv(start0, file="ipoptTEMP/start.csv",row.names = F)
system("python ../Python3/fitNoCompetition.py", ignore.stdout = FALSE)

regtable <- read.csv("ipoptTEMP/regtable.csv",header = F)
conv <- read.csv("ipoptTEMP/convergence.csv",header = F)
v <- read.csv("ipoptTEMP/v.csv",header = F)[,1]
vcov <- read.csv("ipoptTEMP/VCOV1.csv",header = F)

noCompModel <- list(regData=pyData, Trans=Trans, states=states, params=py.params,
                    vcov=vcov,
                       conv=conv, v=v, regtable=regtable)


# Clean up and save
system("rm ipoptTEMP/ -r")
save(list=c("noCompModel"), 
     file="Results/noCompetition.Rdata")

