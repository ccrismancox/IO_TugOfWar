#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Fit the tit for tat model
#' ---
rm(list=ls())
library(data.table)
library(moments)
library(knitr)
load("Results/mainStateActions.Rdata")
load("../../Data/measurement.rdata")


states <- 1:4
#### 
# state 1: no attacks last month
# state 2: F attack last month
# state 3: H attach last month
# state 4: both attack last month
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()

regData[, states:=2*lag.Hattacks+lag.Fattacks+1]

Trans <- rbind(diag(4), diag(4),diag(4),diag(4))

start <- runif(4+length(states)*4)/10

py.params <- data.table(delta=.999,
                        nkappa=2, nbeta=2)

#Same as the old ordering
regData <- data.table(states=regData$states,
                      Hattacks=mainData$Hattacks,
                      Fattacks=mainData$Fattacks,
                      lag.states = regData[,shift(states)],
                      states.discrete=regData$states)



system("mkdir -p ipoptTEMP")
write.csv(regData, file="ipoptTEMP/regData.csv",row.names = F)
write.csv(Trans, file="ipoptTEMP/trans.csv",row.names = F)
write.csv(states, file="ipoptTEMP/statespace.csv",row.names = F)
write.csv(py.params, file="ipoptTEMP/params.csv",row.names = F)
write.csv(start, file="ipoptTEMP/start.csv",row.names = F)

system("python ../Python3/fitMainModel_t4t.py", ignore.stdout = FALSE)

regtable <- read.csv("ipoptTEMP/regtable.csv",header = F)
conv <- read.csv("ipoptTEMP/convergence.csv",header = F)
v <- read.csv("ipoptTEMP/v.csv",header = F)[,1]
V1 <- read.csv("ipoptTEMP/VCOV1.csv",header = F)

model.main.t4t <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                       V1=V1[1:4, 1:4], 
                       conv=conv, v=v, regtable=regtable)



save(list=c("model.main.t4t"), file="Results/t4t_Model.rdata")
system("rm ipoptTEMP/ -r")

