#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: Second stage main output and time limits
#' ---
rm(list=ls())
library(data.table)
library(knitr)
load("Results/mainStateActions.Rdata")
load("Results/firstStageOutput.Rdata")


load("Results/startingvalues.rdata")
start <- c(startmod$regtable$V1, startmod$v)

py.params <- data.table(delta=.999, 
                        nkappa=2, nbeta=2)
regData <- mainData[,.(states, Hattacks, Fattacks, lag.states, states.discrete)]

system("mkdir -p ipoptTEMP")
write.csv(regData, file="ipoptTEMP/regData.csv",row.names = F)
write.csv(Trans, file="ipoptTEMP/trans.csv",row.names = F)
write.csv(states, file="ipoptTEMP/statespace.csv",row.names = F)
write.csv(py.params, file="ipoptTEMP/params.csv",row.names = F)
write.csv(start, file="ipoptTEMP/start.csv",row.names = F)
write.csv(c(mod0$coef, summary(mod0)$sigma), file="ipoptTEMP/gamma.csv",row.names = F)
write.csv(bootOut, file="ipoptTEMP/V1.csv",row.names = F)

system("python ../Python3/fitMainModel.py", ignore.stdout = FALSE)

regtable <- read.csv("ipoptTEMP/regtable.csv",header = F)
regtable2 <- read.csv("ipoptTEMP/regtable2.csv",header = F) #one sided tests
conv <- read.csv("ipoptTEMP/convergence.csv",header = F)
v <- read.csv("ipoptTEMP/v.csv",header = F)[,1]
V1 <- read.csv("ipoptTEMP/VCOV1.csv",header = F)
V2 <- read.csv("ipoptTEMP/VCOV2.csv",header = F)


model.main <- list(regData=regData, Trans=Trans, states=states, params=py.params,
                   V1=V2[1:4, 1:4], V2=V2[1:4, 1:4], 
                   conv=conv, v=v, regtable=regtable, regtable2=regtable2)


startmod <- model.main

save(list=c("model.main"), file="Results/mainModel.rdata")
save(list=c("startmod"), file="Results/startingvalues.rdata")
system("rm ipoptTEMP/ -r")

rm(list=ls())
load("Results/mainModel.rdata")
attach(model.main)
rownames(regtable2) <- rownames(regtable) <-c("beta[H]", "beta[F]", "kappa[H]", "kappa[F]")
output <- cbind(regtable[,1:2], regtable2[,2])
output[1:2, ] <- round(output[1:2,], 4)
output[3:4, ] <- round(output[3:4,], 2)

cat(kable(output, format = "pipe",
            col.names=c("Estimate", "Std. Error 1", "Std. Error 2"),
          caption="Payoff estimates"),
    file="../../Output/Tables/Table2.txt", sep="\n")

cat("T = ", nrow(regData), "\n",  file="../../Output/Tables/Table2.txt", sep="\t", append=TRUE)
cat("LogLik = ", round(conv$V1[2],2), "\n",  file="../../Output/Tables/Table2.txt", sep="\t", append=TRUE)

# Joint hypothesis that betas =0
A <- rbind(c(1,0,0,0),
           c(0,1,0,0))
b <- c(0,0)
theta <- c(regtable$V1)
V2 <- as.matrix(V2)[1:4, 1:4]
wald.stat <- t(A %*% theta -b) %*% solve(A %*% V2 %*% t(A)) %*% (A %*% theta -b)
pchisq(wald.stat, df=2, lower=FALSE)
