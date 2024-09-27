#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: First stage regressions
#' ---

library(stargazer)
library(sandwich)
library(zoo)
library(urca)
library(data.table)
library(stringr)
library(doParallel)
library(doRNG)
library(parallel)
library(moments)
rm(list=ls())
## load the data sets

load("../../Data/measurement.rdata")
load("../../Data/actionsSetup.Rdata")

source("gamma2trans.R")
source("firststageboot.r")
d = 0.05
bound = 0.025

regData[,Hattacks.count:=dat$lag.Hattacks]
regData[,Fattacks.count:=dat$lag.Fattacks]

mod0 <- lm(states~lag.Hattacks + lag.Fattacks + lag.states+
             lag.Hattacks:lag.states+lag.Fattacks:lag.states, data=regData,
           x=T,y=T)

#' COINTEGRATION: need to reject the null of unit root in the residuals
cat("Do we reject the null of unit root in the residuals (i.e., we have cointegration)?\n")
df.test <- summary(ur.df(mod0$residuals))
cat(ifelse(abs(drop(df.test@teststat)) > abs( df.test@cval[2]), TRUE, FALSE), "\n")
cat("Do we reject the null of normal residuals (KS test)?\n")
ks.test(scale(mod0$resid), pnorm)$p.value < 0.05
cat("Do we reject the null of normal residuals (JB test)?\n")
jarque.test(mod0$resid)$p.value < 0.05

bootOut <- firststageboot(mod0, regData)
se0 <- sqrt(diag(bootOut))[1:6]
names(se0) <- names(mod0$coef)

cat("After the parametric bootstrap, do we find that gamma's are approximately normal?\n")
cat("Do we reject the null of normal gammas in any bootstrap iteration (KS test)?\n")
any(apply(bootOut, 2, \(x){ks.test(scale(x), pnorm)$p.val}) < 0.05)
cat("Do we reject the null of normal gammas in any bootstrap iteration (JB test)?\n")
any(apply(bootOut, 2, \(x){jarque.test(x)$p.val}) < 0.05)



#' ECM differences using the cointegration vector of 1
regData[, diff.states := states-lag.states]
regData[,L.diff.states := shift(diff.states)]


#' Test for stationarity in the independent variables
cat("Testing for stationarity in the independent variables\n")
df.test <- function(x){
    dickey <-summary(ur.df(na.omit(x)))
    return(ifelse(abs(drop(dickey@teststat)) > abs(dickey@cval[2]), TRUE, FALSE))
}

print(regData[, lapply(.SD, df.test), .SDcols=c("lag.Hattacks", "lag.Fattacks", "L.diff.states")])

#' ECM
mod1 <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + L.diff.states + lag.Hattacks:lag.states +
             lag.Fattacks:lag.states, data=regData)

#prepare the table
mod.list <- list(mod0, mod1)
#SEs
se.list <- list(se0*NA, sqrt(diag(NeweyWest(mod1))))



tab1 <- capture.output(stargazer(mod.list,
              se=se.list,
               no.space = TRUE,
               omit.stat =  "all",
               notes = c("Note: Newey-West standard errors in parenthesis. No standard errors are reported in column 1 due to unit root."),
               add.lines=list(
                   c("$T",
                     paste(sapply(mod.list, function(m){nrow(m$model)}),
                           sep="")),
                   c("adj. R2",
                     paste(sapply(mod.list, function(m){formatC(summary(m)$adj.r.squared,digits=3, format="f")}),
                           sep="")),
                   c("sigma$",
                     paste(sapply(mod.list, function(m){formatC(summary(m)$sigma,digits=3, format="f")}),
                           sep=""))
               ),
               title="Regressing the state space on terrorist attacks",
               label="tab:firststage",
               dep.var.labels = c("State", rep("Delta State", 10)),
               digits=2,
               align=TRUE,
               header=FALSE,
               notes.append = FALSE,
               star.cutoffs = c(NA),
               notes.label = "",
               notes.align = "l",
               covariate.labels = c("Hamas attack",
                                    "Fatah attacks",
                                    "Lag state",
                                    "Delta Lag state",
                                    "Hamas attacks x lag state",
                                    "Fatah attacks x lag state",
                                    "Constant"),
               out="../../Output/Tables/table1.txt"))

regData[,Date:= (seq.Date(as.Date("1994-1-1"), as.Date("2018-12-1"), by="month"))]
mainData <- regData[,list(Date,states, lag.states)]
mainData$Hattacks <- ifelse(dat$Hattacks>0,1,0)
mainData$Fattacks <- ifelse(dat$Fattacks>0,1,0)





#### set up for second stage ####
Trans <- gamma2trans(mod0$coef, summary(mod0)$sigma,
                     d=d, bound=bound,
                     mars.states = regData$states)
mainData$states.discrete <- sapply(mainData$states,
                                   function(x){return(states[which.min((states-x)^2)])})


save(mainData, file="Results/mainStateActions.Rdata")
save(list=c("Trans", "mod0", "states", "bootOut", "regData"),
     file="Results/firstStageOutput.Rdata")


#' For the ECM we consider two sets of hypothesis tests:
#'
#' 1. Are the combined estimates significant for each side at each state?
#' 2. Does Fatah have an advantage at each state?




marginalsF <- cbind(0,0, 1, 0, 0, states) %*%  mod1$coef
marginalsH <- cbind(0,1,0,  0, states, 0) %*%  mod1$coef
Vmod1 <-  NeweyWest(mod1)
seMF <- sqrt(diag(cbind(0,0, 1, 0, 0, states) %*% Vmod1 %*% t(cbind(0,0, 1, 0, 0, states))))
seMH <- sqrt(diag(cbind(0, 1,0, 0, states, 0) %*% Vmod1 %*% t(cbind(0, 1,0,  0, states, 0))))

pF <- 2*pnorm(abs(marginalsF/seMF), lower=FALSE)
pH <- 2*pnorm(abs(marginalsH/seMH), lower=FALSE)

cat("At what percentage of the states can we reject the null that Fatah and Hamas are capable of moving the state space \n")
colMeans(cbind(Fatah=pF, Hamas=pH) < 0.05) * 100


A <- cbind(0,1, 1, 0, states, states)
cancel <- (marginalsF + marginalsH)
cancelSE <- sqrt(diag(A %*% Vmod1 %*% t(A)))
p.cancel <- 2*pnorm(abs(cancel/cancelSE), lower=FALSE)
cat("At what percentage of the states can we reject the null that Fatah and Hamas are equallty capable of moving the state space \n")
mean(p.cancel< 0.05) * 100


