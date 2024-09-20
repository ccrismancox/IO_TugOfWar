#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: First stage regressions
#' ---
#' Clear workspace and load packages:
rm(list=ls())
library(stargazer)
library(sandwich)
library(zoo)
library(margins)
library(urca)
library(data.table)
library(stringr)
library(doParallel)
library(doRNG)
library(parallel)
library(car)
library(moments)
#'
#' Build the data sets
#'
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
summary(mod0)

#' COINTEGRATION: need to reject the null of unit root in the residuals
summary(ur.df(mod0$residuals)) # good
summary(ur.pp(mod0$residuals)) # good

bootOut <- firststageboot(mod0, regData)
se0 <- sqrt(diag(bootOut))[1:6]
names(se0) <- names(mod0$coef)
any(apply(bootOut, 2, \(x){ks.test(scale(x), pnorm)$p.val}) < 0.05)
any(apply(bootOut, 2, \(x){jarque.test(x)$p.val}) < 0.05)
ks.test(scale(mod0$resid), pnorm)
jarque.test(mod0$resid)

#' ECM differences using the cointegration vector of 1
regData[, diff.states := states-lag.states]
regData[,L.diff.states := shift(diff.states)]

#' Test for stationarity in the differenced variable
summary(ur.df(na.omit(regData$diff.states))) #good
summary(ur.pp(na.omit(regData$diff.states))) #good

#' ECM
mod1 <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + L.diff.states + lag.Hattacks:lag.states +
             lag.Fattacks:lag.states, data=regData)

#prepare the table
mod.list <- list(mod0, mod1)
#SEs
se.list <- list(se0*NA, sqrt(diag(NeweyWest(mod1))))



cat( stargazer(mod.list,
               se=se.list,
               no.space = TRUE,
               omit.stat =  "all",
               notes = c("\\scriptsize\\emph{Note:} $^*p<0.05$. Newey-West standard errors in parenthesis. No standard errors are reported in column 1 due to unit root."),
               add.lines=list(
                   c("$T$",
                     paste("\\multicolumn{1}{c}{",
                           sapply(mod.list, function(m){nrow(m$model)}),
                           "}",sep="")),
                   c("adj. $R^2$",
                     paste("\\multicolumn{1}{c}{",
                           sapply(mod.list, function(m){formatC(summary(m)$adj.r.squared,digits=3, format="f")}),
                           "}",sep="")),
                   c("$\\hat{\\sigma}$",
                     paste("\\multicolumn{1}{c}{",
                           sapply(mod.list, function(m){formatC(summary(m)$sigma,digits=3, format="f")}),
                           "}",sep=""))
               ),
               title="Regressing the state space on terrorist attacks",
               label="tab:firststage",
               dep.var.labels = c("State", rep("$\\Delta$ State", 10)),
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
                                    "$\\Delta$ Lag state",
                                    "Hamas attacks $\\times$ lag state",
                                    "Fatah attacks $\\times$ lag state",
                                    "Constant")),
    file="../../Output/Tables/Table1.tex", sep="\n")
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

# how many states are actually visited
sum(states %in% mainData$states.discrete)



save(mainData, file="Results/mainStateActions.Rdata")
save(list=c("Trans", "mod0", "states", "bootOut", "regData"),
     file="Results/firstStageOutput.Rdata")  


#' For the ECM we consider two sets of hypothesis tests:
#'
#' 1. Are the combined estimates significant for each side at each state?
#' 2. Does Fatah have an advantage at each state?

M1 <- margins(mod1,
              variables = c("lag.Hattacks", "lag.Fattacks"),
              data=mod1$model,
              at=data.frame(lag.states=states),
              vcov=NeweyWest(mod1))
mean(summary(M1)$p < 0.05)
H1 <- rep(0, 23)
for(i in 1:23){
  k <- seq(-10, 12)[i]
  H1[i] <-
    linearHypothesis(mod1,
                     paste("lag.Fattacks+", k,
                           "*lag.Fattacks:lag.states+lag.Hattacks+",
                           k,"lag.Hattacks:lag.states"), vcov=NeweyWest)$Pr[2]
}
table(H1 < 0.05) #we reject the null at no differnce at all these values of the states 

