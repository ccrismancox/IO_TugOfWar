#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: First stage robustness checks
#' ---
library(stargazer)
library(sandwich)
library(zoo)
library(lmtest)
library(data.table)
library(ivreg)
library(stringr)
library(cragg)
library(car)
library(sensemakr)
library(rsq)
library(gridExtra)
library(metR)
library(knitr)
library(tipr)
library(ggplot2)
library(modelsummary)

rm(list=ls())

#### load the data ####
load("../../Data/measurement.rdata")
load("../../Data/actionsSetup.Rdata")
load("../../Data/PalestinianDeaths.rdata")
load("../../Data/ExtraFactors.rdata") 
load("../../Data/actionsSetup_byAttackType.Rdata")

## For additional controls
corrupt <- fread("../../Data/corruption_WBG.csv",header = TRUE)
mortal <- fread("../../Data/mortality_WB.csv",header = TRUE)
corrupt[,(colnames(corrupt)[-1]) := lapply(.SD, as.numeric), .SDcols=colnames(corrupt)[-1]]
mortal[,(colnames(mortal)[-1]) := lapply(.SD, as.numeric), .SDcols=colnames(mortal)[-1]]

## For the IV ## 
load("../../Data/rainData.rdata") 
load("../../Data/otherattacks.rdata") ## need the files for this
source("liml.r")


source("gamma2trans.R")
source("firststageboot.r")
source("helperFunctions.r")
d = 0.05
bound = 0.025

regData[,Hattacks.count:=dat$lag.Hattacks]
regData[,Fattacks.count:=dat$lag.Fattacks]

#unlagged attacks
regData[,Fattack0:=dat$Fattacks]
regData[,Hattack0:=dat$Hattacks]

regData <- data.table(merge( gg.extra[1:300,], regData, by="date"))
regData[,lag.emp := shift(employStates)]
regData[,lag.violence := shift(violence1DStates)]

regData <- merge(regData, palDeaths,by="date", all.x=TRUE)
regData[,lag.killed := log(shift(deaths))]


regData[,elections := ifelse(date %in% as.yearmon(c("May 1996",
                                                    "May 1999",
                                                    "Jan 2003",
                                                    "Mar 2006",
                                                    "Feb 2009",
                                                    "Jan 2013",
                                                    "Mar 2015"), "%b %Y"),
                             1, 0)]
regData[,timesinceelection :=seq_len(.N)-1, by=.(rleid(cumsum(elections)))]
regData[date < "May1996",timesinceelection :=  timesinceelection + 19]

regData[ , second.int := ifelse(ts > 81 & ts<=144, 1,0)]
regData[ , olso.era := ifelse(ts <= 81, 1,0)]
regData[ , post := ifelse(ts > 144, 1,0)]



#' ECM differences using the cointegration vector of about 1
#' Differencing the other I(1) variables
regData[, diff.states := states-lag.states]
regData[,L.diff.states := shift(diff.states)]
regData[,L.diff.emp := lag.emp-shift(lag.emp)]
regData[,L.diff.violence := lag.violence-shift(lag.violence)]


#' Robustness to time changes structure 
mod2t <- lm(diff.states ~ lag.Hattacks*ts + lag.Fattacks*ts +
              L.diff.states+
              lag.Hattacks:lag.states +lag.Fattacks:lag.states,
            data=regData)
cat("Can we reject the null that the time trends here are 0?\n")
linearHypothesis(mod2t,c("lag.Hattacks:ts","ts:lag.Fattacks"), vcov=NeweyWest)


#' Robustness to covariates
mod2a <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + L.diff.states, data=regData)
mod2b <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + L.diff.emp+L.diff.states + lag.Hattacks:lag.states +
              lag.Fattacks:lag.states, data=regData)
mod2c <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + L.diff.emp+L.diff.violence+L.diff.states + lag.Hattacks:lag.states +
              lag.Fattacks:lag.states, data=regData)

mod2d <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + second.int +L.diff.emp+L.diff.violence+L.diff.states+ lag.Hattacks:lag.states +
              lag.Fattacks:lag.states, data=regData)
mod2e <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + second.int +L.diff.emp+L.diff.violence+L.diff.states + timesinceelection+  lag.Hattacks:lag.states +
              lag.Fattacks:lag.states, data=regData, x=T, y=T)
mod2f <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + second.int
            +L.diff.emp+L.diff.violence+L.diff.states + timesinceelection
            + lag.killed + lag.Hattacks:lag.states +
              lag.Fattacks:lag.states, data=regData)





corrupt <- melt(corrupt, "Country Name", 
                variable.name="year",
                value.name = "corruption")
corrupt[, `Country Name` := NULL]
corrupt[,year:=as.integer(as.character(year))]
mortal <- melt(mortal, "Country Name", 
               variable.name="year",
               value.name = "infant")
mortal[, `Country Name` := NULL]
mortal[,year:=as.integer(as.character(year))]
regData<-merge(regData, corrupt, by="year", all.x=TRUE, all.y=F)
regData<-merge(regData, mortal, by="year", all.x=TRUE, all.y=F)
regData[,corruption := -1*corruption]
regData[,post.election := 1*(ts >= 146)]

mod2g <- lm(diff.states ~ lag.Hattacks + lag.Fattacks + second.int
            +L.diff.emp+L.diff.violence+L.diff.states + timesinceelection
            + lag.killed 
            + corruption+corruption:post.election 
            + infant+infant:post.election
            +lag.Hattacks:lag.states +
              lag.Fattacks:lag.states, data=regData)
cat("What are the effects post 2006?\n")
deltaMethod(mod2g, "corruption+`corruption:post.election`", vcov=NeweyWest)
deltaMethod(mod2g, "infant+`post.election:infant`", vcov=NeweyWest)










mod.list.controls <- list(mod2a, mod2t,mod2b, mod2c, mod2d, mod2e, mod2f, mod2g)
se.list.controls <- lapply(mod.list.controls, function(x){sqrt(diag(NeweyWest(x)))})
tabD3 <- capture.output(stargazer(mod.list.controls,
              se=se.list.controls,
              no.space = TRUE,
              omit.stat =  "all",
              notes = c("Note:  Newey-West standard errors in parenthesis."),
              add.lines=list( c("State-attack interactions",
                                paste(ifelse(sapply(mod.list.controls,
                                                    function(x){any(str_detect(names(x$coef), ":"))}),
                                             "Yes", "No"),
                                      sep="")),
                             c("T",
                               paste( sapply(mod.list.controls, function(m){nrow(m$model)}),
                                     sep="")),
                             c("adj. R2",
                               paste(sapply(mod.list.controls, function(m){formatC(summary(m)$adj.r.squared,digits=3, format="f")}),
                                     sep="")),
                             c("sigma",
                               paste(sapply(mod.list.controls, function(m){formatC(summary(m)$sigma,digits=3, format="f")}),
                                     sep=""))
                             ),
              title="Robustness checks for the first-stage model: Specification changes",
              label="tab:firststage.app",
              dep.var.labels = c("Delta State"),
              digits=2,
              align=TRUE,
              omit=c(":lag.states"),
              header=FALSE,
              notes.append = FALSE,
              star.cutoffs = c(NA),
              notes.label = "",
              notes.align = "l",
              order=c("lag.Hattacks", "lag.Fattacks", "L.diff.states",
                      "L.diff.emp",
                      "L.diff.violence"),
              covariate.labels = c("Hamas attacks",
                                   "Hamas attacks x time",
                                   "Fatah attacks",
                                   "Fatah attacks x time",
                                   "Delta Lag state",
                                   "Delta unemployment",
                                   "Delta support for violence",
                                   "Time",
                                   "Second Intifada",
                                   "Time since last election",
                                   "Palestinian fatalities by Israel",
                                   "Corruption",
                                   "Infant mortality",
                                   "Corruption x post-2006 election",
                                   "Infant mortality x post-2006 election",
                                   "Constant"),
         out="../../Output/Tables/tableD3.txt"))




#' Robustness to action measures
mod3a <- lm(diff.states ~ Hattacks.count + Fattacks.count
            +  second.int+L.diff.emp+L.diff.violence+timesinceelection+
              L.diff.states+
              Hattacks.count:lag.states + Fattacks.count:lag.states, data=regData)
regData[,Hkill.attack := ifelse(Hattacks.count >0 ,Hkills/Hattacks.count,0)]
regData[,Fkill.attack := ifelse(Fattacks.count >0 ,Fkills/Fattacks.count,0)]
mod3b <- lm(diff.states ~ lag.Hattacks  + lag.Fattacks  +  second.int+L.diff.emp+L.diff.violence+timesinceelection+
              L.diff.states+
              Fkills+Hkills+
              lag.Hattacks:lag.states + lag.Fattacks:lag.states, data=regData)
mod3c <- lm(diff.states ~ Hkills  + Fkills  +  second.int+L.diff.emp+L.diff.violence+timesinceelection+
              L.diff.states+
              Hkills:lag.states + Fkills:lag.states, data=regData)
mod3d <- lm(diff.states ~ Hkill.attack  + Fkill.attack  +  second.int+L.diff.emp+L.diff.violence+timesinceelection+
              L.diff.states+
              Hkill.attack:lag.states + Fkill.attack:lag.states, data=regData)


### attack targets/types
regData$lag.Hattacks.mil <- (dat.Civs$lag.Hattacks.mil >0)*1
regData$lag.Fattacks.mil <- (dat.Civs$lag.Fattacks.mil >0)*1
regData$lag.Hattacks.notciv <- (dat.Civs$lag.Hattacks.notciv >0)*1
regData$lag.Fattacks.notciv <- (dat.Civs$lag.Fattacks.notciv >0)*1
regData$lag.Hattacks.civ <- (dat.Civs$lag.Hattacks.civ >0)*1
regData$lag.Fattacks.civ  <- (dat.Civs$lag.Fattacks.civ >0)*1

mod0.notciv <- lm(diff.states~lag.Hattacks.notciv + lag.Fattacks.notciv + second.int + 
                      L.diff.emp + L.diff.violence + timesinceelection + L.diff.states + 
                      lag.Hattacks.notciv:lag.states + lag.Fattacks.notciv:lag.states, data=regData,
                  x=T,y=T)

mod0.civ <- lm(diff.states~lag.Hattacks.civ + lag.Fattacks.civ + second.int + 
                   L.diff.emp + L.diff.violence + timesinceelection + L.diff.states + 
                   lag.Hattacks.civ:lag.states + lag.Fattacks.civ:lag.states, data=regData,
               x=T,y=T)









mod.list.measurement <- list(mod3a, mod3b, mod3c, mod3d, mod0.civ, mod0.notciv)
se.list.measurement <- lapply(mod.list.measurement, function(x){sqrt(diag(NeweyWest(x)))})
names(mod.list.measurement) <- names(se.list.measurement) <- c("Counts",
                                                               "Binary & fatalities",
                                                               "Fatalities",
                                                               "Fatalities/attack",
                                                               "Civilian targets",
                                                               "Non-civilian targets")
inters <- rbind.data.frame(c("State-attack interactions", rep("Yes", length(mod.list.measurement))))
attr(inters,"position") <- c(21)
modelsummary(mod.list.measurement, 
             vcov="NeweyWest",
             stars=FALSE,
             output="../../Output/Tables/tableD4.txt",
             fmt=num2str,
             escape=FALSE,
             note="Newey-West standard errors in parenthesis",
             title=c("Robustness checks for the first-stage model: Measurement changes"),
             coef_map= c("Hattacks.count" ="Hamas attacks",
                         "Fattacks.count" ="Fatah attacks",
                         "lag.Hattacks"="Hamas attacks",
                         "lag.Fattacks"="Fatah attacks",
                         "Hkills"="Hamas fatalities",
                         "Fkills"="Fatah fatalities",
                         "Hkill.attack"="Hamas attacks",
                         "Fkill.attack"="Fatah attacks",
                         "lag.Hattacks.civ"="Hamas attacks",
                         "lag.Fattacks.civ"="Fatah attacks",
                         "lag.Fattacks.notciv"="Fatah attacks",
                         "lag.Hattacks.notciv"="Hamas attacks",
                         "second.int"= "Second Intifada",
                         "L.diff.emp"="Delta unemployment",
                         "L.diff.violence"="Delta support for violence",
                         "timesinceelection"="Time since late election",
                         "L.diff.states"="Delta Lag state",
                         "(Intercept)"="Constant"),
             gof_map=list(list("raw"="nobs", "clean"="T", "fmt"=\(x){round(x)}),
             list("raw"="adj.r.squared", "clean"="adj. R2", "fmt"=\(x){round(x,2)}),
list("raw"="rmse", "clean"="sigma", "fmt"=\(x){round(x,2)})),
add_rows=inters)





for(i in 1:length(mod.list.measurement)){
    mod <- mod.list.measurement[[i]]
    V <- NeweyWest(mod)
    strF <- ifelse(i==3 | i==4, "Fkill", "Fatt")
    strH <- ifelse(i==3 | i==4, "Hkill", "Hatt")
    Fidx <- grep(strF,names(mod$coefficients))
    Hidx <- grep(strH,names(mod$coefficients))
    X <-  matrix(0, ncol=length(mod$coef), nrow=440)
    colnames(X) <- names(mod$coefficients)  
    X[,Fidx[1]] <- 1; X[,Fidx[2]] <- states
    X[,Hidx[1]] <- 1; X[,Hidx[2]] <- states
    cancel <- X %*%  mod$coef
    cancelSE <- sqrt(diag(X %*% V %*% t(X)))
    p.cancel <- 2*pnorm(abs(cancel/cancelSE), lower=FALSE)
    reject <- mean(p.cancel< 0.1) * 100
    cat("At what percentage of the states can we reject the null that Fatah and Hamas are equallty capable of moving the state space in model", i, "of Table D.4\n")
    print(reject)
    if(reject < 100){
      cat("In model ", i, "of Table D.4, we fail to reject the null of equally effective at states less than\n")
      print(max(states[which(p.cancel>0.1) ]))
    }
}
  



###### IV ######
datMonth[,month:=NULL]
regData <- merge(regData, datMonth,all=TRUE, by="date")
regData <- merge(regData, other.attack,all=TRUE, by="date")
regData[,L.dev.gaza := shift(dev.gaza,1)]
regData[,L.dev.wb := shift(dev.wb,1)]
regData[,L2.dev.gaza := shift(dev.gaza,2)]
regData[,L2.dev.wb := shift(dev.wb,2)]
regData[,L3.dev.gaza := shift(dev.gaza,3)]
regData[,L3.dev.wb := shift(dev.wb,3)]
regData[,L4.dev.gaza := shift(dev.gaza,4)]
regData[,L4.dev.wb := shift(dev.wb,4)]
regData[,L5.dev.gaza := shift(dev.gaza,5)]
regData[,L5.dev.wb := shift(dev.wb,5)]

regData[,Pattacks := as.numeric(Pattacks > 0)]
regData[,lag.Pattacks := shift(Pattacks,1)]
regData[,lag2.Pattacks := shift(Pattacks,2)]
regData[,lag3.Pattacks := shift(Pattacks,3)]
regData[,lag4.Pattacks := shift(Pattacks,4)]
regData[,lag5.Pattacks := shift(Pattacks,5)]

regData[,lag2.states := shift(lag.states)]

## exact ID with no controls (use ivreg to get the right obs, first stage F, and Sargan)
mod5 <- ivreg(diff.states ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                lag.Fattacks:lag.states+ L.diff.states |
                L.diff.states
              +L2.dev.gaza + L2.dev.wb
              +L2.dev.gaza:lag2.states
              +L2.dev.wb:lag2.states,
              data=regData)
D5a <- summary(mod5, vcovHC)$diagnostics
X1 = model.matrix(~ L.diff.states,
                  data=regData[as.numeric(rownames(mod5$model))])
X2 = model.matrix( ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                     lag.Fattacks:lag.states+0, data=regData[as.numeric(rownames(mod5$model))])
Z = model.matrix(~   L.diff.states +L2.dev.gaza + L2.dev.wb
                 +L2.dev.wb:lag2.states
                 +L2.dev.gaza:lag2.states,
                 data=regData[as.numeric(rownames(mod5$model))])
y = mod5$model$diff.states
mod5a <- liml.fit( X1, X2, Z, y, C=1, vcov=c("bekker"))

mod5A <- mod5
mod5A$coefficients <- mod5a$coef.table[,1][names(mod5A$coefficients)]

## Cragg and Donald test
CD5a <- cragg_donald(X=~L.diff.states, # Control Variables
                     D=~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                       lag.Fattacks:lag.states, # Treatments
                     Z=~L2.dev.gaza + L2.dev.wb
                     +L2.dev.wb:lag2.states
                     +L2.dev.gaza:lag2.states,
                     data = regData[as.numeric(rownames(mod5$model))])
sigma5a <- sqrt(crossprod(y-cbind(X1,X2)[,rownames(mod5a$coef.table)] %*% mod5a$coef.table[,1])/(nrow(X1)-ncol(X1)-ncol(X2)))

## over ID with no controls
mod5 <- ivreg(diff.states ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                lag.Fattacks:lag.states+ L.diff.states |
                L.diff.states+ L2.dev.gaza + L2.dev.wb
              +L2.dev.wb:lag2.states
              +L2.dev.gaza:lag2.states
              +lag2.Pattacks,
              data=regData)
D5b <- summary(mod5, vcovHC)$diagnostics
shpb <- summary(mod5,vcovHC)$diagnostics['Sargan', "p-value"]
X1 = model.matrix(~ L.diff.states,
                  data=regData[as.numeric(rownames(mod5$model))])
X2 = model.matrix( ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                     lag.Fattacks:lag.states+0, data=regData[as.numeric(rownames(mod5$model))])
Z = model.matrix(~   L.diff.states+ L2.dev.gaza + L2.dev.wb
                 +L2.dev.wb:lag2.states
                 +L2.dev.gaza:lag2.states
                 +lag2.Pattacks,
                 data=regData[as.numeric(rownames(mod5$model))])
y = mod5$model$diff.states

mod5b <- liml.fit( X1, X2, Z, y, C=1, vcov=c("bekker"))
mod5B <- mod5
mod5B$coefficients <- mod5b$coef.table[,1][names(mod5B$coefficients)]
CD5b <- cragg_donald(X=~L.diff.states, # Control Variables
                     D=~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                       lag.Fattacks:lag.states, # Treatments
                     Z=~ L2.dev.gaza + L2.dev.wb
                     +L2.dev.wb:lag2.states
                     +L2.dev.gaza:lag2.states
                     +lag2.Pattacks,
                     data = regData[as.numeric(rownames(mod5$model))])
sigma5b <- sqrt(crossprod(y-cbind(X1,X2)[,rownames(mod5b$coef.table)] %*% mod5b$coef.table[,1])/(nrow(X1)-ncol(X1)-ncol(X2)))

## just id'ed with controls
mod5 <- ivreg(diff.states ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                lag.Fattacks:lag.states
              +  second.int+L.diff.emp+L.diff.violence+timesinceelection + L.diff.states|
                L.diff.states+ L2.dev.gaza + L2.dev.wb
              +L2.dev.wb:lag2.states
              +L2.dev.gaza:lag2.states
              +  second.int+L.diff.emp+L.diff.violence+timesinceelection + L.diff.states,
              data=regData)
D5c <- summary(mod5, vcovHC)$diagnostics
X1 = model.matrix(~ L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection,
                  data=regData[as.numeric(rownames(mod5$model))])
X2 = model.matrix( ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                     lag.Fattacks:lag.states+0, data=regData[as.numeric(rownames(mod5$model))])
Z = model.matrix(~ L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection
                 + L2.dev.gaza + L2.dev.wb
                 +L2.dev.wb:lag2.states
                 +L2.dev.gaza:lag2.states                 ,
                 data=regData[as.numeric(rownames(mod5$model))])
y = mod5$model$diff.states
mod5c <- liml.fit( X1, X2, Z, y, C=1, vcov=c("bekker"))
mod5C <- mod5
mod5C$coefficients <- mod5c$coef.table[,1][names(mod5C$coefficients)]
sigma5c <- sqrt(crossprod(y-cbind(X1,X2)[,rownames(mod5c$coef.table)] %*% mod5c$coef.table[,1])/(nrow(X1)-ncol(X1)-ncol(X2)))
CD5c <- cragg_donald(X=~L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection, # Control Variables
                     D=~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                       lag.Fattacks:lag.states, # Treatments
                     Z=~ L2.dev.gaza + L2.dev.wb
                     +L2.dev.wb:lag2.states
                     +L2.dev.gaza:lag2.states,# Instruments
                     data = regData[as.numeric(rownames(mod5$model))])


## over id'ed with controls
mod5 <- ivreg(diff.states ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                lag.Fattacks:lag.states
              +  second.int+L.diff.emp+L.diff.violence+timesinceelection + L.diff.states|
                L.diff.states
              + L2.dev.gaza + L2.dev.wb
              +L2.dev.wb:lag2.states
              +L2.dev.gaza:lag2.states
              +lag2.Pattacks
              +  second.int+L.diff.emp+L.diff.violence+timesinceelection + L.diff.states,
              data=regData)
D5d <- summary(mod5, vcovHC)$diagnostics
shpd <- summary(mod5, vcovHC)$diagnostics['Sargan', "p-value"]
X1 = model.matrix(~ L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection,
                  data=regData[as.numeric(rownames(mod5$model))])
X2 = model.matrix( ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                     lag.Fattacks:lag.states+0, data=regData[as.numeric(rownames(mod5$model))])
Z = model.matrix(~ L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection
                 +L2.dev.gaza + L2.dev.wb
                 +L2.dev.wb:lag2.states
                 +L2.dev.gaza:lag2.states
                 +lag2.Pattacks,
                 data=regData[as.numeric(rownames(mod5$model))])
y = mod5$model$diff.states
mod5d <- liml.fit( X1, X2, Z, y, C=1, vcov=c("bekker"))
mod5D <- mod5
mod5D$coefficients <- mod5d$coef.table[,1][names(mod5D$coefficients)]
sigma5d <- sqrt(crossprod(y-cbind(X1,X2)[,rownames(mod5d$coef.table)] %*% mod5d$coef.table[,1])/(nrow(X1)-ncol(X1)-ncol(X2)))
CD5d <- cragg_donald(X=~L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection, # Control Variables
                     D=~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                       lag.Fattacks:lag.states, # Treatments
                     Z=~ L2.dev.gaza + L2.dev.wb
                     +L2.dev.wb:lag2.states
                     +L2.dev.gaza:lag2.states
                     +lag2.Pattacks,
                     data = regData[as.numeric(rownames(mod5$model))])

## Extra over id'ed with controls
mod5 <- ivreg(diff.states ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                lag.Fattacks:lag.states
              +  second.int+L.diff.emp+L.diff.violence+timesinceelection + L.diff.states|
                L2.dev.gaza + L2.dev.wb
              +L2.dev.wb:lag2.states
              +L2.dev.gaza:lag2.states
              + L3.dev.gaza + L3.dev.wb
              +L3.dev.wb:lag2.states
              +L3.dev.gaza:lag2.states
              + L4.dev.gaza + L4.dev.wb
              +L4.dev.wb:lag2.states
              +L4.dev.gaza:lag2.states
              + L5.dev.gaza + L5.dev.wb
              +L5.dev.wb:lag2.states
              +L5.dev.gaza:lag2.states
              +lag2.Pattacks +lag3.Pattacks + lag4.Pattacks + lag5.Pattacks
              +  second.int+L.diff.emp+L.diff.violence+timesinceelection + L.diff.states,
              data=regData)
D5e <- summary(mod5, vcovHC)$diagnostics
shpe <- summary(mod5, vcovHC)$diagnostics['Sargan', "p-value"]
X1 = model.matrix(~ L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection,
                  data=regData[as.numeric(rownames(mod5$model))])
X2 = model.matrix( ~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                     lag.Fattacks:lag.states+0, data=regData[as.numeric(rownames(mod5$model))])
Z = model.matrix(~ L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection
                 + L2.dev.gaza + L2.dev.wb
                 +L2.dev.wb:lag2.states
                 +L2.dev.gaza:lag2.states
                 + L3.dev.gaza + L3.dev.wb
                 +L3.dev.wb:lag2.states
                 +L3.dev.gaza:lag2.states
                 + L4.dev.gaza + L4.dev.wb
                 +L4.dev.wb:lag2.states
                 +L4.dev.gaza:lag2.states
                 + L5.dev.gaza + L5.dev.wb
                 +L5.dev.wb:lag2.states
                 +L5.dev.gaza:lag2.states
                 +lag2.Pattacks +lag3.Pattacks + lag4.Pattacks + lag5.Pattacks,
                 data=regData[as.numeric(rownames(mod5$model))])
y = mod5$model$diff.states
mod5e <- liml.fit( X1, X2, Z, y, C=1, vcov=c("bekker"))
mod5E <- mod5
mod5E$coefficients <- mod5e$coef.table[,1][names(mod5E$coefficients)]
sigma5e <- sqrt(crossprod(y-cbind(X1,X2)[,rownames(mod5e$coef.table)] %*% mod5e$coef.table[,1])/(nrow(X1)-ncol(X1)-ncol(X2)))
CD5e <- cragg_donald(X=~L.diff.states+  second.int+L.diff.emp+L.diff.violence+timesinceelection, # Control Variables
                     D=~ lag.Hattacks + lag.Fattacks + lag.Hattacks:lag.states +
                       lag.Fattacks:lag.states, # Treatments
                     Z=~L2.dev.gaza + L2.dev.wb
                     +L2.dev.wb:lag2.states
                     +L2.dev.gaza:lag2.states
                     + L3.dev.gaza + L3.dev.wb
                     +L3.dev.wb:lag2.states
                     +L3.dev.gaza:lag2.states
                     + L4.dev.gaza + L4.dev.wb
                     +L4.dev.wb:lag2.states
                     +L4.dev.gaza:lag2.states
                     + L5.dev.gaza + L5.dev.wb
                     +L5.dev.wb:lag2.states
                     +L5.dev.gaza:lag2.states
                     +lag2.Pattacks +lag3.Pattacks + lag4.Pattacks + lag5.Pattacks,
                     data = regData[as.numeric(rownames(mod5$model))])


### table
mod.list.structure <- list(mod5A, mod5B, mod5C, mod5D, mod5E)
se.list.structure <- list(mod5a$coef.table[,2], mod5b$coef.table[,2],
                          mod5c$coef.table[,2], mod5d$coef.table[,2],
                          mod5e$coef.table[,2])
tabD5 <- capture.output(stargazer(mod.list.structure,
              se=se.list.structure,
              no.space = TRUE,
              omit.stat =  "all",
              notes = c("Note: Bekker's robust standard errors (IV models) "),
              add.lines=list( c("Interactions",
                                paste(rep(c("Yes"),5),sep="")),
                             c("Controls",
                               paste(c("No", "No","Yes", "Yes", "Yes"),sep="")),
                             c("$T$",
                               paste(sapply(mod.list.structure, function(m){nrow(m$model)}),sep="")),
                             c("sigma$",
                               paste(formatC(c(sigma5a, sigma5b, sigma5c, sigma5d, sigma5e),digits=3, format="f"),
                                     sep="")),
                             c("Cragg and Donald statistic",
                               paste(c(formatC(c(CD5a$cd_stat, CD5b$cd_stat,
                                                 CD5c$cd_stat, CD5d$cd_stat, CD5e$cd_stat),digits=3, format="f"))
                                    ,sep="")),
                             c("Number of instruments",
                               paste(sapply(mod.list.structure, function(m){length(m$instruments)}),
                                     sep="")),
                             c("Sargan-Hansen p value",
                               paste(c("",
                                       formatC(shpb,digits=3, format="f"),
                                       "",
                                       formatC(shpd,digits=3, format="f"),
                                       formatC(shpe,digits=3, format="f")),
                                     sep=""))
                             ),
              title="Robustness checks for the first-stage model: Time effects and endogeneity",
              label="tab:firststage.app3",
              dep.var.labels = c("Delta State"),
              digits=2,
              order=c("lag.Hattacks", "lag.Fattacks"),
              align=TRUE,
              omit=c(":lag.states", "L.", "second.int", "timesince"),
              header=FALSE,
              notes.append = FALSE,
              star.cutoffs = c(NA),
              notes.label = "",
              notes.align = "l",
              covariate.labels = c("Hamas attacks",
                                   "Fatah attacks",
                                   "Constant"),
              out="../../Output/Tables/tableD5.txt"))

print(min(round(rbind(D5a[1:4,], D5c[1:4,]), 1)[,3])) # first stage F; rainfall only
print(max(round(rbind(D5a[1:4,], D5c[1:4,]), 1)[,3])) # first stage F; over ID 
print(min(round(rbind(D5b[1:4,], D5d[1:4,], D5e[1:4,]), 1)[,3])) # first stage F; rainfall only
print(max(round(rbind(D5b[1:4,], D5d[1:4,], D5e[1:4,]), 1)[,3])) # first stage F; over ID 











##### sensitivity D3 ##### 


mod.sen <- mod2e

#Percentage change needed for symmetry
q <- 1+mod.sen$coefficients["lag.Hattacks"]/mod.sen$coefficients["lag.Fattacks"]
q2 <- 1+abs(mod.sen$coefficients["lag.Fattacks"]/mod.sen$coefficients["lag.Hattacks"])


## For reference; what are best predictors of lagged attacks?
lm.f <- update(mod.sen, lag.Fattacks~.-lag.Fattacks-lag.Fattacks:lag.states-lag.Hattacks:lag.states)
lm.h <- update(mod.sen, lag.Hattacks~.-lag.Hattacks-lag.Fattacks:lag.states-lag.Hattacks:lag.states)
sort(lm.f$coef, decreasing = TRUE)
sort(lm.h$coef, decreasing = TRUE)

sen1 <- sensemakr(model = mod.sen,
                  treatment = "lag.Fattacks",
                  q=q)
sen2 <- sensemakr(model = mod.sen,
                  treatment = "lag.Hattacks",
                  q=q2)
senF0 <- sensemakr(model = mod.sen,
                   treatment = c("lag.Fattacks"))
senH0 <- sensemakr(model = mod.sen,
                   treatment = c("lag.Hattacks"))

sen.stats <- list("Fatah.symm"=sen1, 
                  "Fatah.null"=senF0, 
                  "Hamas.symm"=sen2, 
                  "Hamas.null"=senH0)

sen.out <- sapply(sen.stats, \(x){with(x$sensitivity_stats, c(rv_q)) })
sen.out <- matrix(sen.out, ncol=2)
sen.out<-rbind(sen.out,
               c(summary(lm.f)$r.sq,
                 summary(lm.h)$r.sq),
               c(max(partial_r2(lm.f)),
                 max(partial_r2(lm.h))))
rownames(sen.out) <- c("RV-symmetry",
                       "RV-null",
                       "$R^2$ other covariates",
                       "Max. partial $R^2$ of other covariates")
colnames(sen.out) <- paste0("\\multicolumn{1}{c}{",
                            c("Fatah", "Hamas"),
                            "}")
partial_r2(mod.sen)["lag.Hattacks"]
sort(partial_r2(lm.h),decreasing = TRUE)[1]


cat(kable(sen.out, digits=2,
             caption="Robustness values and sensitivity of gamma_{i,1}."),
     file="../../Output/Tables/tableD6.txt", sep="\n")


grid <- expand.grid(Eff.Fatah= seq(-1, 1, by=.05),
                    Eff.Popularity = seq(-1.5, 1.5, by=.1))
adj.coef <- cbind(grid,adj.est=0)
for(i in 1:nrow(grid)){
  adj.coef[i,3] <- adjust_coef(mod.sen$coef['lag.Fattacks'], 
                               exposure_confounder_effect = grid[i,"Eff.Fatah"],
                               confounder_outcome_effect=grid[i,"Eff.Popularity"],
                               verbose = FALSE)[1]
}


adj.coef$actor <- "Adjusted~estimate~of~gamma[list(F,1)]"
sensitivity.coef <- ggplot(adj.coef, aes(x = Eff.Fatah, y = Eff.Popularity, z=adj.est)) +
  geom_contour2(aes(color=cut(after_stat(level), c(-Inf, -.1,.2, Inf)),
                    label=after_stat(level),
                    label_color="black",
                    linetype = cut(after_stat(level), c(-Inf, -.1, 0,.2, Inf))),
                breaks=seq(-.4,2.6, by=.2), 
                linewidth=1.3,
                label_size=5, 
                skip=0)+
  scale_color_manual(values=c("navyblue", "orangered", "navyblue")) +
  scale_linetype_manual(values=c("solid", "dotted", "dashed", "solid"))+
  theme_bw(16)+
  guides(linetype="none", color="none")+
  xlab(expression(Marginal~effect~of~an~omitted~variable~on~Pr(italic(a[F]^{t-1}==1))))+
  ylab(expression(Marginal~effect~of~an~omitted~variable~on~italic(Delta~tilde(s)^{~t})))+
  facet_wrap(~actor,labeller = label_parsed)








grid2 <- expand.grid(Eff.Hamas= seq(-1, 1, by=.05),
                     Eff.Popularity = seq(-1.5, 1.5, by=.1))
adj.coef2 <- cbind(grid2,adj.est=0)
for(i in 1:nrow(grid2)){
  adj.coef2[i,3] <- adjust_coef(mod.sen$coef['lag.Hattacks'], 
                                exposure_confounder_effect = grid2[i,"Eff.Hamas"],
                                confounder_outcome_effect=grid2[i,"Eff.Popularity"],
                                verbose = FALSE)[1]
}


adj.coef2$actor <- "Adjusted~estimate~of~gamma[list(H,1)]"
sensitivity.coefH <- ggplot(adj.coef2, aes(x = Eff.Hamas, y = Eff.Popularity, z=adj.est)) +
  geom_contour2(aes(color=cut(after_stat(level), c(-Inf,-1.1, -.9, -.1, 0,Inf)),
                    label=after_stat(level),
                    label_color="black",
                    linetype = cut(after_stat(level), c(-Inf, -1.1,-.9,-.1, 0, Inf))),
                breaks=seq(-1.6,1.4, by=.2), 
                linewidth=1.3,
                label_size=5, 
                skip=0)+
  scale_linetype_manual(values=c("solid", "dashed", "solid", "dotted", "solid"))+
  scale_color_manual(values=c("navyblue", "orangered","navyblue", "orangered", "navyblue")) +
  theme_bw(16)+
  guides(linetype="none", color="none")+
  xlab(expression(Marginal~effect~of~an~omitted~variable~on~Pr(italic(a[H]^{t-1}==1))))+
  ylab("")+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  facet_wrap(~actor,labeller = label_parsed)
figD1 <- arrangeGrob(sensitivity.coef, sensitivity.coefH, nrow=1, ncol=2)
ggsave(figD1, file="../../Output/Figures/figureD1.pdf", width=11.5, height=5)





## question does fatality/attack change over time? (Appendix D.5)
dat$ts <- 1:300
mod4a <- lm(Fkills0~Fattacks*ts, data=dat)
mod4b <- lm(Hkills0~Hattacks*ts, data=dat)

dummies <- with(regData[-c(1:12),], olso.era+2*second.int +3*(1-olso.era-second.int))
dat2 <- cbind(dat, dummies)
mod4c <- lm(Fkills0~Fattacks*factor(dummies), data=dat2)
head(sort(hatvalues(mod4c), decreasing=TRUE))
head(sort(cooks.distance(mod4c), decreasing=TRUE))
mod4d <- lm(Fkills0~Fattacks*factor(dummies), data=dat2[-162])
mod4e <- lm(Hkills0~Hattacks*factor(dummies), data=dat2)


## hypothesis that all three eras are equal
cat("Testing the hypotheses that all three eras are equal\n")
linearHypothesis(mod4e,
                 c("Hattacks:factor(dummies)2=Hattacks:factor(dummies)3",
                   "Hattacks:factor(dummies)3=0"),
                 vcov=NeweyWest)
linearHypothesis(mod4c,
                 c("Fattacks:factor(dummies)2=Fattacks:factor(dummies)3",
                   "Fattacks:factor(dummies)3=0"),
                 vcov=NeweyWest)




mod.list.time <- list(mod4a, mod4c, mod4d, mod4b, mod4e )
##Rename all to "fattacks" for ease of tabling
mod.list.time <- lapply(mod.list.time, \(x){names(x$coefficients)<-gsub(x=names(x$coef), "Hattacks", "Fattacks");return(x)})
se.list.time <- lapply(mod.list.time, function(x){sqrt(diag(NeweyWest(x)))})


tabd7 <- capture.output(stargazer(mod.list.time, 
                                  se=se.list.time,
                                  no.space = TRUE,
                                  omit.stat =  "all",
                                  notes = c("Note:  Newey-West standard errors in parenthesis"),
                                  add.lines=list(c("Actor",
                                                    paste(rep(c("Fatah","Hamas"),c(3,2)),sep="")),
                                                 c("Excluding June 2007 (Battle of Gaza)",
                                                   paste(ifelse(sapply(mod.list.time,
                                                                       \(x){nrow(x$model)==300}
                                                                       ),
                                                         "No", "Yes"),
                                                  sep="")),
                                  c("T",
                                    paste(sapply(mod.list.time, \(m){nrow(m$model)}),
                                    sep="")),
                        c("adj. R2",
                          paste(num2str(sapply(mod.list.time, \(m){summary(m)$adj.r.squared})),
                          sep=""))),
title="Average fatalities per attacks by time",
dep.var.labels = rep(c("Fatalities"), length(mod.list.time)),
digits=2,
order=c("Fattacks"),
align=TRUE,
header=FALSE,
notes.append = FALSE,
star.cutoffs = c(NA),
notes.label = "",
notes.align = "l",
covariate.labels = c("Attacks (count)",
                     "Attacks x time",
                     "Fattacks:factor(dummies)2"="Attacks x Second Intifada -- Election",
                     "Fattacks:factor(dummies)3"="Attacks x Post-2006 election",
                     "Time", 
                     "Second Intifada -- Election",
                     "Post-2006 election",
                     "Constant"),
out="../../Output/Tables/tableD7.txt"))



