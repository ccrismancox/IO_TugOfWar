#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Model fit
#' ---
rm(list=ls())
library(data.table)
library(knitr)
library(moments)
load("Results/mainStateActions.Rdata")
load("Results/firstStageOutput.Rdata")
source("gamma2trans.R")
source("helperFunctions.r")

load("Results/mainModel.rdata")
load("Results/noCompetition.Rdata")
load("Results/t4t_Model.rdata")

source("helperFunctions.r")
eras <- factor(rep(1:3, c(81,64,155)),
               labels=c("Oslo", "2nd Intifada", "Post-election"))




#### Main E[% correctly predicted] ####
model.main$regData$ts <- 1:300
state.probs <- cbind(data.table(states=model.main$states),
                     AttackProbs(model.main$v))
stateProbs <- merge(model.main$regData,
                    state.probs, by.y="states", by.x="states.discrete", 
                    all.y=FALSE, all.x=TRUE)
stateProbs <- stateProbs[order(ts)]
ECP.main <- colMeans(with(stateProbs,
                          cbind(Hattacks, Fattacks)*cbind(prAH, prAF)+
                            (1-cbind(Hattacks, Fattacks))*(1-cbind(prAH, prAF))))
ECP.main.byEras <- do.call(rbind, by(with(stateProbs,
                           cbind(Hattacks, Fattacks)*cbind(prAH, prAF)+
                             (1-cbind(Hattacks, Fattacks))*(1-cbind(prAH, prAF))),
                      eras, colMeans))



ECP.main2 <- mean(with(stateProbs,
                       (Hattacks*Fattacks) * (prAH*prAF)+
                         ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                         (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                         ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF)))

ECP.main2.byEras <- c(by(with(stateProbs,
                            (Hattacks*Fattacks) * (prAH*prAF)+
                              ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                              (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                              ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF)),
                       eras, mean))

#### T4T E[% correctly predicted] ####
model.main.t4t$regData$ts <- 1:300
state.probs <- cbind(data.table(states=model.main.t4t$states),
                     AttackProbs(model.main.t4t$v))
stateProbs <- merge(model.main.t4t$regData, state.probs, 
                    by.y="states", by.x="states.discrete", 
                    all.y=FALSE, all.x=TRUE)
stateProbs <- stateProbs[order(ts)]
ECP.t4t <- colMeans(with(stateProbs,
                         cbind(Hattacks, Fattacks)*cbind(prAH, prAF)+
                           (1-cbind(Hattacks, Fattacks))*(1-cbind(prAH, prAF))))

ECP.t4t.byEras <- do.call(rbind, by(with(stateProbs,
                          cbind(Hattacks, Fattacks)*cbind(prAH, prAF)+
                            (1-cbind(Hattacks, Fattacks))*(1-cbind(prAH, prAF))),
                     eras,
                     colMeans))

ECP.t4t2 <- mean(with(stateProbs,
                      (Hattacks*Fattacks) * (prAH*prAF)+
                        ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                        (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                        ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF)))

ECP.t4t2.byEras <- c(by(with(stateProbs,
                           (Hattacks*Fattacks) * (prAH*prAF)+
                             ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                             (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                             ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF)),
                      eras, mean))


#### No contest E[% correctly predicted] ####

Pout <- AttackProbs(noCompModel$v)
state.probs <- data.table(cbind(states,Pout))

stateProbsNoComp <- merge(mainData, state.probs,
                          by.y="states", by.x="states.discrete",  
                          all.y=FALSE, all.x=TRUE, sort=FALSE)
ECP.NoComp <- colMeans(with(stateProbsNoComp,
                            cbind(Hattacks, Fattacks)*cbind(prAH, prAF)+
                              (1-cbind(Hattacks, Fattacks))*(1-cbind(prAH, prAF))))
ECP.NoComp.byEras <- do.call(rbind, by(with(stateProbsNoComp,
                            cbind(Hattacks, Fattacks)*cbind(prAH, prAF)+
                              (1-cbind(Hattacks, Fattacks))*(1-cbind(prAH, prAF))),
                        eras,colMeans))


ECP.NoComp2 <- mean(with(stateProbsNoComp,
                       (Hattacks*Fattacks) * (prAH*prAF)+
                         ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                         (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                         ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF)))

ECP.NoComp2.byEras <- c(by(with(stateProbsNoComp,
                         (Hattacks*Fattacks) * (prAH*prAF)+
                           ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                           (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                           ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF)),
                         eras, mean))


# LR test
## What we have is not the full likelihood
partial.Likelihood.main <- -model.main$conv$V1[2]
partial.Likelihood.noComp <- -noCompModel$conv$V1[2] 

## So we need to include the gamma term
gamma2likelihood_term <- function(Z, gamma, sigma0, Y, given,states, sum=TRUE){
  mu <- Z %*% gamma
  pr <- pr0 <- pr1 <- rep(1, 300)
  d <- diff(states)[3]/2
  
  
  for(i in 1:299){
    hi <- pnorm((Y[2:nrow(Y)]$states.discrete[i]+d-mu[i])/sigma0)
    lo <- pnorm((Y[2:nrow(Y)]$states.discrete[i]-d-mu[i])/sigma0)
    pr[i+1] <- hi-lo
    pr0[i+1] <- hi
    pr1[i+1] <- 1-lo
  }
  lowest <- Y[,1]==min(given$stateActions[,1])
  highest <- Y[,1]==max(given$stateActions[,1])
  trans <- pr*(1-lowest)*(1-highest) + pr0*(lowest)*(1-highest) + pr1*(1-lowest)*(highest)
  if(sum){
    return(sum(log(trans)))
  }else{
    return(log(trans))
  }
  
}




Z <- with(mainData, cbind(1,Hattacks, Fattacks, states.discrete, 
                          Hattacks*states.discrete, Fattacks*states.discrete))
Y <- mainData[,c("states.discrete", "Hattacks", "Fattacks")]
gamma <- mod0$coefficients
given <- genGiven(states, mainData)
sigma0 <- summary(mod0)$sigma
gammaTerm.main <- gamma2likelihood_term(Z, gamma, sigma0, Y, given, states)
FULL.Likelihood.main <- -model.main$conv$V1[2] + gammaTerm.main


# No Competition
modRW <- lm(states~lag.states, data=regData)
gamma <- c(modRW$coef[1],0,0, modRW$coef[2],0,0)
sigma <- summary(modRW)$sigma
gammaTerm.NoComp <- gamma2likelihood_term(Z,
                                          gamma=gamma,
                                          sigma, Y, given, states)
FULL.Likelihood.NoComp <- -noCompModel$conv$V1[2] + gammaTerm.NoComp



LRstat.noComp<- -2*(FULL.Likelihood.NoComp-FULL.Likelihood.main)
full.pval <-    pchisq(LRstat.noComp, lower=FALSE, df=6)



out <- rbind(cbind(ECP.main,ECP.NoComp, ECP.t4t),
      c(ECP.main2,ECP.NoComp2, ECP.t4t2))
rownames(out) <- c("ePCP--Hamas",
                   "ePCP--Fatah",
                   "ePCP--Overall")
colnames(out) <- c("Main model", "No competition", "Tit-for-tat")

cat(kable(out,
          caption="In-sample model fit",
          digits=2,
          format="pipe"),
    file="../../Output/Tables/tableG2.txt",
    sep="\n")
      




## Comparing t4t to the main man
LL_pointwise <- function(p, mod){
  
  state.probs <- cbind(data.table(states=mod$states),
                       data.table(prAH=p[1:length(mod$states)],
                                  prAF=p[(length(mod$states)+1):length(p)]))
  mod$regData$ts <- 1:300
  stateProbs <- merge(mod$regData,
                      state.probs, by.y="states", by.x="states.discrete", all.y=FALSE, all.x=TRUE,
                      sort=FALSE)
  
  return( log(with(stateProbs,
                   (Hattacks*Fattacks) * (prAH*prAF)+
                     ((1-Hattacks)*Fattacks) * ((1-prAH)*prAF)+
                     (Hattacks*(1-Fattacks)) * (prAH*(1-prAF))+
                     ((1-Hattacks)*(1-Fattacks)) * (1-prAH)*(1-prAF))))
}

n <- 300
p1 <- unlist(AttackProbs(model.main$v))
p2 <- unlist(AttackProbs(model.main.t4t$v))
## positive point wise are better for the main model
point.wise <- LL_pointwise(p1, model.main)-LL_pointwise(p2, model.main.t4t)



## Excess kurtosis is kurtosis -3; positive values are leptokurtic and speak well for the Clarke test
kurtosis(point.wise) -3 
skewness(point.wise) #skew in -1 to 1 is good too (not skewed )


## Clarke test with the null the alternative hypothesis that model 1 is better
Clarke.stat <- sum(point.wise>0)
pC <- pbinom(Clarke.stat, size=300, prob = .5, lower=FALSE ) 




#### Table 3####
cat(kable(data.table(`Alternative model`=c("No competition", "Tit-for-tat"),
                     Test=c("Likelihood ratio", "Clarke's test"),
                     `Null distribution`=c("chi^2(6)", "Binomial(300, 0.5)"),
                     Statistic=paste(round(c(LRstat.noComp, Clarke.stat))),
                     `p value`=ifelse(c(full.pval, pC) < 0.01, "< 0.01", "> 0.01")),
          caption="Comparative model tests"),
    file="../../Output/Tables/table3.txt",
    sep="\n")




X <- num2str(cbind(t(rbind(model.main$regtable[1:2,1:2],
                           matrix(NA, 2,2), 
                           model.main$regtable[3:4,1:2])),
                   t(rbind(matrix(NA, 4,2), noCompModel$regtable[,1:2])),
                   t(rbind(matrix(NA, 2,2), model.main.t4t$regtable[,1:2]))
                   ))
X[X==" NA"] <- NA
X[2,] <- ifelse(!is.na(X[2,]),paste0("(",X[2,], ")"),NA)
X <- matrix(X, ncol=3)
X <- rbind(X, num2str(-c(model.main$conv$V1[2], 
                         noCompModel$conv$V1[2],
                         model.main.t4t$conv$V1[2])), 
           c(300,300,300))
X <- cbind(c(rbind(c("beta_H","beta_F",
                     "tau_H", "tau_F",
                     "kappa_H", "kappa_F"), ""),
             "LL", "T"),
           X)
colnames(X) <- c(" ",
                 "Outbidding (main model)", 
                 "No-competition",
                 "Tit-for-tat")
X[is.na(X)] <- "  "
cat(kable(X,
          caption="Comparing outbidding to other theories.",
          format="pipe"),
    file="../../Output/Tables/tableG1.txt",
    sep="\n"
)



######### The Eras Tour  ###### 
## comparing to t4t
clark.byEras <- do.call(rbind,
                        by(point.wise, eras, 
                           \(x){c(length(x),sum(x>0),
                                  pbinom(sum(x>0),
                                         size=length(x), 
                                         prob = .5, lower=FALSE ))}))





## comparing to no competition
p1 <- unlist(AttackProbs(model.main$v))
p2 <- unlist(AttackProbs(noCompModel$v))





#### add in the first stage
Z <- with(mainData, cbind(1,Hattacks, Fattacks, states.discrete, 
                          Hattacks*states.discrete, Fattacks*states.discrete))
Y <- mainData[,c("states.discrete", "Hattacks", "Fattacks")]
gamma <- mod0$coefficients
given <- genGiven(states, mainData)
sigma0 <- summary(mod0)$sigma
gammaTerm.main <- gamma2likelihood_term(Z, gamma, sigma0, Y, given, states,sum = FALSE)
FULL.Likelihood.main <- LL_pointwise(p1, model.main)+ gammaTerm.main


# No Competition
modRW <- lm(states~lag.states, data=regData)
gamma <- c(modRW$coef[1],0,0, modRW$coef[2],0,0)
sigma <- summary(modRW)$sigma
gammaTerm.NoComp <- gamma2likelihood_term(Z,
                                          gamma=gamma,
                                          sigma, Y, given, states, sum=FALSE)
FULL.Likelihood.NoComp <- LL_pointwise(p2, noCompModel)+ gammaTerm.NoComp


erasLL <- by(FULL.Likelihood.main, eras, sum)
erasLL.noComp <- by(FULL.Likelihood.NoComp, eras, sum)


LRstat.noComp<- c(-2*(erasLL.noComp-erasLL))
full.pval <-    c(pchisq(LRstat.noComp, lower=FALSE, df=6))

## eras table 
out <- cbind.data.frame(clark.byEras, LRstat.noComp, full.pval)
colnames(out) <- c("Obs", "Clarke stat", "Clarke $p$ value", "LR stat", "LR $p$ value")
out[,c(1)] <- paste(round(out[,c(1)] ))
out[,c(2)] <- paste(round(out[,c(2)] ))
out[,3] <- round(out[,3], 2)
out[out[,3] < 0.05,][,3] <- "$ < 0.05$"
out[,c(4)] <- paste(round(out[,c(4)] ))
out[out[,5] < 0.05,][,5] <- "$ < 0.05$"
out <- cbind(Era=unique(eras), out)

cat(kable(out,
          caption="Sub-sample model fit",
          format="pipe"),
    file="../../Output/Tables/tableG3.txt",
    sep="\n")

