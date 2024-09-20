#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: 3 Feb 2022
#' title: First stage regressions
#' header-includes:
#'     - \usepackage{rotating}
#'     - \usepackage{dcolumn}
#' output: pdf_document
#' ---
#' Clear workspace and load packages:
rm(list=ls())
library(data.table)
library(vars)
library(urca)
library(knitr)
library(lmtest)
library(sandwich)
library(ggplot2)

load("../../Data/measurement.rdata")
load("../../Data/actionsSetup.Rdata")

source("helperFunctions.r")
regData$Hattack <- as.numeric(dat$Hattacks>0)
regData$Fattack <- as.numeric(dat$Fattacks>0)
regData$Hattack.count <- (dat$Hattack)
regData$Fattack.count <- (dat$Fattack)
regData$ln.Hattack <- log(dat$Hattack + 1)
regData$ln.Fattack <- log(dat$Fattacks+1)

regData$attacks <- dat$Hattacks+dat$Fattacks
regData[,state.deviations := -1*abs(scale(states))]
regData[, diff.states:= states -lag.states]
regData[,lag.Hattack.count := shift(Hattack.count)]
regData[,lag.Fattack.count := shift(Fattack.count)]
regData[,lag.attacks := shift(attacks)]


## deviations version mentioned in the text ##
coeftest(glm(attacks~state.deviations, family=poisson, data=regData),
         NeweyWest)

##### VAR approach#####
## look for stationarity
coint_ca.joA <- ca.jo(regData[,c("diff.states", "Hattack", "Fattack")],
                      type = "trace",
                      spec = "transitory")
summary(coint_ca.joA) # every null rejected so VAR is fine with these variables
## using BIC or SC as it's listed here
VARselect( na.omit(regData[,c("diff.states", "Hattack", "Fattack")]), lag.max = 10)$selection

V1 <- VAR( na.omit(regData[,c("diff.states", "Hattack", "Fattack")]), p =2)




n.ahead=6
IRF <- irf(V1, n.ahead=n.ahead, ortho=FALSE)
irf.dat <- data.frame(est=c(sapply(IRF$irf, c)),
           lo=c(sapply(IRF$Lower, c)),
           hi=c(sapply(IRF$Upper, c)),
           Shock=rep(c("Change in state", "Hamas attacks", "Fatah attacks"), each=3*(n.ahead+1)),
           Period=0:n.ahead,
           Outcome=rep(c("Change in state", "Pr(Hamas attacks)", "Pr(Fatah attacks)"), each=n.ahead+1))
           
irf.plot <- ggplot(irf.dat) + 
  geom_line(aes(x=Period, y=est))+
  geom_ribbon(aes(x=Period, ymin=lo, ymax=hi), alpha=.2)+
  facet_grid(Outcome~Shock,scales = "free_y" , switch="y")+
  theme_bw(13)+
  ylab("Marginal effect")+
  xlab("Months after shock")+
  geom_hline(aes(yintercept=0), linetype="dashed", color="red")
ggsave("../../Output/Figures/figureG1.pdf", irf.plot, height=5, width=7)

VAR.out <- rbind(summary(V1)$varresult$diff.states$coef,
                 summary(V1)$varresult$Hattack$coef,
                 summary(V1)$varresult$Fattack$coef)
VAR.out <- t(num2str(VAR.out[,1:2]))
VAR.out[2,] <- paste0("(",VAR.out[2,],")")
VAR.out <- matrix(VAR.out, ncol=3)
VAR.out  <-cbind( c("Delta state_{t-1}", "",
                    "Hamas attacks_{t-1}","",
                    "Fatah attacks_{t-1}","",
                     "Delta state_{t-2}", "",
                     "Hamas attacks_{t-2}","",
                     "Fatah attacks_{t-2}","",
                    "Const.", ""
                    ),
                  VAR.out)
colnames(VAR.out) <- c(" ",
                       "Delta states",
                        "Hamas Attack",
                        "Fatah Attack")
tab.out <- VAR.out[c(3:4, 9:10,
                     5:6, 11:12, 
                     1:2, 7:8,
                     (nrow(VAR.out)-1):nrow(VAR.out)),]
cat(kable(VAR.out, format="pipe", digits=2,
          caption="VAR regression results"),
    file="../../Output/Tables/tableG4.txt", sep="\n")


