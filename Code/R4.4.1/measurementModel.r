#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: Measurement model
#' ---
#'
#' In this file we produce the state variable and related figures and tables
#'
#' Clear workspace and load packages:
library(MARSS)
library(data.table)
library(ggplot2)
library(knitr)
library(zoo)
rm(list=ls())


load("../../Data/actionsSetup.Rdata")
d <- .05
lower <- 0.025
#' Check the factors that measure s
plot.dat <- melt.data.table(dat[,1:7], id.vars = "date", variable.name = "Survey")
plot.dat[, Survey := factor(Survey,
                          levels=c("trustHamas", "trustFatah",
                                   "supportHamas","supportFatah",
                                   "legisHamas", "legisFatah"),
                          labels=c("Trust in Hamas", "Trust in Fatah",
                                   "Support for Hamas", "Support for Fatah",
                                   "Vote for Hamas", "Vote for Fatah"))]
plot.dat[, Survey := paste(Survey, " (N = ", length(na.omit(value)), ")", sep=""), by=Survey]
stateDataPlot <- ggplot(na.omit(plot.dat))+
  geom_point(aes(x=as.Date(date),y=value),size=.5)+
  facet_wrap(~Survey,dir='v', nrow=2)+
  theme_bw(8)+
  xlab("Date")+
    ylab("Percentage")+
    scale_x_date(date_breaks = "10 year", date_labels = "%Y")
ggsave(stateDataPlot,file="../../Output/Figures/figure1.pdf", height=2, width=4)

CorMat <- cor(dat[,2:7], use="pairwise")
rownames(CorMat) <- colnames(CorMat) <- c("Trust in Hamas", "Trust in Fatah",
                                          "Support for Hamas", "Support for Fatah",
                                          "Vote for Hamas", "Vote for Fatah")

cat(kable(CorMat, format="pipe", digits=2,
          caption="Correlations among survey responses"),
    file="../../Output/Tables/tableC2.txt", sep="\n")

#' Build the first stage.
use <- subset(dat, select=c(trustHamas,trustFatah,
                            supportHamas,supportFatah,
                            legisHamas, legisFatah))

X1 <- data.frame(Constant=1,
                 HA=as.numeric(dat$lag.Hattacks>0),
                 FA = as.numeric(dat$lag.Fattacks>0))
useB <- "unconstrained" #\rho
useR <- "identity" # variance of \xi
useQ <- "identity" # variance of \eta


mod <-  list(Z=matrix(names(use)),
             A="zero",
             R=useR,
             B=useB,
             U="zero",
             Q=useQ,
             x0="zero",
             C="unconstrained",
             c = t(X1),
             V0=matrix(5,1,1),
             tinitx=1)
cntl.list = list(maxit=10000, minit=500, abstol = 1e-5, conv.test.slope.tol = 0.05)
#' Fit the measurement model
mars <- MARSS(t(scale(use)), model=mod, control=cntl.list, silent=T)
regData <- data.table(date=dat$date,
                      states = mars$states[1,],
                      lag.Hattacks = mod$c["HA",],
                      lag.Fattacks = mod$c["FA",],
                      Hkills = dat$Hkills, Fkills=dat$Fkills)
#' MARSS latent variable is not sign identified so let's make sure that it moves
#' in the way it should. Observation 153 is the 2006
if(sign(regData$states[153])==1){
  d <-d

  bounds <- quantile( -mars$states[1,], c(lower, 1-lower))
  states <- seq(from=bounds[1], to=bounds[2], by=d)

  G <- expand.grid(states, c(0,1), c(0,1))
  names(G) = c("state", "aH", "aF")
  G <- G[order(G$state, G$aH),]
  rownames(G) <- c()
  nG <- dim(G)[1]
  regData$states <- -1* regData$states
  mars$par$Z <- -1*mars$par$Z
  mars$par$U <- -1*mars$par$U
}else{
  d <- d

  bounds <- quantile( mars$states[1,], c(lower, 1-lower))
  states <- seq(from=bounds[1], to=bounds[2], by=d)

  G <- expand.grid(states, c(0,1), c(0,1))
  names(G) = c("state", "aH", "aF")
  G <- G[order(G$state, G$aH),]
  rownames(G) <- c()
  nG <- dim(G)[1]
}
regData[,lag.states := shift(states)]
save(mars, file= "Results/MARSSoutput.rdata")

marsOut <- as.data.frame(do.call(rbind,mars$par))
marsOut$Equation <- c("Factor Weights ($\\omega$)", rep("", 5),
                      "AR(1) term ($\\rho$)",
                      "Additional inputs ($\\alpha$)", rep("",2))
marsOut$Variable <- c("Trusts Hamas", "Trust Fatah",
                      "Supports Hamas", "Supports Fatah",
                      "Votes Hamas", "Votes Fatah",
                      "Lagged DV",
                      "Constant",
                      "Hamas attack",
                      "Fatah attack")
marsOut <- marsOut[,c("Equation", "Variable","V1")]

cat(kable(marsOut,
      caption="ML estimates for the factor model", format="pipe",
      digits=2,
      row.names = FALSE, col.names=c("Equation", "Variable", "Estimate")),
    file="../../Output/Tables/tableC3.txt", sep="\n")



## states and stats mentioned in paper
cat("When is the state space most friendly to Hamas?\n")
regData[which.min(regData$states)]$date
cat("When is the state space most friendly to Fatah?\n")
regData[which.max(regData$states)]$date

cat("What's the mean, median, and IQR of the state space\n")
round(summary(regData$states), 2)
cat("What's st. dev. of the state space\n")
round(sd(regData$states),2)


##interpretation
cat("How does the state map back onto the trust, support, and vote measures?\n")
print(round(c(trust=lm(I(dat$trustFatah-dat$trustHamas)~regData$states)$coef[2], #map onto trust
        support=lm(I(dat$supportFatah-dat$supportHamas)~regData$states)$coef[2], #map onto support
        vote=lm(I(dat$legisFatah-dat$legisHamas)~regData$states)$coef[2]),1 #map onto vote
      ))

state.plot <- ggplot(regData)+
  geom_line(aes(y=states, x=as.Date(date)))+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
  geom_vline(aes(xintercept=as.Date("1995-11-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("1995-11-30"), y = 12, label = "Oslo II",
           hjust=0, size=3, angle=-90, vjust=.1) +
  geom_vline(aes(xintercept=as.Date("1999-7-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("1999-7-30"), y = 15, label = "Camp David\nsummit fails",
           hjust=0, size=3, angle=-90, vjust=.5) +
  geom_vline(aes(xintercept=as.Date("2001-6-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2001-6-30"), y = 12, label = "Start of several\nmajor Hamas attacks",
           hjust=0, size=3, angle=-90, vjust=.5) +
  geom_vline(aes(xintercept=as.Date("2006-1-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2006-1-30"), y =13 , label = "Hamas wins 2006\nlegislative elections",
           hjust=0, size=3, angle=-90, vjust=.5) +
  geom_vline(aes(xintercept=as.Date("2006-9-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2006-9-30"), y =0 , label = "Hamas-Fatah\nviolence escalates",
           hjust=0, size=3, angle=-90, vjust=.5) +
  geom_vline(aes(xintercept=as.Date("2007-6-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2007-6-30"), y =5 , label = "Battle of Gaza",
           hjust=0, size=3, angle=-90, vjust=.1) +
  geom_vline(aes(xintercept=as.Date("2011-5-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2011-5-30"), y =15 , label = "Fatah and Hamas sign\n1st unity agreement",
           hjust=0, size=3, angle=-90, vjust=.5) +
  geom_vline(aes(xintercept=as.Date("2012-4-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2012-4-30"), y =-4 , label = "Negotiations\nbreakdown",
           hjust=0, size=3, angle=-90, vjust=.5) +
  geom_vline(aes(xintercept=as.Date("2014-6-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2014-6-30"), y =5 , label = "Unity gov't formed",
           hjust=0, size=3, angle=-90, vjust=.1) +
  geom_vline(aes(xintercept=as.Date("2016-2-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2016-2-29"), y =4.5 , label = "Tensions increase\nover Gaza",
           hjust=0, size=3, angle=-90, vjust=.5) +
  theme_bw(14)+
  ylab("Relative popularity (states)")+
  xlab("Date")+
  theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), "points"),
        axis.text.x = element_text(angle = -90))

ggsave(state.plot, file="../../Output/Figures/figure2.pdf", height=4, width=8)

save(regData, states, file="../../Data/measurement.rdata")







rm(list=ls())

## LOAD DATA
jmcc0 <- subset(read.csv("../../Data/jmcc.csv"), year <= 2018)
jmccWB <- subset(read.csv("../../Data/jmcc_WB.csv"), year <= 2018)
jmccGZ <- subset(read.csv("../../Data/jmcc_GAZA.csv"), year <= 2018)
jmccNew <- subset(read.csv("../../Data/jmcc_2023.csv"), year <= 2018)

cpsrWB <- subset(read.csv("../../Data/cpsr_WB.csv"), year <= 2018)
cpsrGZ <- subset(read.csv("../../Data/cpsr_GAZA.csv"), year <= 2018)
cpsr0 <- subset(read.csv("../../Data/cpsr.csv"), year <= 2018)


## SAMPLE SIZES: descriptive statistics for JMCC
tab <- rbind(c(summary(jmccNew$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))])),
             summary(jmccWB$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))]),
             summary(jmccGZ$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))]))
tab <- round(tab)
rownames(tab) <- c("JMCC:All","JMCC:WB", "JMCC:GZ")
print(tab)
## proportion from west bank on average
cat("How many from the West Bank on average?",
    round(mean(jmccWB$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))]/
               jmccNew$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))])*100),
    "%\n")



## SAMPLE SIZES: descriptive statistics for PCPSR
tab <- rbind(c(summary(cpsr0$size[!is.na(cpsr0$supportHamas)])),
             summary(cpsrWB$size[!is.na(cpsr0$supportHamas)]),
             summary(cpsrGZ$size[!is.na(cpsr0$supportHamas)]))
tab <- round(tab)
rownames(tab) <- c("PCPSR:All","PCPSR:WB", "PCPSR:GZ")
print(tab)
cat("How many from the West Bank on average?\n",
    round(mean(cpsrWB$size[!is.na(cpsr0$supportHamas)]/
               cpsr0$size[!is.na(cpsr0$supportHamas)]*100, na.rm=TRUE)),
    "%\n")

cat("How many from the West Bank min/max?\n",
    round(quantile(cpsrWB$size[!is.na(cpsr0$supportHamas)]/
               cpsr0$size[!is.na(cpsr0$supportHamas)]*100, c(0, 1), na.rm=TRUE)),
    "%\n")
