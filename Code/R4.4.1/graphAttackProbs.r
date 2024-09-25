#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: Graph attack probs
#' ---
#' Clear workspace and load packages:

######################################
# packages i need
library(matrixStats)
library(ggplot2)
library(rootSolve)
library(Matrix)
rm(list=ls())

source("helperFunctions.r")
source("gamma2trans.R")

###############################################
# results
load("Results/firstStageOutput.Rdata")
gammaStar <- coefficients(mod0)
sigmaStar <- summary(mod0)$sigma

load("Results/mainModel.rdata")
states <- model.main[["states"]]
Trans <- model.main[["Trans"]]
thetaEst <- model.main[["regtable"]]$V1
delta <-  as.numeric(model.main[["params"]]$delta)
vEst <- model.main[["v"]]
Obs <-  model.main[["regData"]]
EQCP <- AttackProbs(vEst)
mainData <- model.main$regData
mainData$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
mainData$states.int <- sapply(mainData$states.discrete, function(x){which(x==states)})
Tperiod <- dim(mainData)[1]

###############################################
# create states
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()
nG <- dim(G)[1]

# does everything work?
thetaStar = list(betaH = thetaEst[1], # hamas state payoff
                 betaF = thetaEst[2], # fatah state payoff
                 kappaH = c(thetaEst[3], 0), # hamas cost of attack
                 kappaF = c(thetaEst[4], 0), # fatah cost of attack
                 delta=delta)

cat("Checks to make sure that everything works as expected\n")
max(abs(PSI(vEst, thetaStar, Trans, G) - vEst)) < 1e-8
max(abs(Trans - gamma2trans(gammaStar, sigmaStar, states, discretize=F,d = .05, bound=0.025))) < 1e-10
max(abs(PSI(vEst, thetaStar, gamma2trans(gammaStar, sigmaStar, states, discretize=F,d = .05, bound=0.025), G) - vEst)) < 1e-8

###############################################
# Plot estimated attack probabilites as a function of the state

EQCP <- AttackProbs(vEst)
EQIV <- invarDist(vEst, Trans)
vH <- log(colSums(exp(matrix(vEst[1:(2*length(states))],  ncol=length(states)))))
vF <- log(colSums(exp(matrix(vEst[(2*length(states)+1):(4*length(states))],  ncol=length(states)))))

ggdata <-  data.frame(PrAttack = c(EQCP$prAH, EQCP$prAF),
                      states = rep(states,2),
                      ivd = rep(EQIV,2),
                      Vi = c(vH, vF),
                      dVipos = c(diff(vH), NA, diff(vF), NA),
                      dVineg = c(NA, diff(vH[length(vH):1]), NA, diff(vF[length(vF):1])),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)))


pccp <- ggplot(ggdata, aes(x=states, y=PrAttack, color=actor, linetype=actor)) +
    geom_histogram(data=Obs, aes(x=states.discrete, y=after_stat(density)),
                 inherit.aes=F, fill="grey", color="grey", alpha=.75)+
            geom_line(linewidth=1.5) + theme_bw(16) +
            xlab("Relative Popularity (states)") + ylab("Pr. Attack") +
            scale_color_manual(values=c("navyblue", "orangered"),
                     name="Actor") + scale_linetype_manual(values= c(1, 2), name="Actor") +
          geom_rug(data = subset(Obs, Hattacks==1), aes(x = states.discrete), inherit.aes=F, color="orangered1") +
          geom_rug(data = subset(Obs, Fattacks==1), aes(x = states.discrete), inherit.aes=F, color="navyblue")+
  theme(legend.key.width = unit(.5,"in"), legend.position = "bottom") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0))

ggsave("../../Output/Figures/figure4.pdf", pccp, width=8, height=4)

###############################################
# Plot estimated attack probabilites as a function of the OBSERVED state WITH STATE
lostate <- min(mainData$states)
histate <-max(mainData$states)

EQCP <- AttackProbs(vEst)
EQIV <- invarDist(vEst, Trans)
vH <- log(colSums(exp(matrix(vEst[1:(2*length(states))],  ncol=length(states)))))
vF <- log(colSums(exp(matrix(vEst[(2*length(states)+1):(4*length(states))],  ncol=length(states)))))

ggdata <-  data.frame(PrAttack = c(EQCP$prAH[mainData$states.int], EQCP$prAF[mainData$states.int],
                                   (mainData$states-lostate)/(histate-lostate)),
                      time = as.Date(rep(mainData$Date,3), "%Y-%m-%d"),
                      Variable = rep(c("Hamas Pr. Attack", "Fatah Pr. Attack", "Relative Popularity"), each = Tperiod))
ggdata$Variable <- factor(ggdata$Variable, levels = c( "Relative Popularity","Hamas Pr. Attack", "Fatah Pr. Attack"))

pccptime2yd  <- ggplot(ggdata, aes(x=time, y=PrAttack, color=Variable, linetype=Variable, alpha=Variable)) +
  geom_line(linewidth=1.5) + theme_bw(12) +
  xlab("Time") + ylab("Pr. Attack") +
  scale_y_continuous(sec.axis = sec_axis(~.*(histate-lostate) + lostate, name = "Relative popularity")) +
  scale_x_date(breaks=as.Date(paste(seq(from=1994,to=2018, by=4), "-02-01",sep="")),
               labels = seq(from=1994,to=2018, by=4)) +
      scale_color_manual(values=c("grey50", "orangered","navyblue" ),name="Variable")  +
      scale_alpha_manual(values=c(.8,1,1))+
      scale_linetype_manual(values= c(3,2,1), name="Variable") +
  geom_rug(data = subset(mainData, Hattacks==1), aes(x = Date), inherit.aes=F, color="orangered1") +
  geom_rug(data = subset(mainData, Fattacks==1), aes(x = Date), inherit.aes=F, color="navyblue") +
  theme(legend.key.width = unit(.5,"in"), legend.position = "bottom") +
   theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0))
ggsave("../../Output/Figures/figure3.pdf", pccptime2yd, width=8, height=4.75)


cat("Average attack probabilities by actor\n")
round(by(ggdata, ggdata$Variable, \(x){mean(x$PrAttack)}), 2)[-1]
