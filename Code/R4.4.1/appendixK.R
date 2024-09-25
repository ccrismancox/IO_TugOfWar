#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title:  Interpreting betas
#' ---


######################################
# packages i need
library(ggplot2)
library(rootSolve)
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


max(abs(PSI(vEst, thetaStar, Trans, G) - vEst)) < 1e-8
max(abs(Trans - gamma2trans(gammaStar, sigmaStar, states, discretize=F,d = .05, bound=0.025))) < 1e-10
max(abs(PSI(vEst, thetaStar, gamma2trans(gammaStar, sigmaStar, states, discretize=F,d = .05, bound=0.025), G) - vEst)) < 1e-8

###############################################
# Plot estimated attack probabilites as a function of the state

EQCP <- AttackProbs(vEst)

JV <- PsiDer(vEst, thetaStar, Trans, G)-diag(length(vEst))
JBi <- cbind(matrix(c(rep(states, each =2), rep(0, length(states)*2)), ncol=1),
            matrix(c(rep(0, length(states)*2), rep(states, each =2)), ncol=1))
TJ <- - solve(JV) %*%  JBi  
predicted <- vEst -  matrix(c(thetaStar$betaH, thetaStar$betaF),1,2) %*% t(TJ)

thetaNew <- thetaStar
thetaNew$betaF <- thetaNew$betaH <- 0
EQnobeta <- multiroot(function(V){PSI(V, thetaNew, Trans, G) - V}, as.numeric(predicted),
                      jactype="fullusr",
                maxiter = 400, jacfunc = function(V){PsiDer(V, thetaNew, Trans, G)-diag(length(V))},
                rtol=1e-12, atol=1e-12, ctol = 1e-12)
EQnb <- AttackProbs(EQnobeta$root)

ggdata <-  data.frame(PrAttack = c(EQCP$prAH[mainData$states.int], 
                                  EQCP$prAF[mainData$states.int]),
                     time = as.Date(rep(mainData$Date,2), "%Y-%m-%d"),
                     actor = rep(c("Hamas", "Fatah"), each = Tperiod))

size <- 1.5
plotU <- ggplot(ggdata, aes(x=time, y=PrAttack, color=actor,linetype=actor)) +
  geom_path(linewidth=size) +
  theme_bw(16) + ylim(c(0,0.65)) + 
  xlab("Time") + ylab("Pr. Attack") +  
  scale_color_manual("Actor", values=c("navyblue", "orangered")) +
  scale_linetype_manual(values= c(1, 2), name="Actor") +
  scale_x_date(breaks=as.Date(paste(seq(from=1994,to=2018, by=4), "-02-01",sep="")),
               labels = seq(from=1994,to=2018, by=4)) +
  theme(legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0)) + 
  geom_hline(yintercept=unique(round(EQnb$prAH, digits=10)), linetype=2, color="orangered", linewidth=size,alpha=0.4) + 
  geom_hline(yintercept=unique(round(EQnb$prAF, digits=10)), linetype=1, color="navyblue", linewidth=size,alpha=0.4) +
  annotate("text", x = as.Date("2016-09-01", "%Y-%m-%d"), y = 0.06, label = paste0("'Assuming'~beta[F]==",0),parse = TRUE) +
  annotate("text", x = as.Date("2016-09-01", "%Y-%m-%d"), y = 0.26, label = paste0("'Assuming'~beta[H]==",0),parse = TRUE)

ggsave("../../Output/Figures/figureK1.pdf", plot = plotU,
       width = 7.75, height = 5.25, units = "in")
