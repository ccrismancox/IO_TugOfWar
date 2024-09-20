#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Counterfactual on gamma 
#' ---
rm(list=ls())
###############################################
# SET magnitude

magnitude <- 0.01
nstep <- 5

######################################
# packages i need
library("matrixStats")
library("ggplot2")
library("rootSolve")
library("Matrix")
source("gamma2trans.R")
source("helperFunctions.r")
###############################################
# results
load("Results/firstStageOutput.Rdata")
gammaStar <- coefficients(mod0)
sigmaStar <- summary(mod0)$sigma
rm( mod0, states, Trans)


load("Results/mainModel.rdata")
states <- model.main[["states"]]
Trans <- model.main[["Trans"]]
thetaEst <- model.main[["regtable"]]$V1
delta <-  as.numeric(model.main[["params"]]$delta)
vEst <- model.main[["v"]]
Obs <-  model.main[["regData"]]

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
max(abs(Trans - gamma2trans(gammaStar, sigmaStar, states, discretize=F, d=.05, bound=0.025))) < 1e-10
max(abs(PSI(vEst, thetaStar, gamma2trans(gammaStar, sigmaStar, states, discretize=F, d=.05, bound=0.025), G) - vEst)) < 1e-8

EQCP <- AttackProbs(vEst)
EQIV <- invarDist(vEst, Trans)
vH <- log(colSums(exp(matrix(vEst[1:(2*length(states))],  ncol=length(states)))))
vF <- log(colSums(exp(matrix(vEst[(2*length(states)+1):(4*length(states))],  ncol=length(states)))))


###############################################
# Counterfacult: Change Hamas's ability to pull

steps0 <- seq(from=gammaStar[2], to=(1-magnitude)*gammaStar[2], length.out=nstep)
vCF0 <- matrix(NA, nrow=length(vEst), ncol=length(steps0))
vCF0[,1] <- vEst
gamma <- gammaStar 

for (i in 2:length(steps0)){
  gamma[2] <- steps0[i]
  TransNew <- gamma2trans(gamma, sigmaStar, states, discretize=F, d=.05, bound=0.025)
  EQ <- multiroot(function(V){PSI(V, thetaStar, TransNew, G) - V}, vCF0[,i-1],
                               jactype="fullusr",
                               maxiter = 400,
                               jacfunc = function(V){as.matrix(PsiDer(V, thetaStar, TransNew, G)-diag(length(V)))})
  
  vCF0[,i] <- EQ$root
}
ivd0 <- invarDist(vCF0[,length(steps0)],TransNew)

steps1 <- seq(from=gammaStar[2], to=(1+magnitude)*gammaStar[2], length.out=nstep)
vCF1 <- matrix(NA, nrow=length(vEst), ncol=length(steps1))
vCF1[,1] <- vEst
gamma <- gammaStar

for (i in 2:length(steps1)){
  gamma[2] <- steps1[i]
  TransNew <- gamma2trans(gamma, sigmaStar, states, discretize=F, d=.05, bound=0.025) 
  EQ <- multiroot(function(V){PSI(V, thetaStar, TransNew, G) - V}, vCF0[,i-1],
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){as.matrix(PsiDer(V, thetaStar, TransNew, G)-diag(length(V)))})
  
  vCF1[,i] <- EQ$root
}

CFCP0 <- AttackProbs(vCF0[,length(steps0)])
CFCP1 <- AttackProbs(vCF1[,length(steps1)])

diffH0 <- CFCP0$prAH - EQCP$prAH
diffF0 <- CFCP0$prAF - EQCP$prAF

diffH1 <- CFCP1$prAH - EQCP$prAH
diffF1 <- CFCP1$prAF - EQCP$prAF
 

 ggCF_Hgamma <- data.frame(change = c(diffH0[mainData$states.int], diffH1[mainData$states.int],
                                      diffF0[mainData$states.int], diffF1[mainData$states.int]),
                           time = as.Date(rep(mainData$Date,4), "%Y-%m-%d"),
                           actor = rep(c("Hamas", "Fatah"), each=2*Tperiod),
                           cfb = as.factor(rep(c(0, 2, 0, 2), each=Tperiod)), 
                           cfbt =rep(c(paste0("'Decrease magnitude of'~gamma['H,1']~'by'~",magnitude*100,"*'%'"),
                                       paste0("'Increase magnitude of'~gamma['H,1']~'by'~",magnitude*100,"*'%'")), each = Tperiod))
 ggCF_Hgamma$actor <- factor(ggCF_Hgamma$actor, level=c("Fatah", "Hamas"), ordered=T)

 
 
 
###############################################
# Counterfacult: CHange Fatah's ability to pull

steps0 <- seq(from=gammaStar[3], to=(1-magnitude)*gammaStar[3], length.out=nstep)
vCF0 <- matrix(NA, nrow=length(vEst), ncol=length(steps0))
vCF0[,1] <- vEst
gamma <- gammaStar 

for (i in 2:length(steps0)){
  gamma[3] <- steps0[i]
  TransNew <- gamma2trans(gamma, sigmaStar, states, discretize=F, d=.05, bound=0.025) 
  EQ <- multiroot(function(V){PSI(V, thetaStar, TransNew, G) - V}, vCF0[,i-1],
                  jactype="fullusr",
                  maxiter = 500,
                  jacfunc = function(V){as.matrix(PsiDer(V, thetaStar, TransNew, G)-diag(length(V)))})
  
  vCF0[,i] <- EQ$root
}


steps1 <- seq(from=gammaStar[3], to=(1+magnitude)*gammaStar[3], length.out=nstep)
vCF1 <- matrix(NA, nrow=length(vEst), ncol=length(steps1))
vCF1[,1] <- vEst
gamma <- gammaStar

for (i in 2:length(steps1)){
  gamma[3] <- steps1[i]
  TransNew <- gamma2trans(gamma, sigmaStar, states, discretize=F, d=.05, bound=0.025) 
  EQ <- multiroot(function(V){PSI(V, thetaStar, TransNew, G) - V}, vCF0[,i-1],
                  jactype="fullusr",
                  maxiter = 400,
                  jacfunc = function(V){as.matrix(PsiDer(V, thetaStar, TransNew, G)-diag(length(V)))})
  
  vCF1[,i] <- EQ$root
}

CFCP0 <- AttackProbs(vCF0[,length(steps0)])
CFCP1 <- AttackProbs(vCF1[,length(steps1)])

diffH0 <- CFCP0$prAH - EQCP$prAH
diffF0 <- CFCP0$prAF - EQCP$prAF

diffH1 <- CFCP1$prAH - EQCP$prAH
diffF1 <- CFCP1$prAF - EQCP$prAF

ggCF_Fgamma <- data.frame(change = c(diffH0[mainData$states.int], diffH1[mainData$states.int],
                                     diffF0[mainData$states.int], diffF1[mainData$states.int]),
                          time = as.Date(rep(mainData$Date,4), "%Y-%m-%d"),
                          actor = rep(c("Hamas", "Fatah"), each=2*Tperiod),
                          cfb = as.factor(rep(c(0, 2, 0, 2), each=Tperiod)), 
                          cfbt =rep(c(paste0("'Decrease magnitude of'~gamma['F,1']~'by'~",magnitude*100,"*'%'"),
                                      paste0("'Increase magnitude of'~gamma['F,1']~'by'~",magnitude*100,"*'%'")), each = Tperiod))
ggCF_Fgamma$actor <- factor(ggCF_Fgamma$actor, level=c("Fatah", "Hamas"), ordered=T)



###############################################
# final figure

useLevel <- c("i=='Hamas'", "i=='Fatah'")
ggCF_together <- rbind(ggCF_Hgamma, ggCF_Fgamma)
ggCF_together$CF <- rep(useLevel, each = dim(ggCF_Hgamma)[1])
ggCF_together$CF<- factor(ggCF_together$CF, level=useLevel, ordered=T)



plotCF_together <- ggplot(ggCF_together, aes(x=time, y=change, color=actor)) +
  geom_path(size=1.05) +
  theme_bw(12) +
  xlab("Time") + ylab("Change in Pr. Attack") + 
  facet_wrap(~cfbt, labeller=label_parsed,nrow = 2)+
  scale_color_manual("Effect on",
                     values=c("navyblue", "orangered")) +
  scale_linetype_manual(values= c(1, 2), name="Actor") +
  geom_hline(yintercept=0,size=0.75) +
  scale_x_date(breaks=as.Date(paste(seq(from=1994,to=2018, by=4), "-02-01",sep="")),
               labels = seq(from=1994,to=2018, by=4)) +
  theme(legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0)) 

# sATE 
# needs to be added to plotCF_together
UL <- ggCF_together[ggCF_together$cfb==0 & ggCF_together$CF == useLevel[1],]
UR <- ggCF_together[ggCF_together$cfb==2 & ggCF_together$CF == useLevel[1],]
BL <- ggCF_together[ggCF_together$cfb==0 & ggCF_together$CF == useLevel[2],]
BR <- ggCF_together[ggCF_together$cfb==2 & ggCF_together$CF == useLevel[2],]

sATE <- data.frame(upperleft = c(mean(UL$change[UL$actor=="Hamas"]),
                                 sd(UL$change[UL$actor=="Hamas"]), 
                                 mean(UL$change[UL$actor=="Fatah"]),
                                 sd(UL$change[UL$actor=="Fatah"])),
                   upperright = c(mean(UR$change[UR$actor=="Hamas"]),
                                  sd(UR$change[UR$actor=="Hamas"]), 
                                  mean(UR$change[UR$actor=="Fatah"]),
                                  sd(UR$change[UR$actor=="Fatah"])),
                   bottomleft = c(mean(BL$change[BL$actor=="Hamas"]),
                                  sd(BL$change[BL$actor=="Hamas"]), 
                                  mean(BL$change[BL$actor=="Fatah"]),
                                  sd(BL$change[BL$actor=="Fatah"])),
                   bottomright = c(mean(BR$change[BR$actor=="Hamas"]),
                                   sd(BR$change[BR$actor=="Hamas"]), 
                                   mean(BR$change[BR$actor=="Fatah"]),
                                   sd(BR$change[BR$actor=="Fatah"])))

sATE.df <- data.frame(sATE= paste("Hamas: ", signif(sATE[1,], 2),"\nFatah: ",signif(sATE[3,],2)),
                      CF=rep(useLevel,each=2),
                      actor=NA,
                      cfbt =c(paste0("'Decrease magnitude of'~gamma['H,1']~'by'~",magnitude*100,"*'%'"),
                              paste0("'Increase magnitude of'~gamma['H,1']~'by'~",magnitude*100,"*'%'"),
                              paste0("'Decrease magnitude of'~gamma['F,1']~'by'~",magnitude*100,"*'%'"),
                              paste0("'Increase magnitude of'~gamma['F,1']~'by'~",magnitude*100,"*'%'")),
                      title="bold(sATE)")

plotCF_together <- plotCF_together+
  geom_text(data=sATE.df, aes(x=ggCF_together$time[1], y=Inf,label=title, hjust=0,vjust=1.7),
            size=3, parse = TRUE,
            color="Black")+
  geom_text(data=sATE.df, aes(x=ggCF_together$time[1], y=Inf,label=sATE, hjust=0,vjust=1.8),
            size=3,
            color="Black") 

ggsave("../../Output/Figures/Figure6.pdf", plot = plotCF_together,
       width = 7.75, height = 5.5, units = "in")





