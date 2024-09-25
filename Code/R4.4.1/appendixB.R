#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Numerical example
#' ---

######################################
# packages i need
library("matrixStats")
library("ggplot2")
library("rootSolve")
library("Matrix")

rm(list=ls())

source("helperFunctions.r")
source("gamma2trans.R")
######################################
# functions i need
shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else
    out<-x
  out
}

###############################################
# set up model
states <- seq(from = -50, to = 50, by = 1) 
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()
nG <- dim(G)[1]

thetaStar = list(betaH = -1/500, # hamas state payoff 
                 betaF = 1/500, # fatah state payoff
                 kappaH = c(-2, 0), # hamas cost of attack
                 kappaF = c(-2, 0), # fatah cost of attack
                 delta=0.999)  

gammaStar <- c(0, -1, 1, 1, 0, 0)
sigmaStar <- 2

# transition matrix
Trans <- gamma2trans(gammaStar, sigmaStar, states, d=1, discretize=F)

# per period utilities
uas <- cbind(thetaStar$betaH*G$state + (thetaStar$kappaH[1] + thetaStar$kappaH[2]*G$state)*G$aH,
             thetaStar$betaF*G$state + (thetaStar$kappaF[1] + thetaStar$kappaF[2]*G$state)*G$aF)
colnames(uas) <- c("uH", "uF")

###############################################
# the symmetric eq
vstart <- c(rowMeans(matrix(uas[,'uH'], ncol=2, byrow=T)),
            c(apply(array(uas[,'uF'], c(2,2,length(states))), 3, rowMeans)))

thetaStar$delta <- 0.90
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, vstart, 
                   maxiter = 100, jactype = "fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.99
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
          maxiter = 100, jactype = "fullusr",
          jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.9945
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 200, jactype = "fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.99675
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 200, jactype = "fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

thetaStar$delta <- 0.999
EQsym <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, EQsym$root, 
                   maxiter = 200, jactype = "fullusr",
                   rtol=1e-10, atol=1e-12, ctol=1e-12,
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})

EQCPsym <- AttackProbs(EQsym$root)
ggdatasym <-  data.frame(PrAttack =c(EQCPsym$prAH, EQCPsym$prAF),
                          id = invarDist(EQsym$root,Trans),
                      states = rep(states,2),
                      dPjds = c(0.5*shift(EQCPsym$prAF,1) - 0.5*shift(EQCPsym$prAF,-1),
                                0.5*shift(EQCPsym$prAH,1) - 0.5*shift(EQCPsym$prAH,-1)),
                actor = rep(c("Hamas", "Fatah"), each = length(states)))



###############################################
# the non-symmetric eq

set.seed(1)
err <- rnorm(length(vstart))
EQ1 <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, 
                 vstart+err, 
                   maxiter = 500, rtol=1e-10, atol=1e-12, ctol=1e-12,
                 jactype = "fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})
EQCP1 <- AttackProbs(EQ1$root)
ggdata1 <-  data.frame(PrAttack =c(EQCP1$prAH, EQCP1$prAF),
                      states = rep(states,2),
                      id = invarDist(EQ1$root,Trans),
                      dPjds = c(0.5*shift(EQCP1$prAF,1) - 0.5*shift(EQCP1$prAF,-1),
                                0.5*shift(EQCP1$prAH,1) - 0.5*shift(EQCP1$prAH,-1)),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)))


vHamas <- matrix(EQ1$root[1:(2*length(states))],  ncol=length(states)) 
vFatah <- matrix(EQ1$root[(2*length(states)+1):(4*length(states))],  ncol=length(states))

start3 <- c(vFatah[,length(states):1], vHamas[,length(states):1])
EQ3 <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, start3, 
                 maxiter = 500, rtol=1e-10, atol=1e-12, ctol=1e-12,
                 jactype = "fullusr",
                 jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))})
EQCP3 <- AttackProbs(EQ3$root)
ggdata3 <-  data.frame(PrAttack =c(EQCP3$prAH, EQCP3$prAF),
                      states = rep(states,2),
                      id = invarDist(EQ3$root,Trans),
                      dPjds = c(0.5*shift(EQCP3$prAF,1) - 0.5*shift(EQCP3$prAF,-1),
                                0.5*shift(EQCP3$prAH,1) - 0.5*shift(EQCP3$prAH,-1)),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)))



###############################################
# compare all 3 eq

ggAll <- rbind(ggdata3, ggdatasym, ggdata1)
ggAll$EQ <- rep(c("Fatah-Dominant Eq.", "Symmetric Eq.", "Hamas-Dominant Eq."), each=dim(ggdata1)[1])
ggAll$EQ <- factor(ggAll$EQ, levels= c("Fatah-Dominant Eq.", "Symmetric Eq.", "Hamas-Dominant Eq."), ordered=T)

Trans.symm <- Trans


pEQall <- ggplot(ggAll, aes(x=states, y=PrAttack, color=actor, linetype=actor)) + 
    geom_line(linewidth=1.25) + facet_grid(cols=vars(EQ)) + 
  theme_bw(11) + xlab("Relative Popularity (states)") + ylab("Pr. Attack") + theme(legend.position=c(0.91, 0.8)) +
  scale_color_manual(values=c("navyblue", "orangered"),
  name="Actor") + scale_linetype_manual(values= c(1, 2), name="Actor") +
  theme(legend.text=element_text(size=rel(0.6)),
        legend.title=element_text(size=rel(0.8)),
        legend.key.size = unit(0.5, 'lines')) + 
  theme(legend.background = element_rect(fill="white",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="grey50"))

pIDall <- ggplot(subset(ggAll, actor=="Hamas"),  aes(x=states, y=id)) + 
  geom_bar(stat="identity", color = "seagreen") + facet_grid(cols=vars(EQ)) + 
  theme_bw(11) + xlab("Relative Popularity (states)") + ylab("Invariant Dist.") 

pCorall <- ggplot(subset(ggAll, states <40 & states>-40), aes(x=dPjds, y=PrAttack))+ 
      geom_point(size=1.25) + facet_grid(cols=vars(EQ), rows=vars(actor), scales = "free") +
      theme_bw(11) +  geom_smooth(method='lm', formula= y~x )

ggsave("../../Output/Figures/figureB3.pdf" , height = 2.5, plot = pIDall, width= 6)
ggsave("../../Output/Figures/figureB2.pdf" , height = 2.5, plot = pEQall, width= 6)

## correlations mentioned in the appendix
cat("Correlations between choice probs in \n")
cat("Hamas dominant eq.", round(cor(EQCP1$prAH, EQCP1$prAF), 2), "\n")
cat("Fatah dominant eq.", round(cor(EQCP3$prAH, EQCP3$prAF), 2), "\n")
cat("Sym. dominant eq.", round(cor(EQCPsym$prAH, EQCPsym$prAF), 2), "\n")






###############################################
# comparative statics on eq
# Counterfact: Hamas care more about popularity
EQall <- cbind(EQsym$root, EQ1$root, EQ3$root)

steps0 <- seq(from=thetaStar$betaH, to = 0.5*thetaStar$betaH, length.out=100)
steps1 <- seq(from=thetaStar$betaH, to = 1.5*thetaStar$betaH, length.out=100)
vCF0 <- array(NA, c(length(EQ1$root), length(steps0), 3))
vCF1 <- array(NA, c(length(EQ1$root), length(steps1), 3))

for (e in 1:dim(EQall)[2]){
  vCF0[,1,e] <- EQall[,e]
  theta <- thetaStar 
  for (i in 2:length(steps0)){
    JV <- PsiDer(vCF0[,i-1,e], theta, Trans, G)-diag(length(vCF0[,i-1,e]))
    JBH <- matrix(c(rep(states, each =2), rep(0, length(states)*2)), nrow=1)  
    TJ <- - JBH %*% solve(t(JV))
    predicted <- vCF0[,i-1,e] + as.numeric((steps0[i] - steps0[i-1]) * TJ)
    
    theta$betaH <- steps0[i]
    EQcf <- multiroot(function(V){PSI(V, theta, Trans, G) - V},
                      predicted, maxiter = 400, rtol=1e-10, atol=1e-12, ctol=1e-12,
                      jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))},
                      jactype="fullusr")
    vCF0[,i,e] <- as.numeric(EQcf$root)
  }
  
  vCF1[,1,e] <- EQall[,e]
  theta <- thetaStar 
  for (i in 2:length(steps1)){
    JV <- PsiDer(vCF1[,i-1,e], theta, Trans, G)-diag(length(vCF1[,i-1,e]))
    JBH <- matrix(c(rep(states, each =2), rep(0, length(states)*2)), nrow=1)  
    TJ <- - JBH %*% solve(t(JV))
    predicted <- vCF1[,i-1,e] + as.numeric((steps1[i] - steps1[i-1]) * TJ)
    
    theta$betaH <- steps1[i]
    EQcf <- multiroot(function(V){PSI(V, theta, Trans, G) - V}, 
                      predicted, maxiter = 400, rtol=1e-10, atol=1e-12, ctol=1e-12,
                      jacfunc = function(V){PsiDer(V, theta, Trans, G)-diag(length(V))},
                                            jactype="fullusr")
    vCF1[,i,e] <- EQcf$root
    
  }
}


## par(mfrow=c(2,3))
## plot(vCF0[50,,1])
## plot(vCF0[50,,2])
## plot(vCF0[50,,3])
## plot(vCF1[50,,1])
## plot(vCF1[50,,2])
## plot(vCF1[50,,3])

## Above plots say we can move only 9 steps in the postive direction
## and up to 40 steps in the negative direction. We picked 30 steps.
use0 <- 9
use1 <- 30

ggCF_Hbeta <-  data.frame()

for (e in 1:3){

    CP_eq <-  AttackProbs(EQall[,e])
    CP_cf0 <- AttackProbs(vCF0[,use0,e])
    CP_cf1 <- AttackProbs(vCF1[,use1,e])
    CP_eq$prAny <- CP_eq$prAH +   CP_eq$prAF -  CP_eq$prAH* CP_eq$prAF
    CP_cf0$prAny <- CP_cf0$prAH +   CP_cf0$prAF -  CP_cf0$prAH* CP_cf0$prAF
    CP_cf1$prAny <- CP_cf1$prAH +   CP_cf1$prAF -  CP_cf1$prAH* CP_cf1$prAF
    
    ggCF_Hbeta_e <-  data.frame(prattack = c(CP_cf0$prAH, CP_eq$prAH, CP_cf1$prAH, 
                                             CP_cf0$prAF, CP_eq$prAF, CP_cf1$prAF,
                                             CP_cf0$prAny, CP_eq$prAny, CP_cf1$prAny),
                                states = rep(states, 9),
                                actor = rep(c("Hamas", "Fatah", "Either"), each=3*length(states)),
                                cfb = as.factor(rep(c(0, 1, 2, 0, 1, 2, 0, 1, 2), each= length(states))))  
    ggCF_Hbeta <- rbind(ggCF_Hbeta, ggCF_Hbeta_e)
    
}

ggCF_Hbeta$actor <- factor(ggCF_Hbeta$actor, level=c("Fatah", "Hamas", "Either"), ordered=T)
ggCF_Hbeta$eq <- factor(rep(c("Symmetric Eq.","Hamas-Dominant Eq.", "Fatah-Dominant Eq."), each=dim(ggCF_Hbeta)[1]/3),
                        level =  c("Fatah-Dominant Eq.", "Symmetric Eq.", "Hamas-Dominant Eq."), ordered=T)

pe <- ggplot(ggCF_Hbeta, aes(x=states, y=prattack, color=cfb)) +
  geom_path(linewidth=1.2) + theme_bw(11) + theme(legend.position="bottom") +
  theme(legend.margin=margin(-6,0,0,0), legend.box.margin=margin(-6,0,0,0))+
  xlab("Relative Popularity (states)") + ylab("Pr. Attack") + facet_grid(rows = vars(actor), cols=vars(eq)) +
  scale_color_manual(values=c("#fd8d3c", "#e6550d", "#a63603" ),
                     name=expression(beta[H]),
                     breaks=c("0", "1", "2"),
                     labels=
                       as.character(
                         round(c(steps0[use0], thetaStar$betaH, steps1[use1]), digits=4)) )
  

ggsave("../../Output/Figures/figureB4.pdf", pe, height = 5, width= 6) 



####### single agent comparison #######
###############################################
## Hamas 
gammaStar <- c(0, -1, 0, 1, 0, 0)
sigmaStar <- 2


Trans <- gamma2trans(gammaStar, sigmaStar, states, d=1, discretize=F)
uas <- cbind(thetaStar$betaH*G$state + (thetaStar$kappaH[1] + thetaStar$kappaH[2]*G$state)*G$aH,
             thetaStar$betaF*G$state + (thetaStar$kappaF[1] + thetaStar$kappaF[2]*G$state)*G$aF)
colnames(uas) <- c("uH", "uF")
vstart <- c(rowMeans(matrix(uas[,'uH'], ncol=2, byrow=T)),
            c(apply(array(uas[,'uF'], c(2,2,length(states))), 3, rowMeans)))

thetaStar$delta <- 0.999
EQHamas <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, vstart, 
                   maxiter = 200, 
                   jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))},
                   jactype="fullusr")

EQCPhamas <- AttackProbs(EQHamas$root)
ggdatasym <-  data.frame(PrAttack =c(EQCPhamas$prAH, EQCPhamas$prAF),
                         states = rep(states,2),
                         actor = rep(c("Hamas", "Fatah"), each = length(states)))





## Fatah
gammaStar[2] <- 0 
gammaStar[3] <- 1
Trans <- gamma2trans(gammaStar, sigmaStar, states, d=1, discretize=F)
vHamas <- matrix(EQHamas$root[1:(2*length(states))],  ncol=length(states)) 
vFatah <- matrix(EQHamas$root[(2*length(states)+1):(4*length(states))],  ncol=length(states))
start <- c(vFatah[,length(states):1], vHamas[,length(states):1])

EQFatah <- multiroot(function(V){PSI(V, thetaStar, Trans, G) - V}, start, 
                     maxiter = 200, 
                     jacfunc = function(V){PsiDer(V, thetaStar, Trans, G)-diag(length(V))},
                     jactype="fullusr")
EQCPfatah <- AttackProbs(EQFatah$root)

ggdata <- data.frame(PrAttack =c(EQCPhamas$prAH, EQCPfatah$prAF),
                     states = rep(states,2),
                     actor = rep(c("Hamas's Single Agent Problem", "Fatah's Single Agent Problem"), each = length(states)))

pSAP <- ggplot(ggdata, aes(x=states, y=PrAttack)) + geom_line(linewidth=1.25) + theme_bw(11) +
  xlab("Relative Popularity (states)") + ylab("Pr. Attack")  + 
  facet_grid(cols=vars(actor)) + ylim(c(0.095, 0.71))

ggsave("../../Output/Figures/figureB1.pdf" , height = 2.5, plot = pSAP, width= 6)


### correlation in text 
cor(ggdata$PrAttack[ggdata$actor=="Hamas's Single Agent Problem"], ggdata$PrAttack[ggdata$actor=="Fatah's Single Agent Problem"])

