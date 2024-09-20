#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Single agent counterfactual
#' ---
rm(list=ls())

######################################
# packages and functions
library("matrixStats")
library("ggplot2")
library("rootSolve")
library("Matrix")

source("helperFunctions.r")
source("gamma2trans.R")
# ######################################

###############################################
# results
load("Results/firstStageOutput.Rdata")
gammaStar <- coefficients(mod0)
sigmaStar <- summary(mod0)$sigma
rm(states, Trans)


load("Results/mainModel.rdata")
states <- model.main[["states"]]
Trans <- model.main[["Trans"]]
thetaEst <- model.main[["regtable"]]$V1
delta <-  as.numeric(model.main[["params"]]$delta)
vEst <- model.main[["v"]]

mainData <- model.main$regData
mainData$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
mainData$states.int <- sapply(mainData$states.discrete, function(x){which(x==states)})

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

EQCP <- AttackProbs(vEst)
EQIV <- invarDist(vEst, Trans)
vH <- log(colSums(exp(matrix(vEst[1:(2*length(states))],  ncol=length(states)))))
vF <- log(colSums(exp(matrix(vEst[(2*length(states)+1):(4*length(states))],  ncol=length(states)))))




###############################################
# Hamas single agent problem

# new transition matrix
gammaHSA <- gammaStar
gammaHSA['lag.Fattacks'] <- 0
gammaHSA['lag.Fattacks:lag.states'] <- 0
TransHSA <- gamma2trans(gammaHSA, sigmaStar, states, discretize=F,d = .05, bound=0.025)

start <- vEst
EQHSA <- multiroot(function(V){PSI(V, thetaStar, TransHSA, G) - V}, start,
                   jactype="fullusr",
              maxiter = 200, rtol=1e-8, atol=1e-8,
              jacfunc = function(V){PsiDer(V, thetaStar, TransHSA, G)-diag(length(V))})
EQCPHSA <- AttackProbs(EQHSA$root)

###############################################
# Fatah single agent problem

# new transition matrix
gammaFSA <- gammaStar
gammaFSA['lag.Hattacks'] <- 0
gammaFSA['lag.Hattacks:lag.states'] <- 0
TransFSA <- gamma2trans(gammaFSA, sigmaStar, states,  discretize=F,d = .05, bound=0.025)

start <- vEst
EQFSA <- multiroot(function(V){PSI(V, thetaStar, TransFSA, G) - V}, start,
                   maxiter = 200, rtol=1e-8, atol=1e-8,jactype="fullusr",
                   jacfunc = function(V){PsiDer(V, thetaStar, TransFSA, G)-diag(length(V))})
EQCPFSA <- AttackProbs(EQFSA$root)
###############################################
diffHamas <- EQCP$prAH - EQCPHSA$prAH
diffFatah <- EQCP$prAF - EQCPFSA$prAF

ggdif2time <- data.frame(change = c(diffHamas[mainData$states.int],
                                    diffFatah[mainData$states.int]),
                         time = as.Date(rep(mainData$Date,2), "%Y-%m-%d"),
                         actor = rep(c("Hamas", "Fatah"),each=dim(mainData)[1]))


pd2time <- ggplot(ggdif2time) +
  geom_line(aes(y=change, x=time, color = as.factor(actor), linetype = as.factor(actor)), size=1.15) +
  theme_bw(11) +
  xlab("Time") + ylab("Change in Pr. Attack") +
  scale_color_manual(values=c("navyblue", "orangered"),
                     name="Actor") + scale_linetype_manual(values= c(1, 2), name="Actor") +
  geom_hline(yintercept=0,size=0.9) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  theme(legend.position=c(0.125, 0.8)) +
  theme(legend.text=element_text(size=rel(0.8)),
        legend.title=element_text(size=rel(0.9)),
        legend.key.width = unit(.5, 'inches')) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="grey50")) +
  geom_vline(aes(xintercept=as.Date("2000-09-28")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2000-09-28"), y = 0.21, label = "2nd Intifada",
           hjust=0, size=3, angle=-90, vjust=.1) +
  geom_vline(aes(xintercept=as.Date("2006-1-1")), linetype="dotted", alpha=.3)+
  annotate("text", x = as.Date("2006-1-30"), y =.3 , label = "Hamas wins 2006\nlegislative elections",
           hjust=0, size=3, angle=-90, vjust=.5) +
  scale_x_date(breaks=as.Date(paste(seq(from=1994,to=2018, by=4), "-02-01",sep="")),
               labels = seq(from=1994,to=2018, by=4)) +
  theme(panel.grid.major.x = element_blank(),
        strip.text=element_text(size=10))
ggsave("../../Output/Figures/figure5.pdf", pd2time, width=7, height=4)



# compute average treatment effect before 2nd Intifada
T1 <- ggdif2time[ggdif2time$time < as.Date("2000-09-01"),]
T2 <- ggdif2time[ggdif2time$time >= as.Date("2000-09-01") & ggdif2time$time < as.Date("2006-01-01"),]
T3 <- ggdif2time[ggdif2time$time >= as.Date("2006-01-01"),]


sATE <- data.frame(
  TimePeriod1 = c(mean(T1[T1$actor=="Hamas","change"]),
                  sd(T1[T1$actor=="Hamas","change"])/sqrt(dim(T1)[1]),
                  mean(T1[T1$actor=="Fatah","change"]),
                  sd(T1[T1$actor=="Fatah","change"])/sqrt(dim(T1)[1])),
  TimePeriod2 = c(mean(T2[T2$actor=="Hamas","change"]),
                  sd(T2[T2$actor=="Hamas","change"])/sqrt(dim(T2)[1]),
                  mean(T2[T2$actor=="Fatah","change"]),
                  sd(T2[T2$actor=="Fatah","change"])/sqrt(dim(T2)[1])),
   TimePeriod3 = c(mean(T3[T3$actor=="Hamas","change"]),
                  sd(T3[T3$actor=="Hamas","change"])/sqrt(dim(T3)[1]),
                  mean(T3[T3$actor=="Fatah","change"]),
                  sd(T3[T3$actor=="Fatah","change"])/sqrt(dim(T3)[1])))

row.names(sATE) <- c("Hamas", "", "Fatah", " ")
sATE.out <- num2str(as.matrix(sATE))
sATE.out[3,1] <- paste(round(sATE[3,1], 3)) #take it out 1 more decimal
sATE <- sATE.out
sATE[c(2,4),] <- paste0("(", sATE[c(2,4),], ")")
colnames(sATE) <- c("Oslo era",
                    "2nd Intifada}",
                    "post 2006 election")


cat(kable(sATE, digits=2, format="pipe"),
    file="../../Output/Tables/table4.txt",
    sep="\n")


#######################################
## substantive effects mentioned in text
## Fatah uses 34% more violence b/c of competition
round(mean((EQCP$prAF[mainData$states.int] -
            EQCPFSA$prAF[mainData$states.int])/EQCPFSA$prAF[mainData$states.int]),3)*100
## Hamas uses 37% more violence b/c of competition
round(mean((EQCP$prAH[mainData$states.int] -
            EQCPHSA$prAH[mainData$states.int])/ EQCPHSA$prAH[mainData$states.int]), 3)*100




# During the Oslo lull, Hamas would use 5% more violence absent competition from Fatah
round(mean(EQCPHSA$prAH[mainData$states.int[1:81]] / EQCP$prAH[mainData$states.int[1:81]])   -1, 3)*100



# Hamas  used 4-5% less violence during Oslo lull than if it thought Fatah would never attack.
mean((EQCP$prAH[mainData$states.int][1:81] -EQCPHSA$prAH[mainData$states.int][1:81])/ EQCPHSA$prAH[mainData$states.int][1:81])



### largest effect size with more violence in counterfactual
m <-  max(EQCPHSA$prAH[mainData$states.int[1:81]]/EQCP$prAH[mainData$states.int[1:81]]) 
round(100*(m-1),1) # max percentage change about 9
round(sum(mainData$Hattacks[1:81]) * (m-1), 2) #change in # of months w/ Hammas terrorism




ggdif2 <- data.frame(change = c(EQCP$prAH - EQCPHSA$prAH,
                                   EQCP$prAF - EQCPFSA$prAF),
                     states = rep(states, 2),
                     actor = rep(c("Hamas", "Fatah"),each=length(states)))

pd2states <- ggplot(ggdif2) +
  geom_line(aes(y=change, x=states, color = as.factor(actor), linetype = as.factor(actor)), size=1.15) +
  theme_bw(16) +
  xlab("Relative Popularity (states)") + ylab("Change in Pr. Attack") +
  scale_color_manual(values=c("navyblue", "orangered"),
                     name="Actor") + scale_linetype_manual(values= c(1, 2), name="Actor") +
  geom_rug(data = mainData, aes(x = states.discrete), inherit.aes=F, color="grey50") +
  geom_hline(yintercept=0,size=0.9) +
  theme(legend.position=c(0.85, 0.8)) +
  theme(legend.text=element_text(size=rel(0.6)),
        legend.title=element_text(size=rel(0.8)),
        legend.key.width = unit(.5, 'inches')) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="grey50"))
ggsave("../../Output/Figures/figureA2.pdf", pd2states, width=7, height=4)
