#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title:  Appendix J part 2
#' ---
rm(list=ls())
library(data.table)
library(stringr)
library(knitr)
library(ggplot2)

load("Results/mainStateActions.Rdata")
load("Results/firstStageOutput.Rdata")
load("Results/MARSSoutput.rdata")
load("Results/CoarseResults.rdata")
source("gamma2trans.R")
source("firststageboot.r")

D <- c(.75 , 1, 1.5)
BOUND <- c(0.025)
grid <- expand.grid(D, BOUND)

# grid <- cbind(D, BOUND)
for(i in 1:nrow(grid)){
  set.seed(5)
  
  d <- grid[i,1]
  bound <- grid[i,2]
  cat("Now on d = ", d, " and bound = ", bound, "\n")
  bounds <- quantile( mars$states[1,], c(bound, 1-bound))
  states <- seq(from=bounds[1], to=bounds[2], by=d)
  #### set up for second stage ####
  Trans <- gamma2trans(mod0$coef, summary(mod0)$sigma,
                       d=d, bound=bound,
                       mars.states = regData$states)
  mainData$states.discrete <- sapply(mainData$states, 
                                     function(x){return(states[which.min((states-x)^2)])})
  
   mod <- try(eval(as.name(paste("model.main.d", d, "_bound", bound, sep=""))), silent = TRUE)
  if(!"try-error" %in% class(mod)){
    start <- c(mod$regtable$V1,mod$v)
  }else{
    start <- runif(4+length(states)*4)
    
  }  
  py.params <- data.table(delta=.999, 
                          nkappa=2, nbeta=2,
                          twostep=FALSE)
  regData <- mainData[,.(states, Hattacks, Fattacks, lag.states, states.discrete)]
  
  system("mkdir -p ipoptTEMP")
  write.csv(regData, file="ipoptTEMP/regData.csv",row.names = F)
  write.csv(Trans, file="ipoptTEMP/trans.csv",row.names = F)
  write.csv(states, file="ipoptTEMP/statespace.csv",row.names = F)
  write.csv(py.params, file="ipoptTEMP/params.csv",row.names = F)
  write.csv(start, file="ipoptTEMP/start.csv",row.names = F)
  write.csv(c(mod0$coef, summary(mod0)$sigma), file="ipoptTEMP/gamma.csv",row.names = F)
  write.csv(bootOut, file="ipoptTEMP/V1.csv",row.names = F)
  
  system("python ../Python3/fitMainModel.py", ignore.stdout = FALSE)
  
  regtable <- read.csv("ipoptTEMP/regtable.csv",header = F)
  conv <- read.csv("ipoptTEMP/convergence.csv",header = F)
  v <- read.csv("ipoptTEMP/v.csv",header = F)[,1]
  
  assign(
    paste("model.main.d", d, "_bound", bound, sep=""),
    list(regData=regData, Trans=Trans, states=states, params=py.params,
         conv=conv, v=v,  regtable=regtable)
  )
  cat("Done with  d = ", d, " and bound = ", bound, "\n")
  save(list=ls()[str_detect(ls(), "model.main")], file="Results/CoarseResults.rdata")
  system("rm ipoptTEMP/ -r")
}




rm(list=ls())
load("Results/CoarseResults.rdata")
names <- ls()
models.d <- lapply(names,\(x){eval(as.name(x))})
names(models.d) <- names

source("helperFunctions.r")


D <- c(.75 , 1, 1.5)
BOUND <- c(0.025)
grid <- expand.grid(D, BOUND)
grid <- grid[order(grid$Var1),] #match the order in the model list

X <- lapply(models.d, function(mod){
  outTab <- c(matrix(t(mod$regtable[, c("V1", "V2")]), ncol=1, byrow=T), unlist(-1*mod$conv))
  CONV <- outTab[9]
  temp <- num2str(outTab[-9], digits=2)
  temp <- c(temp, CONV)
  return(temp)
})
X <- do.call(cbind.data.frame,X)
K <- sapply(models.d, \(x){length(x$states)})
X <- as.data.frame(rbind(as.matrix(X),
                         t(as.matrix(grid)),
                         K))

X[c(2,4,6,8),] <-paste("(",unlist(X[c(2,4,6,8),]),")", sep="")
X[9,] <- ifelse(X[10,]=="0", X[9,], paste(X[9,],"^*",sep=""))

X <- X[-10,]

` ` <- c('beta_H', "",
                 'beta_F', "",
                 'kappa_H', "",
                 'kappa_F', "",
                 "LL", "2d", "End percentile", "K")
X <- cbind(` `, X)
colnames(X) <- NULL
X <- rbind(as.matrix(X[10:12,]), as.matrix(X[1:9,]))
X["K",] <- paste0(X["K",])




## Pull in the main man
load("Results/mainModel.rdata")
Y <- num2str(t(model.main$regtable[,1:2]))
Y[2,] <- paste0("(", Y[2,], ")")
Y <- c(Y)
Y <- c(0.05, 0.025, "440", Y, num2str(-1*model.main$conv$V1[2]))
tabOut <- cbind(X[,1], Y, X[,-1])
colnames(tabOut) <- rownames(tabOut) <- NULL
cat(kable(tabOut,
          caption="Estimates at different discretizations of tilde{s}"),
    file="../../Output/Tables/tableJ2.txt", sep="\n")
            



###### Figures ######
##  d = 0.75
vEst <- model.main.d0.75_bound0.025$v
states <- model.main.d0.75_bound0.025$states; K <- length(states)
Obs075 <- model.main.d0.75_bound0.025$regData
Obs075$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
Obs075$states.int <- sapply(Obs075$states.discrete, function(x){which(x==states)})
Tperiod <- dim(Obs075)[1]
Obs075 <- subset(Obs075, select = c("Hattacks", "Fattacks", "states.discrete", "states.int", "Date"))
Obs075$grid <- paste("K = ", K, sep="")



EQCP <- AttackProbs(vEst)
ggdata075 <-  data.frame(PrAttack = c(EQCP$prAH, EQCP$prAF),
                      states = rep(states,2),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)),
                      grid = paste("K = ", K, sep=""))
gtdata075 <-  data.frame(PrAttack = c(EQCP$prAH[Obs075$states.int], EQCP$prAF[Obs075$states.int]),
                      time = as.Date(rep(Obs075$Date,2), "%Y-%m-%d"), 
                      actor = rep(c("Hamas", "Fatah"), each = Tperiod),
                      grid = paste("K = ", K, sep=""))


#######################################
#  d = 1
vEst <- model.main.d1_bound0.025$v
states <- model.main.d1_bound0.025$states; K <- length(states)
Obs1 <- model.main.d1_bound0.025$regData
Obs1$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
Obs1$states.int <- sapply(Obs1$states.discrete, function(x){which(x==states)})
Tperiod <- dim(Obs1)[1]
Obs1 <- subset(Obs1, select = c("Hattacks", "Fattacks", "states.discrete", "states.int", "Date"))
Obs1$grid <- paste("K = ", K, sep="")

EQCP <- AttackProbs(vEst)
ggdata1 <-  data.frame(PrAttack = c(EQCP$prAH, EQCP$prAF),
                      states = rep(states,2),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)),
                      grid = paste("K = ", K, sep=""))
gtdata1 <-  data.frame(PrAttack = c(EQCP$prAH[Obs1$states.int], EQCP$prAF[Obs1$states.int]),
                        time = as.Date(rep(Obs1$Date,2), "%Y-%m-%d"), 
                        actor = rep(c("Hamas", "Fatah"), each = Tperiod),
                        grid = paste("K = ", K, sep=""))


#######################################
#  d = 1.5
vEst <- model.main.d1.5_bound0.025$v
states <- model.main.d1.5_bound0.025$states; K <- length(states)
Obs15 <- model.main.d1.5_bound0.025$regData
Obs15$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
Obs15$states.int <- sapply(Obs15$states.discrete, function(x){which(x==states)})
Tperiod <- dim(Obs15)[1]
Obs15 <- subset(Obs15, select = c("Hattacks", "Fattacks", "states.discrete", "states.int", "Date"))
Obs15$grid <- paste("K = ", K, sep="")

EQCP <- AttackProbs(vEst)
ggdata15 <-  data.frame(PrAttack = c(EQCP$prAH, EQCP$prAF),
                      states = rep(states,2),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)),
                      grid = paste("K = ", K, sep=""))
gtdata15 <-  data.frame(PrAttack = c(EQCP$prAH[Obs15$states.int], EQCP$prAF[Obs15$states.int]),
                         time = as.Date(rep(Obs15$Date,2), "%Y-%m-%d"), 
                         actor = rep(c("Hamas", "Fatah"), each = Tperiod),
                         grid = paste("K = ", K, sep=""))



#######################################
#  BASLINE
states <- model.main[["states"]]
vEst <- model.main[["v"]]
Obsbase <-  model.main[["regData"]]
Obsbase$Date <- seq(as.Date("1994/1/1"), by = "month", length.out = 300)
Obsbase$states.int <- sapply(Obsbase$states.discrete, function(x){which(x==states)})
Tperiod <- dim(Obsbase)[1]
Obsbase <- subset(Obsbase, select = c("Hattacks", "Fattacks", "states.discrete", "states.int", "Date"))
Obsbase$grid <-  "Baseline"

EQCP <- AttackProbs(vEst)
ggdatabase <-  data.frame(PrAttack = c(EQCP$prAH, EQCP$prAF),
                      states = rep(states,2),
                      actor = rep(c("Hamas", "Fatah"), each = length(states)),
                      grid = "Baseline")
gtdatabase <-  data.frame(PrAttack = c(EQCP$prAH[Obsbase$states.int], EQCP$prAF[Obsbase$states.int]),
                          time = as.Date(rep(Obsbase$Date,2), "%Y-%m-%d"), 
                          actor = rep(c("Hamas", "Fatah"), each = Tperiod),
                          grid = "Baseline")



#######################################
## Figure J1 ##
ggdata <- rbind(ggdata075, ggdata1, ggdata15, ggdatabase)
ggdata$grid <- factor(ggdata$grid, levels = c("Baseline", "K = 30", "K = 22", "K = 15"), ordered=T)
Obs <- rbind(Obs075, Obs1, Obs15, Obsbase) 
Obs$grid <- factor(Obs$grid, levels = c("Baseline", "K = 30", "K = 22", "K = 15"), ordered=T)

pccp <- ggplot(ggdata, aes(x=states, y=PrAttack, color=actor, linetype=actor)) +
  geom_line(linewidth=1.5) + facet_wrap(~ grid, nrow=2) + 
  theme_bw(11) + 
  xlab("Relative Popularity (states)") + ylab("Pr. Attack") +
  scale_color_manual(values=c("navyblue", "orangered"),
                     name="Actor") + scale_linetype_manual(values= c(1, 2), name="Actor") +
  geom_rug(data = subset(Obs, Hattacks==1), aes(x=states.discrete, y=.3), sides = 'b', position=position_jitter(width=.15, height=0), inherit.aes=F, color="orangered1") +
  geom_rug(data = subset(Obs, Fattacks==1), aes(x=states.discrete,y=.3), sides = 'b',  position=position_jitter(width=.15, height=0),inherit.aes=F, color="navyblue")+
  theme(legend.key.width = unit(.5,"in"), legend.position = "bottom") +
  ylim(.04,.68)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0))  

ggsave("../../Output/Figures/figureJ1.pdf", pccp, width=6, height=4)

## Figure J2 ##
ggdata <- rbind(gtdata075, gtdata1, gtdata15, gtdatabase)
ggdata$grid <- factor(ggdata$grid, levels = c("Baseline", "K = 30", "K = 22", "K = 15"), ordered=T)
Obs <- rbind(Obs075, Obs1, Obs15, Obsbase) 
Obs$grid <- factor(Obs$grid, levels = c("Baseline", "K = 30", "K = 22", "K = 15"), ordered=T)

pccp <- ggplot(ggdata, aes(x=time, y=PrAttack, color=actor, linetype=actor)) +
  geom_line(linewidth=1.5) + theme_bw(11) + facet_wrap(~ grid, nrow=2) + 
  xlab("Time") + ylab("Pr. Attack") + 
  scale_color_manual(values=c("navyblue", "orangered"),name="Actor") + 
  scale_linetype_manual(values= c(1, 2), name="Actor") +
  geom_rug(data = subset(Obs, Hattacks==1), aes(x = Date), inherit.aes=F, color="orangered1") +
  geom_rug(data = subset(Obs, Fattacks==1), aes(x = Date), inherit.aes=F, color="navyblue") +
  theme(legend.key.width = unit(.75,"in"), legend.position = "bottom") +
  scale_x_date(breaks=as.Date(paste(seq(from=1994,to=2018, by=4), "-02-01",sep="")),
               labels = seq(from=1994,to=2018, by=4)) +
  theme(legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

ggsave("../../Output/Figures/figureJ2.pdf", pccp, width=6, height=4)

