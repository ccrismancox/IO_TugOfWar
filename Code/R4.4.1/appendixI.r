#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Appendix I: Choice of discount factors
#' ---

library(data.table)
library(stringr)
library(knitr)
library(ggplot2)
rm(list=ls())
load("Results/firstStageOutput.Rdata")
load("Results/mainModel.rdata")
load("Results/mainStateActions.Rdata")
load("Results/ChangingDeltas.rdata")
load("Results/startingvalues.rdata")
source("helperFunctions.r")
set.seed(1)
Delta <- c(0,.9,.925,.95,.975, .99, .999,.9999)
for(d in Delta){
  cat("Now on Delta = ", d, "\n")
  py.params <- data.table(delta=d,
                          nkappa=2, nbeta=2)
  regData <- mainData[,list(states, Hattacks, Fattacks, lag.states, states.discrete)]
  mod <- try(eval(as.name(paste("model.main.d", d,sep=""))), silent = TRUE)
  if("try-error" %in% class(mod)){
    mod <- try(eval(as.name("startmod")))
    start <- c(startmod$regtable$V1,startmod$v)
  }else{
    start <- c(mod$regtable$V1,mod$v$V1)

  }
  system("mkdir -p ipoptTEMP_delta")
  write.csv(regData, file=paste("ipoptTEMP_delta/regData.csv", sep=""),row.names = F)
  write.csv(Trans, file=paste("ipoptTEMP_delta/trans.csv", sep=""),row.names = F)
  write.csv(states, file=paste("ipoptTEMP_delta/statespace.csv", sep=""),row.names = F)
  write.csv(py.params, file=paste("ipoptTEMP_delta/params.csv", sep=""),row.names = F)
  write.csv(start, file=paste("ipoptTEMP_delta/start.csv", sep=""),row.names = F)
  write.csv(c(mod0$coef, summary(mod0)$sigma), file=paste("ipoptTEMP_delta/gamma.csv", sep=""),row.names = F)
  write.csv(bootOut, file=paste("ipoptTEMP_delta/V1.csv", sep=""),row.names = F)

  cat("IPOPT running, please wait"\n)  
  system(paste("python ../Python3/fitChangingDeltas.py"), ignore.stdout = TRUE)
  cat("IPOPT complete\n")
  
  regtable <- read.csv(paste("ipoptTEMP_delta/regtable.csv", sep=""),header = F)
  conv <- read.csv(paste("ipoptTEMP_delta/convergence.csv", sep=""),header = F)
  
  v <- read.csv(paste("ipoptTEMP_delta/v.csv", sep=""),header = F)
  assign(
    paste("model.main.d", d, sep=""),
    list(regData=regData, Trans=Trans, states=states, params=py.params,
         conv=conv, v=v,  regtable=regtable)
  )
  save(list=ls()[str_detect(ls(), "model.main")], file="Results/ChangingDeltas.rdata")
  system("rm ipoptTEMP_delta/ -r")
  
}

rm(list=ls())
load("Results/ChangingDeltas.rdata")
source("helperFunctions.r")
Delta <- c(0,.9,.925,.95,.975, .99, .999,.9999)


#### table #### 
X <- lapply(Delta, function(d){
  mod <- eval(as.name(paste("model.main.d", d,sep="")))
  out <- c(matrix(t(mod$regtable[, c("V1", "V2")]), ncol=1, byrow=T), unlist(-1*mod$conv))
  out <- round(out, 2)
  return(out)
})
X <- do.call(cbind.data.frame,X)

X[c(2,4,6,8),] <-paste("(",unlist(X[c(2,4,6,8),]),")", sep="")
X[10,] <- ifelse(X[9,]=="0", X[10,], paste(X[10,],"^*",sep=""))

X[1:4,1] <- " " # Betas not Identified when delta =0

names(X) <- Delta
X  <- X[-9,] # remove convergence codes
delta <- c('beta_H', "",
           'beta_F', "",
           'kappa_H', "",
           'kappa_F', "",
           'LL')
X <- cbind(delta, X)

cat(kable(X, caption="Discount factors and model fit", format="pipe"),
    file="../../Output/Tables/tableI1.txt", sep="\n")
cat("Note: ^* Model failed to converge.", file="../../Output/Tables/tableI1.txt", 
    append=TRUE)



### Figure #### 

vEstd0.9 <- c(unlist(model.main.d0.9[["v"]]))
vEstd0.925 <- c(unlist(model.main.d0.925[["v"]]))
vEstd0.95 <- c(unlist(model.main.d0.95[["v"]]))
vEstd0.975 <- c(unlist(model.main.d0.975[["v"]]))
vEstd0.99 <- c(unlist(model.main.d0.99[["v"]]))
vEstd0.999 <- c(unlist(model.main.d0.999[["v"]]))
vEstd0.9999 <- c(unlist(model.main.d0.9999[["v"]]))
states <- model.main.d0.999[["states"]]

###############################################
# create states
G <- expand.grid(states, c(0,1), c(0,1))
names(G) = c("state", "aH", "aF")
G <- G[order(G$state, G$aH),]
rownames(G) <- c()
nG <- dim(G)[1]

###############################################
# Choice probabilites
EQCP9999 <- AttackProbs(vEstd0.9999)
EQCP999 <- AttackProbs(vEstd0.999)
EQCP99 <- AttackProbs(vEstd0.99)
EQCP975 <- AttackProbs(vEstd0.975)
EQCP95 <- AttackProbs(vEstd0.95)
EQCP925 <- AttackProbs(vEstd0.925)
EQCP9 <- AttackProbs(vEstd0.9)


# distances
DHab  <- abs(EQCP999$prAH - EQCP99$prAH)
DFab  <- abs(EQCP999$prAF - EQCP99$prAF)

cat("Average difference in choice probs for Hamas", round(100*mean(DHab),1), "percentage points\n")
cat("Average difference in choice probs for Fatah", round(100*mean(DFab),1), "percentage points\n")
cat("Max difference in choice probs for Hamas", round(100*max(DHab),1), "percentage points\n")
cat("Max difference in choice probs for Fatah", round(100*max(DFab),1), "percentage points\n")


# figure
ggdata <-  data.frame(PrAttack = c(EQCP9999$prAH, EQCP9999$prAF,
                                   EQCP999$prAH, EQCP999$prAF,
                                   EQCP99$prAH, EQCP99$prAF,
                                   EQCP975$prAH, EQCP975$prAF,
                                   EQCP95$prAH, EQCP95$prAF,
                                   EQCP925$prAH, EQCP925$prAF,
                                   EQCP9$prAH, EQCP9$prAF),
                      states = rep(states,2*7),
                      delta = as.factor(rep(c(0.9999, 0.999, 0.99, 
                                              0.975, 0.95, 0.925, 
                                              0.9), 
                                        each =2*length(states))),
                      actor = rep(rep(c("Hamas", "Fatah"), each = length(states)),7))


pccp_delta <- ggplot(ggdata, aes(x=states, y=PrAttack, color=actor, linetype=actor)) + 
  geom_line(linewidth=1.5) + theme_bw(18) +
  facet_wrap(vars(delta), nrow=2) +
  xlab("Relative Popularity (States)") + ylab("Pr. Attack") + 
  scale_color_manual(values=c("navyblue", "orangered"),
                     name="Actor") + 
  scale_linetype_manual(values= c(1, 2), name="Actor") + 
  theme(legend.position="bottom",
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18),
        legend.key.width = unit(.75,"in"))

ggsave("../../Output/Figures/figureI1.pdf" , width = 11, height= 8.5)  
