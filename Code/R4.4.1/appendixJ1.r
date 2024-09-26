#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Changes to the discretization part 1
#' ---

library(data.table)
library(stringr)
library(knitr)
rm(list=ls())
load("Results/mainStateActions.Rdata")
load("Results/firstStageOutput.Rdata")
load("Results/MARSSoutput.rdata")
load("Results/changingd_B.rdata")
source("gamma2trans.R")
source("firststageboot.r")

D <- c(0.025, .05, 0.075)
BOUND <- c(0.0125, 0.025, 0.0375)
grid <- expand.grid(D, BOUND)

# grid <- cbind(D, BOUND)
for(i in 1:nrow(grid)){
  set.seed(1)

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
                          twostep=FALSE,
                          maxit=1000)
  regData <- mainData[,.(states, Hattacks, Fattacks, lag.states, states.discrete)]

  system("mkdir -p ipoptTEMP")
  write.csv(regData, file="ipoptTEMP/regData.csv",row.names = F)
  write.csv(Trans, file="ipoptTEMP/trans.csv",row.names = F)
  write.csv(states, file="ipoptTEMP/statespace.csv",row.names = F)
  write.csv(py.params, file="ipoptTEMP/params.csv",row.names = F)
  write.csv(start, file="ipoptTEMP/start.csv",row.names = F)
  write.csv(c(mod0$coef, summary(mod0)$sigma), file="ipoptTEMP/gamma.csv",row.names = F)
  write.csv(bootOut, file="ipoptTEMP/V1.csv",row.names = F)

  

  cat("IPOPT running, please wait\n")
  system(paste("python ../Python3/fitMainModel.py"), ignore.stdout = TRUE)
  cat("IPOPT complete\n")

  regtable <- read.csv("ipoptTEMP/regtable.csv",header = F)
  conv <- read.csv("ipoptTEMP/convergence.csv",header = F)
  v <- read.csv("ipoptTEMP/v.csv",header = F)[,1]

  assign(
    paste("model.main.d", d, "_bound", bound, sep=""),
    list(regData=regData, Trans=Trans, states=states, params=py.params,
         conv=conv, v=v,  regtable=regtable)
  )
  cat("Done with  d = ", d, " and bound = ", bound, "\n")
  save(list=ls()[str_detect(ls(), "model.main")], file="Results/changingd_B.rdata")
  system("rm ipoptTEMP/ -r")
}

rm(list=ls())
load("Results/changingd_B.rdata")
names <- ls()
models.d <- lapply(names,\(x){eval(as.name(x))})
names(models.d) <- names
source("helperFunctions.r")

D <- c(0.025, .05, 0.075)
BOUND <- c(0.0125, 0.025, 0.0375)
grid <- expand.grid(D, BOUND)
grid <- grid[order(grid$Var1),] #match the order in the model list

X <- lapply(models.d, function(mod){
  out <- c(matrix(t(mod$regtable[, c("V1", "V2")]), ncol=1, byrow=T), unlist(-1*mod$conv))
  temp <- c(num2str(out[1:8]), round(out[9], 1), num2str(out[10]))
  return(temp)
})
X <- do.call(cbind.data.frame,X)
K <- sapply(models.d, \(x){length(x$states)})
X <- as.data.frame(rbind(as.matrix(X),
                         t(as.matrix(grid)),
                         K))

X[c(2,4,6,8),] <- paste0("(",unlist(X[c(2,4,6,8),]),")")
X[10,] <- ifelse(X[9,]=="0", X[10,], paste(X[10,],"^*",sep=""))

X <- X[-9,]
X[nrow(X),] <- paste0(X[ nrow(X),])

` ` <- c('beta_H', "",
         'beta_F', "",
         'kappa_H', "",
         'kappa_F', "",
         "LL", "2d",
         "End percentile", "K")
X <- cbind(` `, X)
colnames(X) <- NULL
X <- rbind(as.matrix(X[10:12,]), as.matrix(X[1:9,]))
rownames(X) <- NULL
cat(kable(X, caption="Estimates at different discretizations of tilde{s}",
          format="pipe"),
    file="../../Output/Tables/tableJ1.txt", sep="\n")
