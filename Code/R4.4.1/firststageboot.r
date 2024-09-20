firststageboot <- function(mod, regData, B=1000, workers=detectCores()-1, withCor=FALSE){
  
  mod0 <- mod
  if(workers >0){
    cl <- makeCluster(workers)
    registerDoParallel(cl)
    
    set.seed(1)
    
    bootOut <- foreach(b= 1:B, .combine = 'rbind', 
                       #.export=c("mod0", "regData"),
                       .packages = c("data.table"), 
                       .errorhandling = "remove") %dorng%{
                         # k <- 20000+b
                         
                         state <- c(regData$states[1], rep(0, 299))
                         if(withCor){
                           u <- c(rnorm(1, mean=0, sd=summary(mod0)$sigma), rep(0,299))
                           rho.model <- lm(mod0$residuals[2:299]~mod0$residuals[1:298])
                           rho.hat <- rho.model$coef[2]
                           sigma.rho <- summary(rho.model)$sigma
                           
                           
                           for(t in 2:300){
                             u[t] <- u[t-1]*rho.hat + rnorm(1, mean=0, sd=sigma.rho)
                             state[t] <- predict(mod0, newdata=data.frame(lag.states=state[t-1],
                                                                          lag.Hattacks=regData$lag.Hattacks[t],
                                                                          lag.Fattacks=regData$lag.Fattacks[t])) + u[t]
                             
                           }
                         }else{
                           for(t in 2:300){
                             state[t] <- rnorm(1,mean=predict(mod0, 
                                                              newdata=data.frame(lag.states=state[t-1],
                                                                                 lag.Hattacks=regData$lag.Hattacks[t],
                                                                                 lag.Fattacks=regData$lag.Fattacks[t])),
                                               sd=summary(mod0)$sigma)
                           }
                         }
                         regData.boot <- data.table(states=state,
                                                    lag.Hattacks=regData$lag.Hattacks,
                                                    lag.Fattacks=regData$lag.Fattacks)
                         regData.boot[,lag.states := shift(states)]
                         mod0a <- lm(states~lag.Hattacks +lag.Fattacks + lag.states+
                                       lag.Hattacks:lag.states+lag.Fattacks:lag.states, data=regData.boot)
                         
                         c(mod0a$coef,summary(mod0a)$sigma)
                         
                         # setTxtProgressBar(pb, i)
                         # save(bootOut,file="../../results/parametricbootTestA.rdata")
                         
                       }
    stopCluster(cl)
  }else{
    bootOut <- matrix(0, nrow=B, ncol=7)
    for(b in 1:B){
      state <- c(regData$states[1], rep(0, 299))
      if(withCor){
        u <- c(rnorm(1, mean=0, sd=summary(mod0)$sigma), rep(0,299))
        rho.model <- lm(mod0$residuals[2:299]~mod0$residuals[1:298])
        rho.hat <- rho.model$coef[2]
        sigma.rho <- summary(rho.model)$sigma
        
        
        for(t in 2:300){
          u[t] <- u[t-1]*rho.hat + rnorm(1, mean=0, sd=sigma.rho)
          state[t] <- predict(mod0, newdata=data.frame(lag.states=state[t-1],
                                                       lag.Hattacks=regData$lag.Hattacks[t],
                                                       lag.Fattacks=regData$lag.Fattacks[t])) + u[t]
          
        }
      }else{
        for(t in 2:300){
          state[t] <- rnorm(1,mean=predict(mod0, newdata=data.frame(lag.states=state[t-1],
                                                                    lag.Hattacks=regData$lag.Hattacks[t],
                                                                    lag.Fattacks=regData$lag.Fattacks[t])),
                            sd=summary(mod0)$sigma)
        }
      }
      regData.boot <- data.table(states=state,
                                 lag.Hattacks=regData$lag.Hattacks,
                                 lag.Fattacks=regData$lag.Fattacks)
      regData.boot[,lag.states := shift(state)]
      mod0a <- lm(states~lag.Hattacks +lag.Fattacks + lag.states+
                    lag.Hattacks:lag.states+lag.Fattacks:lag.states, data=regData.boot)
      
      bootOut[b,] <-  c(mod0a$coef,summary(mod0a)$sigma)
    }
  }
  colnames(bootOut) <- c("gamma0", "gamma_Hattacks", "gamma_Fattacks",
                         "gamma_lag", "gamma_Hattacks:lag", "gamma_Fattacks:lag",
                         "sigma")
  
  cov(bootOut)
}
