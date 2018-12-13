# Step 7. Using fixed parameters and kernel, re-train Coef (alphas) on a subset training data.
# Compute the new errors of prediction: RMSE_train and RMSE_test, and new performances.
# Repeat 20 times to have 20 diff Perf/RMSE for train/test, so 20 diff Coef (alphas), and 20 perfs.
# Get mean and std dev of the 20 Perf_train/test. Conclusions on prediction reliability.

gamma1 <- 14788
lambda1 <- unname(eigen(expFunc(gamma1))$values[444]) 

set.seed(1)
index.train <- sample(1:nrow(x.ceo), 0.85*nrow(x.ceo))
index.test <- c(1:nrow(x.ceo))[-index.train]  # Test set

# Function for predictions and perf
'============================='
subset=index.test
'============================='
mtdist <- function(data){
  res <- matrix(0,nrow(data),nrow(data))
  for (i in 1:nrow(data)){
    for (j in 1:nrow(data)){
      res[i,j] <- crossprod((data[i,] - data[j,]))
    }
    
  }
  return(res)
}
#
expFunc.FIX <- function(gam=gamma1){
  return(exp(-gam * mtdist(x.ceo[subset,])^2))
}
#
alphaFunc.FIX <- function(G, lam){
  return(solve(G + lam * diag(1, nrow(G)), y.ceo[subset]))
}
#
predFunc.FIX <- function(G, lam){
  pred <- alphaFunc.FIX(G, lam) %*% G -mean(alphaFunc.FIX(G, lam) %*% G)
  RMSE2 <- mean(abs(y.ceo[subset]-pred)^2)
  c <- mean(abs(ceo[, 1]))
  perf <- sqrt(RMSE2)/c
  #return(c(RMSE=round(sqrt(RMSE2),4), Perf=round(perf,4)))
  return(perf)
}

predFunc.FIX(G=expFunc.FIX(14788), lam=lambda1)

perftest.res <- c()
for (i in c(1:20)){
  set.seed(i)
  index.train <- sample(1:nrow(x.ceo), 0.85*nrow(x.ceo))
  index.test <- c(1:nrow(x.ceo))[-index.train] 
  
  subset = index.train
  perftrain = predFunc.FIX(G=expFunc.FIX(14788), lam=lambda1)
  
  subset = index.test
  perftest = predFunc.FIX(G=expFunc.FIX(14788), lam=lambda1)
  perftest.res <- rbind(perftest.res, c(Perf.Train=perftrain, Perf.Test=perftest))
}

perftest.res <- data.frame(perftest.res)
perftest.res

avg.perftest <- apply(perftest.res, 2, mean)
avg.perftest

std.perftest <- apply(perftest.res, 2, sd)
std.perftest

write.csv(perftest.res, file='perftest_result.csv')



eigen.exp14800 <- eigen(expFunc(14790))$values
lam.14800 <- cutoff(eigen.exp14800); lam.14800
which(eigen.exp14800==lam.14800)
predFunc(expFunc(14787), lam.14800)

cutoff.res <- eigen(expFunc(14800))$values[444]/eigen(expFunc(14800))$values[1]; cutoff.res  # >5% cutoff
exp.TryFunc(14800, 9.99, plot=T)

exp.TryFunc(15000, 9.99, plot=T)








