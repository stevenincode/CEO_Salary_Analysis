# Final project for ML
setwd("C:/Users/steven/Dropbox/UH Study/AAA Fall 2018/2. Topic-Data Clustering & ML/Final Project")
ceo <- read.csv('ceo.csv', row.names = 1)
#head(ceo)

# Standardization function
stdFunc <- function(data){
  res <- c()
  for (i in 1:ncol(data)){
    col <- (data[,i]-mean(data[,i]))/sd(data[,i])
    res <- cbind(res, col)
    colnames(res)[i] <- colnames(data)[i]
  }
  return(res)
}

x.ceo <- stdFunc(ceo[, -c(1, 2)])
y.ceo <- ceo[, 1] - mean(ceo[, 1])

'=================================================='
gamma1 <- 0.2
'=================================================='

# 3 kernel methods
mtdist <- function(data){
  res <- matrix(0,nrow(data),nrow(data))
  for (i in 1:nrow(data)){
    for (j in 1:nrow(data)){
      res[i,j] <- crossprod((data[i,] - data[j,]))
    }
    
  }
  return(res)
}

expFunc <- function(gam=gamma1){
  return(exp(-gam * mtdist(x.ceo)^2))
}

# Gramian matrix for 3 mtds
G.poly2 <- (1 + crossprod(t(x.ceo)))^2
#G.poly2[1:7,1:7]
G.poly3 <- (1 + crossprod(t(x.ceo)))^3
#G.poly3[1:7,1:7]
G.exp1 <- expFunc()
#G.exp1[1:7,1:7]

eigen.poly2 <- eigen(G.poly2)$value
#plot(eigen.poly2[1:20], main='Lambda of Poly2')
#lines(eigen.poly2[1:20])

eigen.poly3 <- eigen(G.poly3)$values
#plot(eigen.poly3[1:20], main='Lambda of Poly3')
#lines(eigen.poly3[1:20])

eigen.exp1 <- eigen(G.exp1)$values
#plot(eigen.exp1[1:20], main='Lambda of exp1')
#lines(eigen.exp1[1:20])

# Choosing cutoff of lambda (eigenvalue)
'=================================================='
cutoff.perc <- 9.99

expStart=1
expResMax=30
'=================================================='

cutoff <- function(lam, cutoff.perc=9.99){
  #lamsum <- sum(lam)
  cut <- cutoff.perc/100
  for (i in 1:length(lam)){
    perc <- lam[i]/lam[1]
    if (perc < cut){
      return(lam[i-1])
    } else {
      next
    }
  }
}

cutoff.before <- function(lam, cutoff.perc=9.99){
  #lamsum <- sum(lam)
  cut <- cutoff.perc/100
  lamList <- c()
  limit <- 0
  for (i in 1:length(lam)){
    perc <- lam[i]/lam[1]
    lamList <- c(lamList, lam)
    if (limit>=expResMax){
      print('Exceed expResMax limit')
      return(lamList[expStart:i])
    } else {
      if (perc < cut){
        return(lamList[expStart:i-1])
      } else {
        limit <- limit + 1
        next
      }
    }
  }
}

lam2 <- cutoff(eigen.poly2); lam2
lam3 <- cutoff(eigen.poly3); lam3
lam1 <- cutoff(eigen.exp1); lam1

plot(eigen.poly2[14:28], main='Lambda of Poly2')
lines(eigen.poly2[14:28])
abline(h=lam2, col='red')

plot(eigen.poly3[1:20], main='Lambda of Poly3')
lines(eigen.poly3[1:20])
abline(h=lam3, col='red')

plot(eigen.exp1[1:20], main='Lambda of exp')
lines(eigen.exp1[1:20])
abline(h=lam1, col='red')

# Coef Alpha of f(X) = A1 K(X, X_1) + A2 K(X, X_2) + ... + AN K(X, X_N)
alphaFunc <- function(G, lam){
  return(solve(G + lam * diag(1, nrow(G)), y.ceo))
}

a2 <- alphaFunc(G.poly2, lam2)
a3 <- alphaFunc(G.poly3, lam3)
a1 <- alphaFunc(G.exp1, lam1)

# Get b coef of Z = f(X) + b
#fx.2 <- apply(a2*G.poly2, 1, sum)
fx.2 <- a2 %*% G.poly2
b2 <- -mean(fx.2); b2

z.poly2 <- fx.2 + b2
mean(z.poly2)
z.poly2[1:6]

errorlist <- abs(y.ceo-(alphaFunc(G.poly2, lam2) %*% G.poly2 -mean(alphaFunc(G.poly2, lam2) %*% G.poly2)))

# Root mean squared error

#RMSE2 <- mean(abs(y.ceo-z.poly2)^2); sqrt(RMSE2)
#c2 <- mean(abs(ceo[, 1]))
#Perf2 <- sqrt(RMSE2)/c2; Perf2


# Function for predictions and perf
predFunc <- function(G, lam){
  pred <- alphaFunc(G, lam) %*% G -mean(alphaFunc(G, lam) %*% G)
  RMSE2 <- mean(abs(y.ceo-pred)^2)
  c <- mean(abs(ceo[, 1]))
  perf <- sqrt(RMSE2)/c
  #return(c(RMSE=round(sqrt(RMSE2),4), Perf=round(perf,4)))
  return(perf)
}

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

# Trying lambda for poly
#predFunc(G.poly2, lam2)

poly.TryFunc <- function(degree=2, cut=9.99, pick=-1, plot=F){
  if (degree==2){
    G.polyTry <- G.poly2
  } else {
    G.polyTry <- G.poly3
  }
  
  if (pick <= 0){
    lamTry <- cutoff(eigen(G.polyTry)$values, cut)
  } else {
    degree <- 3
    lamTry <- eigen(G.polyTry)$values[pick]
  }
  #print(lamTry)

  ans <- predFunc(G.polyTry, lamTry)
  
  if (plot==T){
    plot(eigen(G.polyTry)$values[1:20], main=paste("Lambda of Poly D=", toString(degree), sep=""), 
         ylab='Lambda')
    lines(eigen(G.polyTry)$values[1:20])
    abline(h=lamTry, col='red')
    text(10, lamTry, toString(round(lamTry,4)), col = "red")
  } 
  return(ans)
}

poly.TryFunc()
poly.TryFunc(2, pick=21)
cutoff.res <- eigen(G.poly2)$values[18]/eigen(G.poly2)$values[1]; cutoff.res  # >8% cutoff
poly.TryFunc(2, pick=4)

poly.TryFunc(3, pick=56)
cutoff.res <- eigen(G.poly3)$values[56]/eigen(G.poly3)$values[1]; cutoff.res  # >5% cutoff
poly.TryFunc(3, pick=3)

tryall.poly2 <- c()
for (i in c(1:21)){  # Max 21, limit at 22
  tryall.poly2 <- c(tryall.poly2, poly.TryFunc(degree=2, pick=i))
}
#tryall.poly2
plot(tryall.poly2, ylab='Performance', xlab='Cutoff Index', main='Poly 2 Performance Plot')
lines(tryall.poly2, col='red')

tryall.poly3 <- c()
for (i in c(1:21)){  # Max 56, limit at 57
  tryall.poly3 <- c(tryall.poly3, poly.TryFunc(degree=3, pick=i))
}
#tryall.poly3
plot(tryall.poly3, ylab='Performance', main='Performance Comparison')
lines(tryall.poly3, col='green')
lines(tryall.poly2, col='red')


'--------------------------------------------------------'
# Trying lambda and gamma for exp
#gamma.Try <- 0.2
#cutoff.Try <- 9.99
#G.expTry <- expFunc(gamma.Try)
#eigen(G.expTry)$values[1:10] # show lambda selections
#lamTry <- cutoff(eigen(G.expTry)$values, cutoff.Try)
#lamTry <- eigen(G.expTry)$values[27]
#predFunc(G.expTry, lamTry)

exp.TryFunc <- function(gamma=0.2, cut=9.99, pick=-1, plot=F){
  G.expTry <- expFunc(gamma)
  if (pick <= 0){
    lamTry <- cutoff(eigen(G.expTry)$values, cut)
  } else {
    lamTry <- eigen(G.expTry)$values[pick]
  }
  ans <- predFunc(G.expTry, lamTry)
  
  if (plot==T){
    plot(eigen(G.expTry)$values[430:447], main=paste("Lambda of Exp G=", toString(gamma), sep=""), 
         ylab='Lambda')
    lines(eigen(G.expTry)$values[430:447])
    abline(h=lamTry, col='red')
    text(10, lamTry, labels=toString(round(lamTry,4)), col = "red") 
    abline(v=9, col='green')
  }
  return(ans)
}

exp.TryFunc(14800, 9.99, plot=T)
exp.TryFunc(0.2, cut=7, plot=T)

tryall.exp <- c()
for (i in c(1:16)){  # No.17 is very high -> outlier
  tryall.exp <- c(tryall.exp, exp.TryFunc(gamma=0.2, pick=i))
}
#tryall.exp
plot(tryall.exp, ylab='Performance', main='Exp Performance Plot')
lines(tryall.exp, col='blue')
lines(tryall.poly3, col='green')
lines(tryall.poly2, col='red')


exp.TryFunc(0.2, pick=2)

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'





# Step 6. Record all Perf for each kernel mtds and diff Gammas & Lambdas in a Table
# Choose the best ones.
# A copy from above
poly.TryFunc()
poly.TryFunc(2, pick=21)
cutoff.res <- eigen(G.poly2)$values[4]/eigen(G.poly2)$values[1]; cutoff.res  # >8% cutoff
poly.TryFunc(2, pick=4)

poly.TryFunc(3, pick=56)
cutoff.res <- eigen(G.poly3)$values[3]/eigen(G.poly3)$values[1]; cutoff.res  # >5% cutoff
poly.TryFunc(3, pick=3)
'==================================================<<<<<<<<<<<<<<<<<<<<<'
eigen.exp14800 <- eigen(expFunc(14788))$values
lam.14800 <- cutoff(eigen.exp14800); lam.14800
which(eigen.exp14800==lam.14800)
predFunc(expFunc(14788), lam.14800)

cutoff.res <- eigen(expFunc(14800))$values[444]/eigen(expFunc(14800))$values[1]; cutoff.res  # >5% cutoff
exp.TryFunc(14800, 9.99, plot=T)

exp.TryFunc(14788, 9.99, plot=T)

#df.EGout <- c()

#df.EGout <- cbind(df.EGout, G.14790=eigen(expFunc(14790))$values[430:447])
#df.EGout <- data.frame(df.EGout)


#df.PFout <- c()

#perflist <- c()
for (i in eigen(expFunc(14790))$values[430:447]){
  perflist <- c(perflist, predFunc(expFunc(14790), i))
}

plot(df.PFout$G.14788, ylab='Performance', main='Exp Performance Plot')
lines(df.PFout$G.14788, col='blue')

#df.PFout <- cbind(df.PFout, G.14790=perflist)
#df.PFout <- data.frame(df.PFout)
#rownames(df.PFout) <- c(430:447)

#write.csv(df.EGout, file='EigenOut.csv')
#write.csv(df.PFout, file='PerfOut.csv')


#df.EGpoly <- cbind(Poly2=eigen.poly2[1:20], Poly3=eigen.poly3[1:20])
#df.EGpoly <- data.frame(df.EGpoly)

#df.PFpoly <- cbind(Poly2=tryall.poly2[1:20], Poly3=tryall.poly3[1:20])
#df.PFpoly <- data.frame(df.PFpoly)

#write.csv(df.EGpoly, file='Eigen_Poly.csv')
#write.csv(df.PFpoly, file='Perf_Poly.csv')

'==================================================<<<<<<<<<<<<<<<<<<<<<<'
'Starts here...'
perf.table.cut <- c()
cutlist <- c(9.99, 5)  # Cutoff percent
for (i in cutlist){  # Max 21, limit at 22
  tryall.poly2 <- c(tryall.poly2, poly.TryFunc(degree=2, cut=i))
}
perf.table.cut <- cbind(perf.table.cut, poly2=tryall.exp)
rownames(perf.table.cut) <- cutlist
tryall.poly3 <- c()
for (i in cutlist){  # Max 56, limit at 57
  tryall.poly3 <- c(tryall.poly3, poly.TryFunc(degree=3, cut=i))
}
perf.table.cut <- cbind(perf.table.cut, poly3=tryall.poly3)
rownames(perf.table.cut) <- cutlist
perf.table.cut  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Exp loop 1, by pick i
perf.table <- c()
picklist <- c(1:10)  # Pick position
for (gam in c(0.0001, 0.01)){
  tryall.exp <- c()
  for (i in picklist){  # Num of rows
    tryall.exp <- c(tryall.exp, exp.TryFunc(gamma=gam, pick=i))
  }
  perf.table <- cbind(perf.table, gam=tryall.exp)
  colnames(perf.table)[colnames(perf.table)=='gam'] <- gam
  rownames(perf.table) <- picklist
}
perf.table

# Exp loop 2, by cutoff j
perf.table.cut <- c()
cutlist <- c(9.99)  # Cutoff percent
for (gam in c(3100, 3050, 3210)){
  tryall.exp <- c()
  for (j in cutlist){  # Num of rows
    tryall.exp <- c(tryall.exp, exp.TryFunc(gamma=gam, cut=j))
  }
  perf.table.cut <- cbind(perf.table.cut, gam=tryall.exp)
  colnames(perf.table.cut)[colnames(perf.table.cut)=='gam'] <- gam
  #rownames(perf.table.cut) <- cutlist
}
perf.table.cut  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Exp loop 3, everything before cutoff j (ONE VALUE ONLY!!!)
perf.table.cut.before <- c()
cutlist <- c(9.99)  # Cutoff percent
for (gam in c(500, 1000, 1500, 2000)){
  tryall.exp <- c()
  for (j in cutlist){  # Num of rows
    tryall.exp <- c(tryall.exp, exp.TryFunc.list(gamma=gam, cut=j))
  }
  perf.table.cut.before <- cbind(perf.table.cut.before, gam=tryall.exp)
  colnames(perf.table.cut.before)[colnames(perf.table.cut.before)=='gam'] <- gam
  #rownames(perf.table.cut) <- cutlist
}
perf.table.cut.before 




# Step 7. Using fixed parameters and kernel, re-train Coef (alphas) on a subset training data.
# Compute the new errors of prediction: RMSE_train and RMSE_test, and new performances.
# Repeat 20 times to have 20 diff Perf/RMSE for train/test, so 20 diff Coef (alphas), and 20 perfs.
# Get mean and std dev of the 20 Perf_train/test. Conclusions on prediction reliability.

set.seed(1)
train <- sample(1:nrow(x.ceo), 0.85*nrow(x.ceo))
test <- c(1:nrow(x.ceo))[-train]  # Test set















#expResMax=20


exp.TryFunc.list <- function(gamma=0.2, cut=9.99, pick=-1, plot=F){
  G.expTry <- expFunc(gamma)
  if (pick <= 0){
    lampix <- cutoff(eigen(G.expTry)$values, cut)
    lamTry <- cutoff.before(eigen(G.expTry)$values, cut)
  } else {
    lampix <- eigen(G.expTry)$values[pick]
    lamTry <- eigen(G.expTry)$values[1:pick]
  }
  ans <- c()
  for (ai in lamTry){
    ans <- c(ans, predFunc(G.expTry, ai))
  }
  if (length(ans) <= expResMax){
    ans <- c(ans, rep(0, expResMax-length(ans)))
  }
  
  if (plot==T){
    plot(eigen(G.expTry)$values[1:20], main=paste("Lambda of Exp G=", toString(gamma), sep=""), 
         ylab='Lambda')
    lines(eigen(G.expTry)$values[1:20])
    abline(h=lampix, col='red')
    text(10, lampix, labels=toString(round(lampix,4)), col = "red") 
  }
  return(ans)
}

exp.TryFunc.list(0.2, 9.99)
exp.TryFunc.list(0.2, cut=7, plot=T)

'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
z <- solve(G - lambda * diag(447), y.ceo)

# Leave one out cv
errs <- c()
preds <- c()
for (i in 1:447){
  Gprime <- kernelMatrix(x.ceo[-i,], kernel=polydot(degree=2))
  z <-  solve(Gprime - lambda * diag(446), y.ceo[-i])
  pred <- 0;
  j <- 1;
  while (j < i){
    pred <- pred + z[j] * G[i,j];
    j <- j+1
  }
  while (j < 447){
    pred <- pred + z[j] * G[i,j+1];
    j <- j+1
  }
  errs <- c(errs, pred - y.ceo[i])
  preds <- c(preds, pred)
}
errs
preds <- preds + mean(ceo$salary)
preds
mse <- mean(errs^2)
mse


































































