#required package to compute Newey and West variance estimator
library(sandwich)

############################################################################################################## 
############################################################################################################## 
#functions for Gini's mean difference linearization
ffun<-function(x,y){abs(x-y)}
#Delta_M computes the matrix of values |x_i-x_j| for all i,j necessary to perform linearization of Delta
Delta_M <- function(x){
  outer(x, x, ffun);
}

############################################################################################################## 
###########    inference on K2  
############################################################################################################## 

#K2_TS_mu: starting from a time series x, this functions computes:
# 1) the estimator K2hat of the kurtosis index K2
# 2) the estimator of asymptotic variance of K2hat under the IID assumption (\widehat{Vhat}_{Kn}^M in the paper)
# 3) the estimator of asymptotic variance of K2hat estimator (\widehat{Vhat}_{Kn} in the paper) based on the Newey and West technique 
#    applied to linearized statistics (formulas 20 and 21 in the main paper)
# 4) the bandwidth used in performing estimation at point 3)

K2_TS_mu <- function(x){
  #sample size
  n <- length(x)
  #POINT ESTIMATION
  m1 <- mean(x)
  delta <- mean(abs(x-m1))
  matrice <- Delta_M(x)
  Delta <- sum(matrice) / (n * (n-1))
  #index estimate
  r1<-Delta/delta-1
  
  #LINEARIZATIONS FOR VARIANCE ESTIMATION 
  #linearization of numerator (Gini's Delta) (g_1(x) in the proof of theorem 4 in the paper)
  y_n<-2 * rowSums(matrice) / n
  #linearization of denominator (mean absolute deviation) (g_2(x) in the proof of theorem 4 in the paper)
  y_d <- 2 * (x - m1) * ((length(x[x<=m1]) / n) - (x <= m1))
  #linearization of the statistics (formula 13, h_K(x) in the main paper). Constant
  #d and Delta are not included in the linearizations because they are constants and
  #then do not influence variance computation
  y <- (1/delta) * y_n  - (Delta / delta^2)* y_d
  
  #VARIANCE ESTIMATION UNDER IID ASSUMPTION
  r2 <- var(y) 
  
  #VARIANCE ESTIMATION
  #creation of the object of class "lm". That object is created since function NeweyWest of 
  #package "sandwich" works only with arguments of such a class.
  #A constant regression is perform in order to produce residuals which are simple translation of original data
  fm <- lm(y ~ 1)
  #Newey-West variance estimates with prewhitening and optimal bandwidth computation (see the details in the NeweyWest function documentation)
  r3 <- NeweyWest(fm)*n
  #optimal bandwith used in producing estimate r3
  r4<-bwNeweyWest(fm)
  
  res<-c(r1,r2,r3,r4)
  names(res)<-c("K2hat","Var_IID", "Var_NW","BW_NW")
  return(res)
}


############################################################################################################## 
###########    inference on Beta2  
############################################################################################################## 

#Beta2_TS: starting from a time series x, this functions computes:
# 1) the estimator Beta2hat of the kurtosis index Beta2
# 2) the estimator of asymptotic variance of Beta2hat under the IID assumption (\widehat{Vhat}_{Bn}^M in the paper)
# 3) the estimator of asymptotic variance of Beta2hat estimator (\widehat{Vhat}_{Bn} in the paper) based on the Newey and West technique 
#    applied to linearized statistics (formulas 18 and 29 in the main paper)
# 4) the bandwidth used in performing estimation at point 3)

Beta2_TS <- function(x){
  #sample size
  n <- length(x)
  #POINT ESTIMATION
  m1 <- mean(x)
  m2c<-mean((x-m1)^2)
  sigma<-sqrt(m2c)
  m3c<-mean((x-m1)^3)
  m4c<-mean((x-m1)^4)
  #index estimate
  r1<-m4c/m2c^2
  
  #LINEARIZATIONS FOR VARIANCE ESTIMATION 
  #linearization of numerator (fourth moment estimator) (g_2(x) in the proof of theorem 3 in the paper)
  y_n<-(x-m1)^4-4*m3c*(x-m1)
  #linearization of denominator (second moment estimator) (g_1(x) in the proof of theorem 3 in the paper)
  y_d <- (x - m1)^2
  #linearization of the statistics (formula 9, h_B(x) in the main paper). Constant
  #mu_4 and sigma^2 are not included in the linearizations because they are constants and
  #then do not influence variance computation
  y <- (1/m2c^2) * y_n  - (2*m4c / m2c^3) * y_d
  
  #VARIANCE ESTIMATION UNDER IID ASSUMPTION
  r2 <- var(y) 
  
  #VARIANCE ESTIMATION
  #creation of the object of class "lm". That object is created since function NeweyWest of 
  #package "sandwich" works only with arguments of such a class.
  #A constant regression is perform in order to produce residuals which are simple translation of original data
  fm <- lm(y ~ 1);
  #Newey-West variance estimates with prewhitening and optimal bandwidth computation (see the details in the NeweyWest function documentation)
  r3 <- NeweyWest(fm)*n;
  #optimal bandwith used in producing estimate r3
  r4<-bwNeweyWest(fm)
  
  res<-c(r1,r2,r3,r4)
  names(res)<-c("Beta2hat","Var_IID", "Var_NW","BW_NW")
  return(res)
}
