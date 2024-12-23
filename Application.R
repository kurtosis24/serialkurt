
################################################################################
#############      APPLICATION        ##########################################
################################################################################
source("main_functions.r")

#### evaluation of kurtosis on financial returns

library(quantmod)

#ticker of financial assets analized (all tiker in US S&P500)
Ticker<-as.vector(read.table(file="Ticker.txt"))
#1-alpha= confidence level of computed confidence intervals
#alpha=test level
alpha<-0.05
#matrix in which results are reported: 
#one row for each Ticker
#colums: index estimator, variance estiamator under iid assumption, variance estimator using Newey West, bandwidth of Newey west estimator,
#lower bound of ocnfidence interval, upper bound of confidence interval, sample size
K2Res<-K2Res_Adj<-matrix(0,length(Ticker$V1),7)
Beta2Res<-Beta2Res_Adj<-matrix(0,length(Ticker$V1),7)
#Ticker with outliers
outliers<-NULL
#position in the list of Tickers of assets with outliers
outliers_index<-NULL
#number of outliers in the time series for Ticket with outliers
outliers_number<-NULL
for(i in 1:length(Ticker$V1)){
  #Ticker ididentificator
  name<-Ticker$V1[i]
  #downloading of Ticker data
  getSymbols(name, src = "yahoo", from = "2022-12-03", to = "2024-12-02",return.class='matrix')
  #extrapolation of Adjusted Close Prices
  AdjClose<-get(name)[,6]
  #removing missing values 
  rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)),class) == "data.frame"])
  #log-returns computation
  returns<-diff(log(na.omit(AdjClose)))        
  #log-returns computation with outliers correction
  #log-returns greater than 15% and lower than -15% are substituted by the average
  #between their antecedent and subsequent returns. If the outlier coincides whit the
  #first return, it is set equal to the second. If the outlier coincides whit the
  #last return, it is set equal to the second-last.  
  returns_Adj<-returns
  if(sum(abs(returns_Adj)>0.15)>0){
    outliers<-cbind(outliers, name) 
    outliers_index<-cbind(outliers_index,i) 
    outliers_number<-cbind(outliers_number,sum(abs(returns_Adj)>0.15)) 
  }
  while(sum(abs(returns_Adj)>0.15)>0){
    pos<-which(abs(returns_Adj)>0.15)
    for(j in pos){
      if(j==1){returns_Adj[j]<-returns_Adj[j+1]} else{
        if(j==length(returns_Adj)){returns_Adj[j]<-returns_Adj[j-1]} else{
          returns_Adj[j]<-(returns_Adj[j-1]+returns_Adj[j+1])/2
        }
      }
    }
  }
  #print to video in order to contol the state of progress of computation
  print(c(i,sum(abs(returns)>0.15),sum(abs(returns_Adj)>0.15)))  
  #save sample size
  n<-K2Res[i,7]<-Beta2Res[i,7]<-K2Res_Adj[i,7]<-Beta2Res_Adj[i,7]<-length(returns)
  #point estimates and varince estimation on observed series of log-returns
  K2Res[i,1:4]<-K2_TS_mu(returns)
  Beta2Res[i,1:4]<-Beta2_TS(returns)
  #lower and upper bounds of confidence intervals from the observed series of log-returns 
  K2Res[i,5]<-K2Res[i,1]-qnorm(1-alpha/2)*(K2Res[i,3]/n)^0.5
  K2Res[i,6]<-K2Res[i,1]+qnorm(1-alpha/2)*(K2Res[i,3]/n)^0.5
  Beta2Res[i,5]<-Beta2Res[i,1]-qnorm(1-alpha/2)*(Beta2Res[i,3]/n)^0.5
  Beta2Res[i,6]<-Beta2Res[i,1]+qnorm(1-alpha/2)*(Beta2Res[i,3]/n)^0.5
  #point estimates and varince estimation on the series of log-returns with outliers correction
  K2Res_Adj[i,1:4]<-K2_TS_mu(returns_Adj)
  Beta2Res_Adj[i,1:4]<-Beta2_TS(returns_Adj)
  #lower and upper bounds of confidence intervals from the series of log-returns with outliers correction
  K2Res_Adj[i,5]<-K2Res_Adj[i,1]-qnorm(1-alpha/2)*(K2Res_Adj[i,3]/n)^0.5
  K2Res_Adj[i,6]<-K2Res_Adj[i,1]+qnorm(1-alpha/2)*(K2Res_Adj[i,3]/n)^0.5
  Beta2Res_Adj[i,5]<-Beta2Res_Adj[i,1]-qnorm(1-alpha/2)*(Beta2Res_Adj[i,3]/n)^0.5
  Beta2Res_Adj[i,6]<-Beta2Res_Adj[i,1]+qnorm(1-alpha/2)*(Beta2Res_Adj[i,3]/n)^0.5
}
rownames(Beta2Res)<-rownames(K2Res)<-rownames(Beta2Res_Adj)<-rownames(K2Res_Adj)<-Ticker$V1
colnames(Beta2Res)<-colnames(Beta2Res_Adj)<-c("Beta2hat","Var_IID", "Var_NW","BW_NW", "IC_lowerbound", "IC_upperbound", "sample size")
colnames(K2Res)<-colnames(K2Res_Adj)<-c("K2hat","Var_IID", "Var_NW","BW_NW", "IC_lowerbound", "IC_upperbound", "sample size")

TAB<-cbind(Beta2Res[,1],Beta2Res_Adj[,1], K2Res[,1], K2Res_Adj[,1])
colnames(TAB)<-c("Beta2", "Beta2_Adj","K2", "K2_Adj")
rownames(TAB)<-Ticker$V1
write.table(TAB, file="Stime_puntuali.txt", sep=";", dec=",")

################################################################################
#Analisys of GL (GLOBAL LIFE INC.)
################################################################################ 

getSymbols(Ticker$V1[469], src = "yahoo", from = "2022-12-03", to = "2024-12-02",return.class='matrix')
AdjClose<-get(Ticker$V1[469])[,6]
returns<-diff(log(na.omit(AdjClose)))  
write.table(returns, file="rendimenti_GL.txt", sep=";", dec=",")
which(abs(returns)>0.15)
returns_Adj<-returns
#data of observed outliers
rownames(get(Ticker$V1[469]))[which(abs(returns_Adj)>0.15)+1]
#outliers value
returns_Adj[which(abs(returns_Adj)>0.15)]
#outliers correction
while(sum(abs(returns_Adj)>0.15)>0){
  pos<-which(abs(returns_Adj)>0.15)
  for(j in pos){
    if(j==1){returns_Adj[j]<-returns_Adj[j+1]} else{
      if(j==length(returns_Adj)){returns_Adj[j]<-returns_Adj[j-1]} else{
        returns_Adj[j]<-(returns_Adj[j-1]+returns_Adj[j+1])/2
      }
    }
  }
}
write.table(returns_Adj, file="rendimenti_Adj_GL.txt", sep=";", dec=",")
#results table
RES<-matrix(0,7,4)
RES<-cbind(K2Res[469,],K2Res_Adj[469,],Beta2Res[469,],Beta2Res_Adj[469,])
RES[7,]<-RES[6,]-RES[5,]
colnames(RES)<-c("K2", "K2_Adj","Beta2", "Beta2_Adj")
rownames(RES)<-c("index estimate", "variance estimate_IID", "variance estimate", "Bw_NW", "lower bound", "upper bound", "length")

################################################################################
#Table 11 in the paper
RES
################################################################################


################################################################################
#plot of observed returns (figure 2-left in the paper)
x11()
plot(returns, type="l", xlab=" ", ylab=" ")
#plot of adjusted returns (figure 2-right in the paper)
x11()
plot(returns_Adj, type="l", xlab=" ", ylab=" ")
################################################################################

################################################################################
#Table 12 in the paper
p<-outliers_index[which(outliers=="TRMB")]
RES<-matrix(0,7,4)
RES<-cbind(K2Res[p,],K2Res_Adj[p,],Beta2Res[p,],Beta2Res_Adj[p,])
RES[7,]<-RES[6,]-RES[5,]
colnames(RES)<-c("K2", "K2_Adj","Beta2", "Beta2_Adj")
rownames(RES)<-c("index estimate", "variance estimate_IID", "variance estimate", "Bw_NW", "lower bound", "upper bound", "length")

RES
################################################################################

################################################################################
#number of Ticker with outliers
length(outliers)
#diustribution of the number of outliers in the series
table(outliers_number)
################################################################################

################################################################################
#effects on outliers adjustment on poit estimates of kurtosis indexes
#Plots in Figure 3 of the main paper
################################################################################
op<-par()

x11()
par(mai=c(1,1,0.5,0.5))
plot(Beta2Res[-469,1], Beta2Res_Adj[-469,1],col=0,
     xlab=expression(paste(hat(beta)[paste(2,"n")], " without outliers adjustement")),
     ylab=expression(paste(hat(beta)[paste(2,"n")], " with outliers adjustement")), 
     xlim=c(min(Beta2Res[,1]), max(Beta2Res[,1])+20),      
     ylim=c(min(Beta2Res_Adj[,1]), max(Beta2Res_Adj[,1])),
     pch=20, cex=0.8)
text(Beta2Res[469,1], Beta2Res_Adj[469,1],labels=Ticker$V1[469], pos=4)
abline(a=0, b=1, col="darkgray")
points(Beta2Res[-469,1], Beta2Res_Adj[-469,1],col=1,pch=20, cex=0.8)

op<-par()

x11()
par(mai=c(1,1,0.5,0.5))
plot(K2Res[-469,1], K2Res_Adj[-469,1],col=0,
     xlab=expression(paste(hat(K)[paste(2,"n")], " without outliers adjustement")),
     ylab=expression(paste(hat(K)[paste(2,"n")], " with outliers adjustement")), 
     xlim=c(min(K2Res[,1]), max(K2Res[,1])+0.02),      
     ylim=c(min(K2Res_Adj[,1]), max(K2Res_Adj[,1])),
     pch=20, cex=0.8)
text(K2Res[469,1], K2Res_Adj[469,1],labels=Ticker$V1[469], pos=4)
abline(a=0, b=1, col="darkgray")
points(K2Res[-469,1], K2Res_Adj[-469,1],col=1,pch=20, cex=0.8)


################################################################################
#Analisys of Test Results when outliers are present
################################################################################
#Ticker with outliers
ii<-outliers_index

#p value computation
#value under the null hypotesis
k0<-2^0.5-1
b0<-3

#computation of test statistics on observed returns
ST_K2<-abs(K2Res[ii,1]-k0)/(K2Res[ii,3]/K2Res[ii,7])^0.5
ST_Beta2<-abs(Beta2Res[ii,1]-b0)/(Beta2Res[ii,3]/Beta2Res[ii,7])^0.5

#computation of test statistics on returns with outliers substitution
ST_K2_Adj<-abs(K2Res_Adj[ii,1]-k0)/(K2Res_Adj[ii,3]/K2Res_Adj[ii,7])^0.5
ST_Beta2_Adj<-abs(Beta2Res_Adj[ii,1]-b0)/(Beta2Res_Adj[ii,3]/Beta2Res_Adj[ii,7])^0.5

#p_value computation on observed returns
p_value_K2<-2*(1-pnorm(ST_K2))
p_value_Beta2<-2*(1-pnorm(ST_Beta2))

#p_value computation on returns with outliers substitution
p_value_K2_Adj<-2*(1-pnorm(ST_K2_Adj))
p_value_Beta2_Adj<-2*(1-pnorm(ST_Beta2_Adj))

cbind(p_value_Beta2,p_value_Beta2_Adj,p_value_K2, p_value_K2_Adj,(p_value_Beta2<=0.05)*(p_value_Beta2_Adj>0.05),(p_value_K2<=0.05)*(p_value_K2_Adj>0.05) )

#Agreement table at level lim for K2
lim<-0.05
t<-(p_value_K2<=lim)
t_Adj<-(p_value_K2_Adj<=lim)
TAB_K2<-matrix(0,3,3)
TAB_K2[1,1]<-sum((t==0)*(t_Adj==0))
TAB_K2[1,2]<-sum((t==0)*(t_Adj==1))
TAB_K2[2,1]<-sum((t==1)*(t_Adj==0))
TAB_K2[2,2]<-sum((t==1)*(t_Adj==1))
TAB_K2[,3]<-rowSums(TAB_K2)
TAB_K2[3,]<-colSums(TAB_K2)
colnames(TAB_K2)<-c("DNR_Adj", "R_Adj", "Total")
rownames(TAB_K2)<-c("DNR", "R", "Total_Adj")

t<-(p_value_Beta2<=lim)
t_Adj<-(p_value_Beta2_Adj<=lim)
TAB_Beta2<-matrix(0,3,3)
TAB_Beta2[1,1]<-sum((t==0)*(t_Adj==0))
TAB_Beta2[1,2]<-sum((t==0)*(t_Adj==1))
TAB_Beta2[2,1]<-sum((t==1)*(t_Adj==0))
TAB_Beta2[2,2]<-sum((t==1)*(t_Adj==1))
TAB_Beta2[,3]<-rowSums(TAB_Beta2)
TAB_Beta2[3,]<-colSums(TAB_Beta2)
colnames(TAB_Beta2)<-c("DNR_Adj", "R_Adj", "Total")
rownames(TAB_Beta2)<-c("DNR", "R", "Total_Adj")

################################################################################
#Table 13 in the paper
TAB_K2
TAB_Beta2
################################################################################
