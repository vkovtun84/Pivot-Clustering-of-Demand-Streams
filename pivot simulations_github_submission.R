
library(dse)
library(matrixcalc)
library(tseries)
library(polynom)
library(PolynomF)

nstreams=25 ## number of streams.
## Simulating ARMA model coefficients. All streams are ARMA(1,1)
phi_list=list()
theta_list=list()

for(k in 1:nstreams){
  phi_1=runif(1,-1,1)
  theta_1=runif(1,-1,1)
  
  phi_list=append(phi_list,list(c(1,-phi_1,0)))
  theta_list=append(theta_list,list(c(1,-theta_1,0)))
}

covarmat=matrix(0)
while(!is.positive.definite(covarmat)){
  n <- nstreams 
  A <- matrix(runif(n^2)*2-1, ncol=n) 
  covarmat <- t(A) %*% A
}

####################

TrueMSFEind=sum(covarmat) ## MSFE when using individual streams with TRUE ARMA coefficients.
TrueMSFEind

zspec=zspec3(phi_list=phi_list,theta_list=theta_list,covarmat=covarmat) ### zspec3
res=InnVar(zspec$Pz,zspec$PzDeg,zspec$Qz,zspec$QzDeg) ## roots of Pz and Qz unstable (complex root conjugates are not same)
aggARMA=getARMAfromInnVar(res) ## Getting equal AR and MA degrees but think this is OK due to rule.
aggARMAall=aggARMA
TrueMSFEagg=res$var
TrueMSFEagg ## MSFE when using full aggregate of streams with TRUE ARMA coefficients.

###########################

Nobs=1500 ## number of observations to be simulated.

arvec=c()
for(i in 1:length(phi_list)) arvec=c(arvec,phi_list[[i]],rep(rep(0,3),length(phi_list)))
arvec=arvec[-(length(arvec)):-(length(arvec)-length(rep(rep(0,3),length(phi_list)))+1)]
AR=array(arvec,c(3,length(phi_list),length(phi_list)))

mavec=c()
for(i in 1:length(theta_list)) mavec=c(mavec,theta_list[[i]],rep(rep(0,3),length(theta_list)))
mavec=mavec[-(length(mavec)):-(length(mavec)-length(rep(rep(0,3),length(theta_list)))+1)]
MA=array(mavec,c(3,length(theta_list),length(theta_list)))

#AR=array(c(phi1,rep(rep(0,3),length(phi_list)),phi2,rep(rep(0,3),length(phi_list)),phi3,rep(rep(0,3),length(phi_list)),phi4),c(3,length(phi_list),length(phi_list)))
#MA=array(c(theta1,rep(rep(0,3),length(theta_list)),theta2,rep(rep(0,3),length(theta_list)),theta3,rep(rep(0,3),length(theta_list)),theta4),c(3,length(theta_list),length(theta_list)))
#for 3 processes: AR=array(c(phi1,rep(0,3),rep(0,3),rep(0,3),phi2,rep(0,3),rep(0,3),rep(0,3),phi3),c(3,3,3))

model=ARMA(AR,MA)
SIM=simulate(model,sampleT=Nobs,Cov=covarmat)  ### SIMULATES Multivariate ARMA data as well as shocks corresponding to the models and covariance matrix above.
ProcessData=SIM$output

shocks=SIM$noise
ShockData=shocks$w
ShockCOV=var(ShockData) 
sum(ShockCOV)
############################


res_clusters=c()
order=c(1,0,1)
#for(i in 1:ncol(ProcessData)) res_clusters=cbind(res_clusters,arma(ProcessData[,i],order=order)$residuals[(max(1,order)+1):Nobs])
for(i in 1:ncol(ProcessData)) res_clusters=cbind(res_clusters,arima(ProcessData[,i],order=order)$residuals) ## arima() seems to do a better job for small order 
EstMSFEind=sum(cov(res_clusters))
EstMSFEind


n.clusters=1
Groups=sample(1:n.clusters,nstreams,replace=TRUE)
ClData=t(rowsum(t(ProcessData),Groups))

res_clusters=c()
order=c(5,0,5) #order=c(5,5)
#for(i in 1:ncol(ClData)) res_clusters=cbind(res_clusters,arma(ClData[,i],order=order)$residuals[(max(1,order)+1):Nobs])
for(i in 1:ncol(ClData)) res_clusters=cbind(res_clusters,arima(ClData[,i],order=order,method="CSS")$residuals)
EstMSFEagg=sum(cov(as.matrix(res_clusters)))
EstMSFEagg


#################### 50 iterations from different random group assignments: ###############
n.clusters=4
RandGroupsList=list()

TrueMSFEsub.rand=c()
TrueMSFEsub.piv=c()
TrueGroup.piv=list()

EstMSFEsub.rand=c()
EstMSFEsub.piv=c()
EstGroup.piv=list()

start=proc.time()
for(k in 1:2){ ## 1:50 MAKE THIS LARGER
  RandGroups=sample(1:n.clusters,nstreams,replace=TRUE)
  RandGroupsList[[k]]=RandGroups
  NewRandMSFEtrue=MSFEsubaggTRUE(phi_list,theta_list,covarmat,RandGroups)
  TrueMSFEsub.rand=c(TrueMSFEsub.rand,NewRandMSFEtrue)
  
  ### ALL arma(ClData[,i],order=order)$residuals[(max(1,order)+1):Nobs] changed to arima(ClData[,i],order=order)$residuals on 12/11/18. Also changed in pivotMSFE2().
  ClData=t(rowsum(t(ProcessData),RandGroups))
  res_clusters=c()
  order=c(5,0,5)
  #for(i in 1:ncol(ClData)) res_clusters=cbind(res_clusters,arma(ClData[,i],order=order)$residuals[(max(1,order)+1):Nobs])
  for(i in 1:ncol(ClData)) res_clusters=cbind(res_clusters,arima(ClData[,i],order=order,method="CSS")$residuals)
  NewRandMSFEest=sum(cov(as.matrix(res_clusters)))
  EstMSFEsub.rand=c(EstMSFEsub.rand,NewRandMSFEest)
  
  #### Local minimum (TRUE MSFE):   ### may be a problem on next line: sum(cov(res_clusters)) is estimated MSFE!
  pivold=pivotMSFE_TRUE(phi_list,theta_list,covarmat,NewRandMSFEtrue,RandGroups,unfin.clust=1:max(RandGroups)) ## MSFE of latest random: sum(cov(res_clusters))
  newMSFE_true=pivold[[2]]
  newGroups=pivold[[1]]
  oldMSFE_true=0
  while(abs(oldMSFE_true-newMSFE_true)>.0000001){
    oldMSFE_true=newMSFE_true
    res=pivotMSFE_TRUE(phi_list,theta_list,covarmat,newMSFE_true,newGroups,unfin.clust=1:max(newGroups))
    newMSFE_true=res[[2]]
    newGroups=res[[1]]
  }
  TrueGroup.piv[[k]]=newGroups
  TrueMSFEsub.piv=c(TrueMSFEsub.piv,newMSFE_true)
  
  
  #### Local minimum (ESTIMATED MSFE):
  
  Estpivold=pivotMSFE2(ProcessData,MSFEcurr=NewRandMSFEest,RandGroups,unfin.clust=1:max(RandGroups),ar.order=5,ma.order=5)
  newMSFE_est=Estpivold[[2]]
  newGroups_est=Estpivold[[1]]
  oldMSFE_est=0
  
  while(abs(oldMSFE_est-newMSFE_est)>.0000001){
    oldMSFE_est=newMSFE_est
    res=pivotMSFE2(ProcessData,newMSFE_est, newGroups_est,unfin.clust=1:max(newGroups_est),ar.order=5,ma.order=5)
    newMSFE_est=res[[2]]
    newGroups_est=res[[1]]
  }
  EstGroup.piv[[k]]=newGroups_est
  EstMSFEsub.piv=c(EstMSFEsub.piv,newMSFE_est)
  
}
proc.time()-start


TrueMSFEsub.rand
TrueMSFEsub.piv
TrueGroup.piv

EstMSFEsub.rand
EstMSFEsub.piv
EstGroup.piv

## Compare to True Disaggregated MSFE:
TrueMSFEind


########### The smallest MSFE groups found in both:
nsims=length(TrueGroup.piv)
TrueGroup.piv[[ min((1:nsims)[(TrueMSFEsub.piv-min(TrueMSFEsub.piv))^2<.0001])]]

EstGroup.piv[[ min((1:nsims)[(EstMSFEsub.piv-min(EstMSFEsub.piv))^2<.0001])]]


######################
##Estimated MSFE of 1st True-Group Pivot assignment:
Groups=TrueGroup.piv[[1]]
ClData=t(rowsum(t(ProcessData),Groups))
res_clusters=c()
order=c(5,0,5)
for(i in 1:ncol(ClData)) res_clusters=cbind(res_clusters,arma(ClData[,i],order=order)$residuals[(max(1,order)+1):Nobs])
sum(cov(res_clusters))
#########################
Groups=EstGroup.piv[[1]]
MSFEsubaggTRUE(phi_list,theta_list,covarmat,Groups)

##############################
#### The following defines a distance measure based on how often streams are clustered together
#### in the pivot groupings above. 
#### We then use heirarchical clustering to cluster streams based on this distance measure as a final step.

simmat_true=matrix(rep(0,nstreams^2),ncol=nstreams)

for(i in 1:nsims){
  for(k in 1:n.clusters){
    Clust=(1:nstreams)[TrueGroup.piv[[i]]==k]
    simmat_true[Clust,Clust]=simmat_true[Clust,Clust]+1 ## contains the number of times the streams were clustered together in the 40 simulations.
  }  
}
dist_true=nsims-simmat_true

Clusters=hclust(as.dist(dist_true),method="complete")
True_Fin_clust=cutree(Clusters,k=3)


simmat_est=matrix(rep(0,nstreams^2),ncol=nstreams)
for(i in 1:nsims){
  for(k in 1:4){
    Clust=(1:nstreams)[EstGroup.piv[[i]]==k]
    simmat_est[Clust,Clust]=simmat_est[Clust,Clust]+1 ## contains the number of times the streams were clustered together in the 40 simulations.
  }  
}
dist_est=nsims-simmat_est
Clusters=hclust(as.dist(dist_est),method="complete")
Est_Fin_clust=cutree(Clusters,k=4)

True_Fin_clust
Est_Fin_clust



############################


############### Alternative way to use 


(1:nstreams)[simmat_true[,1]>15]   
(1:nstreams)[simmat_true[,2]>15]   
(1:nstreams)[simmat_true[,5]>15]   
(1:nstreams)[simmat_true[,10]>15]   


(1:nstreams)[simmat_est[,1]>15]   
(1:nstreams)[simmat_est[,2]>15]   
(1:nstreams)[simmat_est[,4]>15]   
(1:nstreams)[simmat_est[,5]>15]     


##########################################
################################################# Use this:






