library(polynom)
library(PolynomF)


polymult2=function(poly1,deg1,poly2,deg2){
  
  ## This function multiplies two Laurent polynomials together.
  ## poly1 and poly2 are vectors of coefficients (according to ascending powers). 
  ## deg1 and deg2 are the degree of poly1 and poly2
  ## Note that the output may contain zeros at end. The actual degree of the product is (Output Degree)-(number of zeros at end).
  #####ie. if polynomial1=x^(-1)+2+x then poly1=c(1,2,1) and deg1=1
  
  while(poly1[length(poly1)]==0 & length(poly1)>1){poly1=poly1[-length(poly1)] ; deg1=deg1-1} ## Added Feb8, 2019
  while(poly2[length(poly2)]==0 & length(poly2)>1){poly2=poly2[-length(poly2)] ; deg2=deg2-1}
  
  deg=deg1+deg2
  if(deg1>deg2){
    poly2=c(poly2,rep(0,deg1-deg2))
    deg2=deg1
  }
  
  else if(deg1<deg2){
    poly1=c(poly1,rep(0,deg2-deg1))
    deg1=deg2
  }
  
  if(length(poly1)>length(poly2)) poly2=c(rep(0,length(poly1)-length(poly2)),poly2)
  if(length(poly2)>length(poly1)) poly1=c(rep(0,length(poly2)-length(poly1)),poly1)
  
  X=poly1%o%poly2
  
  Y=c()
  Y=X[nrow(X):1,]
  Y=as.matrix(Y)
  
  #for(i in nrow(X):1){
  #Y=rbind(Y,X[i,])
  #}
  
  res=c()
  
  for(k in -(nrow(X)-1):(nrow(X)-1)){
    res=c(res,sum(Y[row(Y) == (col(Y) - k)]))
  }
  
  if(sum(res==0) == length(res)) { ### added Feb 8, 2019
    res=0
    deg=0
  }
  
  else{
    while(res[length(res)]==0){ ### added Apr17,2018
      res=res[c(-length(res))]
      #deg=deg-1 ## added in Apr17, 2018
    }
    while(res[1]==0){ ### added Apr17,2018
      res=res[c(-1)]
      #deg=deg-1 ## added in Apr17, 2018
    }
  }
  
  list(res,deg) ### res gives the coefficients of the product of the two polynomials
  ## a list of size 2: polynomial and its degree
}






############################################################################################################################


polyadd=function(poly1,deg1,poly2,deg2){
  ## This function adds two polynomials together.
  ## poly1 and poly2 are vectors of coefficients (according to ascending powers). 
  ## deg1 and deg2 are the degree of poly1 and poly2
  ## Note that the output may contain zeros at end. The actual degree of the product is (Output Degree)-(number of zeros at end).
  #####ie. if polynomial1=x^(-1)+2+x then poly1=c(1,2,1) and deg1=1
  
  if(deg1>deg2){
    poly2=c(poly2,rep(0,deg1-deg2))
    #poly2=c(rep(0,length(poly1)-length(poly2)),poly2)
    deg2=deg1
  }
  if(deg1<deg2){
    poly1=c(poly1,rep(0,deg2-deg1))
    #poly1=c(rep(0,length(poly2)-length(poly1)),poly1)
    deg1=deg2
  }
  
  if(length(poly1)<length(poly2)){
    poly1=c(rep(0,length(poly2)-length(poly1)),poly1)
    #deg1=deg1+length(poly2)-length(poly1)
  }
  
  if(length(poly2)<length(poly1)){
    poly2=c(rep(0,length(poly1)-length(poly2)),poly2)
    #deg2=deg2+length(poly1)-length(poly2)
  }
  
  sum=poly1+poly2
  list(sum,deg1) ## a list of size 2: polynomial and its degree.
}

############################################################################################################################

ratioadd=function(num1,num1deg,den1,den1deg,num2,num2deg,den2,den2deg){
  
  
  if(sum(num1!=0)==0){
    numfin=list(num2,num2deg)
    denfin=list(den2,den2deg)
  }
  else if(sum(num2!=0)==0){
    numfin=list(num1,num1deg)
    denfin=list(den1,den1deg)
  }
  
  else{
    res1=polymult2(num1,num1deg,den2,den2deg)
    res2=polymult2(den1,den1deg,num2,num2deg)
    res3=polymult2(den1,den1deg,den2,den2deg)
    
    ## For each of the above [[1]] contains coefficients in polynomial and [[2]] contains the degree
    
    numfin=polyadd(res1[[1]],res1[[2]],res2[[1]],res2[[2]])
    denfin=res3
  }
  
  append(numfin,denfin) ## a list of size 4: numerator polynomial, its degree, denominator polynomial, its degree
}


############################################################################################################################


zspec2=function(phi_list, theta_list, covarmat){
  ######## DO NOT USE THIS. USE zspec3() instead!!!!
  
  ## theta_list contains theta1, theta2, theta3, etc. for as many ARMAs as we are adding together.
  ## phi_list contains phi1, phi2, phi3, etc. for as many ARMAs as we are adding.
  ## These are lists containing vectors of coefficients according to convention:
  ##### the polynomials in the ARMA expression use the convention phi(z)=phi[1]+phi[2]z+phi[3]z^2+...
  ##### Thus if phi(z)=1+2z then phi1=c(1,2)
  ## covarmat is the variance-covariance matrix of the white-noise sequences.
  
  
  currratio=list(c(0),c(0),c(1),c(0))
  currratio2=list(c(0),c(0),c(1),c(0))
  for(curr in 1:length(theta_list)){
    for(i in 1:length(theta_list)){
      
      #if(curr==2 & i==2) browser()
      if(curr==i){
        thetas=polymult2(covarmat[curr,i]*theta_list[[curr]],deg1=length(theta_list[[curr]])-1,rev(theta_list[[i]]),deg2=0)
        phis=polymult2(phi_list[[curr]],deg=length(phi_list[[curr]])-1,rev(phi_list[[i]]),deg2=0)
        
        currratio=ratioadd(num1=thetas[[1]],num1deg=thetas[[2]],den1=phis[[1]],den1deg=phis[[2]], num2=currratio[[1]], num2deg=currratio[[2]], den2=currratio[[3]], den2deg=currratio[[4]])
        
        Pz=currratio[[1]]
        PzDeg=currratio[[2]]
        Qz=currratio[[3]]
        QzDeg=currratio[[4]]
        
        if(length(Pz)>1){
          while(Pz[length(Pz)]==0){
            Pz=Pz[c(-length(Pz))]
            PzDeg=PzDeg-1
          }
          while(Pz[1]==0){
            Pz=Pz[c(-1)]
          }
        }
        
        if(length(Qz)>1){
          while(Qz[length(Qz)]==0){
            Qz=Qz[c(-length(Qz))]
            QzDeg=QzDeg-1
          }
          while(Qz[1]==0){
            Qz=Qz[c(-1)]
          }
        }
        
        currratio=list(Pz,PzDeg,Qz,QzDeg) ## This is for the non-cross spectra
      }
      
      ##################
      nums=list()
      if(curr!=i){ #### #THIS NEEDS TO BE FIXED!!!! THIS IS NOT DOING THE RIGHT THING. 
        ### NEED A consice expression for the ratio in the cross- covariance.
        
        
        #thetas=polymult2(covarmat[curr,i]*theta_list[[curr]],deg1=length(theta_list[[curr]])-1,rev(theta_list[[i]]),deg2=0)
        #for(j in (i+1):length(theta_list))
        #phis=
        
        #nums=append(nums, polymult2(theta_list[[curr]],theta_list[curr]        )
        
        
        thetas=polymult2(covarmat[curr,i]*theta_list[[curr]],deg1=length(theta_list[[curr]])-1,rev(theta_list[[i]]),deg2=0)
        phis=polymult2(phi_list[[curr]],deg=length(phi_list[[curr]])-1,rev(phi_list[[i]]),deg2=0)
        
        currratio2=ratioadd(num1=thetas[[1]],num1deg=thetas[[2]],den1=phis[[1]],den1deg=phis[[2]], num2=currratio2[[1]], num2deg=currratio2[[2]], den2=currratio2[[3]], den2deg=currratio2[[4]])
        
        Pz2=currratio2[[1]]
        Pz2Deg=currratio2[[2]]
        Qz2=currratio2[[3]]
        Qz2Deg=currratio2[[4]]
        
        if(length(Pz2)>1){
          while(Pz2[length(Pz2)]==0){
            Pz2=Pz2[c(-length(Pz2))]
            Pz2Deg=Pz2Deg-1
          }
          while(Pz2[1]==0){
            Pz2=Pz2[c(-1)]
          }
        }
        
        if(length(Qz2)>1){
          while(Qz2[length(Qz2)]==0){
            Qz2=Qz2[c(-length(Qz2))]
            Qz2Deg=Qz2Deg-1
          }
          while(Qz2[1]==0){
            Qz2=Qz2[c(-1)]
          }
        }
        
        currratio2=list(Pz2,Pz2Deg,Qz2,Qz2Deg) ## This is for the cross spectra
      }
      
    }
  }
  
  Numadd=polyadd(currratio[[1]],currratio[[2]],currratio2[[1]],currratio2[[2]]) ### UPDATED 4/13/18, accounts for equal denominators in cross and non-cross spectra 
  Pz=Numadd[[1]]
  PzDeg=Numadd[[2]]
  
  if(length(Pz)>1){
    while(Pz[length(Pz)]==0){
      Pz=Pz[c(-length(Pz))]
      PzDeg=PzDeg-1
    }
    while(Pz[1]==0){
      Pz=Pz[c(-1)]
    }
  }
  
  Qz=currratio[[3]]
  Qzdeg=currratio[[4]]
  
  list(Pz=Pz,PzDeg=PzDeg,Qz=Qz,QzDeg=QzDeg) ## If output is Pz=p_{-1}z^{-1}+p_0+p_1z then Pz=c(p_{-1},p_0,p_1)
}



zspec3=function(phi_list, theta_list, covarmat){
  
  ## theta_list contains theta1, theta2, theta3, etc. for as many ARMAs as we are adding together.
  ## phi_list contains phi1, phi2, phi3, etc. for as many ARMAs as we are adding.
  ## These are lists containing vectors of coefficients according to convention:
  ##### the polynomials in the ARMA expression use the convention phi(z)=phi[1]+phi[2]z+phi[3]z^2+...
  ##### Thus if phi(z)=1+2z then phi1=c(1,2)
  ## covarmat is the variance-covariance matrix of the white-noise sequences.
  
  #if(sum(covarmat)==0) browser()
  
  covarmat2=covarmat
  savermv=c()
  for(i in 1:ncol(covarmat)){
    if(covarmat[i,i]<.00000000000000001) {
      savermv=c(savermv,i)
    }
  }
  if(length(savermv)>0){
    phi_list=phi_list[-savermv]
    theta_list=theta_list[-savermv]
    covarmat2=matrix(covarmat2[-savermv,],ncol=ncol(covarmat2))
    covarmat2=matrix(covarmat2[,-savermv],ncol=ncol(covarmat2)-length(savermv))
    covarmat=covarmat2
  }
  
  Pz=list(c(0),c(1))
  for(curr in 1:length(theta_list)){
    for(i in 1:length(theta_list)){
      
      currnum=polymult2(covarmat[curr,i]*theta_list[[curr]],deg1=length(theta_list[[curr]])-1,rev(theta_list[[i]]),deg2=0)
      for(p in (1:length(theta_list))[-curr]){
        currnum=polymult2(currnum[[1]],currnum[[2]],phi_list[[p]],length(phi_list[[p]])-1)
      } 
      for(p in (1:length(theta_list))[-i]){
        currnum=polymult2(currnum[[1]],currnum[[2]],rev(phi_list[[p]]),0)
      }
      
      Pz= polyadd(Pz[[1]],Pz[[2]],currnum[[1]],currnum[[2]])
    }
  }
  PzDeg=Pz[[2]]
  Pz=Pz[[1]]
  
  if(length(Pz)>1){
    while(Pz[length(Pz)]==0){
      Pz=Pz[c(-length(Pz))]
      PzDeg=PzDeg-1
    }
    while(Pz[1]==0){
      Pz=Pz[c(-1)]
    }
  }
  
  Qz=list(c(1),c(0))
  for(curr in 1:length(theta_list)){
    phis=polymult2(phi_list[[curr]],length(phi_list[[curr]])-1,rev(phi_list[[curr]]),0)
    Qz=polymult2(Qz[[1]],Qz[[2]],phis[[1]],phis[[2]])
  }
  
  
  
  QzDeg=Qz[[2]]
  Qz=Qz[[1]]
  
  if(length(Qz)>1){
    while(Qz[length(Qz)]==0){
      Qz=Qz[c(-length(Qz))]
      QzDeg=QzDeg-1
    }
    while(Qz[1]==0){
      Qz=Qz[c(-1)]
    }
  }
  
  list(Pz=Pz,PzDeg=PzDeg,Qz=Qz,QzDeg=QzDeg) ## If output is Pz=p_{-1}z^{-1}+p_0+p_1z then Pz=c(p_{-1},p_0,p_1)
}


############################################################################################################################

half_deg_palindrome=function(coefs=c()){
  ## Function to convert palindromic polynomial in x to smaller degree polynomial in y having "same" roots.
  ## example: (cx^4+bx^3+ax^2+bx+c) OR (x^2+bx+a+bx^(-1)+x^(-2))   ->    y^2+by+a-2
  # first element in coefs is "a", second element in coefs is "b", etc.!
  
  if(length(coefs)==2) coefs= coefs
  
  else{
    for(i in length(coefs):3){
      
      pascal=choose(i-1, 0:(i-1))  ## pascals triangle row which gives binomial coefficients
      pascal=pascal[2:ceiling(length(pascal)/2) ]
      
      for(j in 1:(length(pascal))){
        coefs[i-2*j]=coefs[i-2*j]-pascal[j]*coefs[i] 
      }
    }
  }  
  coefs
}

#### Finding roots of palindromic polynomials. First use the half_deg_palindrome() function!

### Note that for every root y_i of the polynomial in y (whose coefficients are obtained using half_deg_palindrome() function)
### we obtain two roots of the original function by finding the roots of x^2-y_i+1.

roots_of_pal=function(palen_coefs){
  
  coefs=  palen_coefs[ceiling(length(palen_coefs)/2):length(palen_coefs)]
  ycoefs=half_deg_palindrome(coefs)
  roots_y=polyroot(ycoefs)
  
  roots=c()
  for(yi in roots_y){
    roots=c(roots,polyroot(c(1,-yi,1)))
  }
  roots
}




InnVar=function(Pz,PzDeg,Qz,QzDeg){
  
  ### This function computes the innovation variance given a certain z-spectrum stated as a ratio of two Laurent Polynomials.
  ### Pz,PzDeg,Qz,QzDeg will typically be the output of zspec()
  ### In the event of round-off error, this function forces each root z_0 of Pz and Qz inside the unit circle to be paired with root 1/z_0 outside the unit circle.
  
  rootsPz=roots_of_pal(Pz)
  if(is.null(rootsPz)) rootsPz=numeric(0)
  
  rootsQz=roots_of_pal(Qz)
  if(is.null(rootsQz)) rootsQz=numeric(0)
  
  
  rp=Re(Pz[length(Pz)]/((-1)^(PzDeg)*(prod(rootsPz[Mod(rootsPz)<1])))) ## This can have issues when root-finding is imprecise.
  
  rq=Re(Qz[length(Qz)]/((-1)^(QzDeg)*(prod(rootsQz[Mod(rootsQz)<1]))))
  
  list(var=rp/rq,rp=rp,rootsPz=rootsPz,rq=rq,rootsQz=rootsQz) ### added rootsQz to output temporarily.
  
}

############################################################################################################################

getARMAfromInnVar=function(res){
  
  res$rootsPz=res$rootsPz[order(Mod(res$rootsPz))]
  res$rootsQz=res$rootsQz[order(Mod(res$rootsQz))]
  
  agg_theta=list(c(1),0)
  if(length(res$rootsPz)>0){
    for(i in 1:(length(res$rootsPz)/2)){
      agg_theta=polymult2(agg_theta[[1]],agg_theta[[2]],c(1,-res$rootsPz[i]),deg2=1) 
      while(agg_theta[[1]][length(agg_theta[[1]])]==0){ ## this while-loop was moved inside for() loop on 3/27/18
        agg_theta[[1]]=agg_theta[[1]][-length(agg_theta[[1]])]
        agg_theta[[2]]=agg_theta[[2]]-1
      }
    }
    
    
  }
  
  agg_phi=list(c(1),0)
  if(length(res$rootsQz)>0){
    for(i in 1:(length(res$rootsQz)/2)){
      agg_phi=polymult2(agg_phi[[1]],agg_phi[[2]],c(1,-res$rootsQz[i]),deg2=1) 
      while(agg_phi[[1]][length(agg_phi[[1]])]==0){ ## this while-loop was moved inside for() loop on 3/27/18
        agg_phi[[1]]=agg_phi[[1]][-length(agg_phi[[1]])]
        agg_phi[[2]]=agg_phi[[2]]-1
      }
    }
    
    
  }
  
  ############################## outputs agg_phi(z)=1+phi_1z+phi_2z+ ... check this with InnVar!!!!!!!!!!!####################################
  
  agg_phi=polynomial(Re(agg_phi[[-2]])) 
  agg_theta=polynomial(Re(agg_theta[[-2]]))
  
  gcd=GCD(agg_phi,agg_theta)
  
  agg_phi=agg_phi/gcd
  agg_theta=agg_theta/gcd
  
  agg_phi=coef(agg_phi)
  agg_theta=coef(agg_theta)
  
  agg_phi=agg_phi/agg_phi[1]
  agg_theta=agg_theta/agg_theta[1]
  
  list(phi=agg_phi, theta=agg_theta)
  
}


############## This find the Covariance matrix between the shocks in the aggregate sequence:
covaragg=function(agg1phi,agg1theta,agg2phi,agg2theta,phi_list,theta_list,covarmat){
  
  ## phi_list and theta_list should be of size 4 where 1st TWO vectors in each belong to armas in the 1st aggregate
  ## covarmat should correspond to this.
  ## Is only capable of handling subaggregates of TWO ARMAS
  #agg1phi=aggARMAsub1$phi
  #agg1theta=aggARMAsub1$theta
  
  
  library(PolynomF)
  num=polynom(agg1phi)*polynom(theta_list[[1]])
  den=polynom(agg1theta)*polynom(phi_list[[1]])
  MA=coef(num)
  AR=coef(den)
  MAinf1=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf1 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  ## Note that the way ARMAtoMA works is if theta1=1+.5z and phi1=1+.3z, then ar=-.5 and ma=.3
  
  num=polynom(agg2phi)*polynom(theta_list[[3]])
  den=polynom(agg2theta)*polynom(phi_list[[3]])
  MA=coef(num)
  AR=coef(den)
  MAinf2=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf2 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  covarPARTIAL1=sum(MAinf1*MAinf2)*covarmat[1,3]
  
  num=polynom(agg1phi)*polynom(theta_list[[1]])
  den=polynom(agg1theta)*polynom(phi_list[[1]])
  MA=coef(num)
  AR=coef(den)
  MAinf1=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf1 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  num=polynom(agg2phi)*polynom(theta_list[[4]])
  den=polynom(agg2theta)*polynom(phi_list[[4]])
  MA=coef(num)
  AR=coef(den)
  MAinf2=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf2 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  covarPARTIAL2=sum(MAinf1*MAinf2)*covarmat[1,4]
  
  num=polynom(agg1phi)*polynom(theta_list[[2]])
  den=polynom(agg1theta)*polynom(phi_list[[2]])
  MA=coef(num)
  AR=coef(den)
  MAinf1=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf1 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  num=polynom(agg2phi)*polynom(theta_list[[3]])
  den=polynom(agg2theta)*polynom(phi_list[[3]])
  MA=coef(num)
  AR=coef(den)
  MAinf2=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf2 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  covarPARTIAL3=sum(MAinf1*MAinf2)*covarmat[2,3]
  
  num=polynom(agg1phi)*polynom(theta_list[[2]])
  den=polynom(agg1theta)*polynom(phi_list[[2]])
  MA=coef(num)
  AR=coef(den)
  MAinf1=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf1 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  num=polynom(agg2phi)*polynom(theta_list[[4]])
  den=polynom(agg2theta)*polynom(phi_list[[4]])
  MA=coef(num)
  AR=coef(den)
  MAinf2=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500)) ## MAinf2 on document MAILED TO vladimir.kovtun@yu.edu on 4/17/2018. 
  covarPARTIAL4=sum(MAinf1*MAinf2)*covarmat[2,4]
  
  covar=covarPARTIAL1+covarPARTIAL2+covarPARTIAL3+covarPARTIAL4
  covar ### THIS IS THE COVARIANCE BETWEEN THE AGGREGATE SHOCKS. 
}



covaragg2=function(agg1phi,agg1theta,agg2phi,agg2theta,phi_list,theta_list,firstagg= 1:(length(phi_list)/2),secondagg=(length(phi_list)/2+1):length(phi_list),covarmat){
  
  ## phi_list and theta_list contain the original ARMAs
  ## firstagg and secondagg are vectors of "positions" of ARMAs belonging to 1st Aggregate and 2nd Aggregate.
  ## by default these split up in the middle of length of phi_list
  ## covarmat should correspond to this.
  ## Is only capable of handling subaggregates of TWO ARMAS
  
  #agg1phi=aggARMAsub1$phi
  #agg1theta=aggARMAsub1$theta
  #agg2phi=aggARMAsub2$phi
  #agg2theta=aggARMAsub2$theta
  
  library(PolynomF)
  
  covar=0
  for(i in firstagg){
    for(j in secondagg){
      num=polynom(agg1phi)*polynom(theta_list[[i]])
      den=polynom(agg1theta)*polynom(phi_list[[i]])
      MA=coef(num)
      AR=coef(den)
      MAinf1=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500))
      
      num=polynom(agg2phi)*polynom(theta_list[[j]])
      den=polynom(agg2theta)*polynom(phi_list[[j]])
      MA=coef(num)
      AR=coef(den)
      MAinf2=c(1,ARMAtoMA(ar=-AR[-1],ma=MA[-1],lag.max=500))
      
      covar=sum(covar,sum(MAinf1*MAinf2)*covarmat[i,j])
    }
  }
  covar
}



MSFEsubaggTRUE=function(phi_list,theta_list,covarmat,Groups){
  
  if(min(Groups)>1) Groups=Groups-1
  res1_list=list()
  subaggPHI_list=list()
  subaggTHETA_list=list()
  n.clusters=max(Groups)
  nstreams=length(Groups)
  
  for(k in 1:n.clusters){
    zspec= zspec3(phi_list=phi_list[Groups==k],theta_list=theta_list[Groups==k],covarmat=as.matrix(covarmat[Groups==k,Groups==k]))
    res1_list[[k]]=InnVar(zspec$Pz,zspec$PzDeg,zspec$Qz,zspec$QzDeg)
    subaggARMA=getARMAfromInnVar(res1_list[[k]])
    subaggPHI_list[[k]]=Re(subaggARMA$phi)
    subaggTHETA_list[[k]]=Re(subaggARMA$theta)
  }
  
  covarmatsubagg=matrix(rep(NA,n.clusters^2),ncol=n.clusters)
  for(k in 1:n.clusters){
    for(i in (1:n.clusters)[-k]){
      covarmatsubagg[k,i]=round(covaragg2(agg1phi=subaggPHI_list[[k]],agg1theta=subaggTHETA_list[[k]],agg2phi=subaggPHI_list[[i]],agg2theta=subaggTHETA_list[[i]],phi_list,theta_list,firstagg=(1:nstreams)[Groups==k],secondagg=(1:nstreams)[Groups==i],covarmat=covarmat),14)
    }
  }
  
  for(k in 1:n.clusters){
    covarmatsubagg[k,k]=res1_list[[k]]$var
  }
  
  sum(covarmatsubagg)
  
}





objfun=function(thetacoef=c(), clusters=c()){
  ## computes value of objective function for given MA(1) coefficients of each stream and clusters
  ## clusters assigns streams to clusters
  
  obj=0
  for(cl in unique(clusters)){
    coefs=thetacoef[clusters==cl]
    b_k=sum(coefs^2)+sum(clusters==cl)
    a_k=sum(coefs)
    obj_clust= sqrt((b_k+2*a_k)*(b_k-2*a_k))
    obj=obj+obj_clust 
  }
  obj
}


half_deg_palindrome=function(coefs=c()){
  ## Function to convert palindromic polynomial in x to smaller degree polynomial in y having "same" roots.
  ## example: (x^4+bx^3+ax^2+bx+1) OR (x^2+bx+a+bx^(-1)+x^(-2))   ->    y^2+by+a-2
  # first element in coefs is "a", second element in coefs is "b", last element in coefs should be 1!
  
  for(i in length(coefs):3){
    
    pascal=choose(i-1, 0:(i-1))  ## pascals triangle row which gives binomial coefficients
    pascal=pascal[2:ceiling(length(pascal)/2) ]
    
    for(j in 1:(length(pascal))){
      coefs[i-2*j]=coefs[i-2*j]-pascal[j]*coefs[i] 
    }
  }
  coefs 
}


#coefs=c(2,4,8,5,10,1)
#half_deg_palindrome(coefs)

#### Finding roots of palindromic polynomials. First use the half_deg_palindrome() function!

### Note that for every root y_i of the polynomial in y (whose coefficients are obtained using half_deg_palindrome() function)
### we obtain two roots of the original function by finding the roots of x^2-y_i+1.

roots_of_pal=function(palen_coefs){
  
  if(length(palen_coefs)<=3) roots=polyroot(palen_coefs)
  
  else{
    coefs=  palen_coefs[ceiling(length(palen_coefs)/2):length(palen_coefs)]
    ycoefs=half_deg_palindrome(coefs)
    roots_y=polyroot(ycoefs)
    
    roots=c()
    for(yi in roots_y){
      roots=c(roots,polyroot(c(1,-yi,1)))
    }
  }
  roots
}




objfun_y_arma=function(yvec){
  ###### ALL CLUTERS WILL HAVE ALL STREAMS IN THEM, BUT USING FRACTIONAL y_i's
  ### Assumed leadtime of ZERO.
  
  ymat=matrix(yvec,ncol=nstreams)
  aggARMA_phi_list=list()
  aggARMA_theta_list=list()
  
  psi_star_list=list()
  omega_star_list=list()
  
  psi_tilde_list=list()
  
  if( sum(rowSums(ymat)^2< .000000000000000001)>0) obj=1000000000 ## This happens when less clusters then specified are identified.
  
  else{
    for(k in 1:n.clusters){
      covarmat_y=covarmat 
      for(i in 1:nrow(covarmat)){
        for(j in 1: ncol(covarmat)){
          if(i==j) covarmat_y[i,j]=covarmat[i,j]*ymat[k,i]
          else covarmat_y[i,j]=covarmat[i,j]*ymat[k,i]*ymat[k,j]
        }
      }
      
      
      zspec=zspec3(phi_list=phi_list,theta_list=theta_list,covarmat=covarmat_y) ### zspec3
      res=InnVar(zspec$Pz,zspec$PzDeg,zspec$Qz,zspec$QzDeg) ## roots of Pz and Qz unstable (complex root conjugates are not same)
      aggARMA=getARMAfromInnVar(res) ## arma representation of each cluster using all streams and fractional y_i (30)
      for(i in 1:length(aggARMA)) aggARMA[[i]]=Re(aggARMA[[i]])
      aggARMA_phi_list=append(aggARMA_phi_list,list(aggARMA$phi))
      aggARMA_theta_list=append(aggARMA_theta_list,list(aggARMA$theta))
      
      
      ## create psi_stars (only needed for leadtime>0):
      #psi_star=c(1,ARMAtoMA(ar=aggARMA_phi_list[[k]][-1], ma=aggARMA_theta_list[[k]][-1],MAinfLAG))
      #psi_star_list[[k]]=psi_star
      #omega_star=psi_star ## assumes a leadtime of 0.
      #omega_star_list[[k]]=omega_star
      
      ### create psi_tilde:
      psi_tilde_eachclust_list=list()
      for(i in 1:nstreams){
        
        num=polynomial(aggARMA_phi_list[[k]])*polynomial(theta_list[[i]])
        num=num[1:length(num)]
        den=polynomial(aggARMA_theta_list[[k]])*polynomial(phi_list[[i]])
        den=den[1:length(den)]
        
        #####################################################################################
        ### The following gets rid of common roots between numerator and denominator in (26):
        numroots=polyroot(num)
        denroots=polyroot(den)
        
        sameroots=as.character(numroots)[as.character(numroots)%in%as.character(denroots)]
        remvnum=c()
        remvden=c()
        for(root in sameroots){
          remvnum=c(remvnum,(1:length(numroots))[(as.character(numroots)==root)][1])
          remvden=c(remvden,(1:length(denroots))[(as.character(denroots)==root)][1])
        }
        
        if(length(remvnum)>0) numroots=numroots[-remvnum]
        if(length(remvden)>0) denroots=denroots[-remvden]
        
        numrootsRE=Re(numroots[(Im(numroots))^2<0.000000001])
        numrootsIM=numroots[(Im(numroots))^2>0.000000001]
        numrootsIM=numrootsIM[sort.list(Mod(numrootsIM))]
        
        denrootsRE=Re(denroots[(Im(denroots))^2<0.000000001])
        denrootsIM=denroots[(Im(denroots))^2>0.000000001]
        denrootsIM=denrootsIM[sort.list(Mod(denrootsIM))]
        
        if(length(numroots)>0){
          num=polynomial(c(1))
          if(length(numrootsRE)>0){ for(p in 1:length(numrootsRE)) num=num* polynomial(c(1, -1/numrootsRE[p]))}
          if(length(numrootsIM)>0){ for(p in 2*(1:(length(numrootsIM)/2))) num=num* polynomial( c(1,Re(-1/numrootsIM[p]-1/numrootsIM[p-1]), Re(1/(numrootsIM[p]*numrootsIM[p-1]))))}
          num=num[1:length(num)]
        }
        else num=1
        
        if(length(denroots)>0){
          den=polynomial(c(1))
          if(length(denrootsRE)>0){ for(p in 1:length(denrootsRE)) den=den* polynomial(c(1, -1/denrootsRE[p]))}
          if(length(denrootsIM)>0){ for(p in 2*(1:(length(denrootsIM)/2))) den=den* polynomial( c(1,Re(-1/denrootsIM[p]-1/denrootsIM[p-1]), Re(1/(denrootsIM[p]*denrootsIM[p-1]))))}
          den=den[1:length(den)]
        }
        else den=1
        ### common root elimination stops here.
        ############################################################################################
        
        psi_tilde_eachclust_list= append(psi_tilde_eachclust_list,list(c(1,ARMAtoMA(ar=-den[-1],ma=num[-1],MAinfLAG))))
        
      }
      
      psi_tilde_list[[k]]=psi_tilde_eachclust_list
      
    }
    
    ## The following computes equation (28):
    
    obj=0
    for(alpha in 1:n.clusters){
      for(beta in 1:n.clusters){
        
        #for leadtime other than 1, we need an extra loop with  omega=omega_star_list[[alpha]][leadtime+1]*omega_star_list[[beta]][leadtime+1]
        ## Note that omega_{beta,0}=1 and that is all we need for leadtime of zero
        omega=1
        
        for(i in 1:nstreams){
          for(j in 1:nstreams){
            for(k in 1:MAinfLAG){
              obj=obj+omega*ymat[alpha,i]*ymat[beta,j]*covarmat[i,j]*psi_tilde_list[[alpha]][[i]][k]*psi_tilde_list[[beta]][[j]][k]
            }
          }
        }
        
        
      }
    }
  }
  obj
  
}  





yvec_from_Groups=function(Groups,n.clusters,nstreams){
  
  ymat=c()
  for(i in 1:length(Groups)){
    add=rep(0,n.clusters)
    add[Groups[i]]=1
    ymat=cbind(ymat,add)
  }
  yvec=c(ymat)
  yvec
}


pivotMSFE_TRUE=function(phi_list,theta_list,covarmat,MSFEcurr,Groups,unfin.clust=1:max(Groups)){
  
  n.groups=max(Groups)
  stream=1
  
  while(length(unfin.clust)>0){
    
    if(stream == 1) clust=unfin.clust[1] ## this makes sure we only consider clusters not yet exhausted
    # changed from clust=sample(unfin.clust,1) on 11/13/18, after INFORMS
    newclusters=(1:n.groups)[-clust]
    MSFEtrack=c()
    
    if(sum(Groups==clust)>1){ ## makes sure that we do not try reducing number of clusters.
      for(j in 1:length(newclusters)){
        newclust=newclusters[j]
        GroupsPiv=Groups
        GroupsPiv[Groups==clust][stream]=newclust ## now clusters have been updated
        
        MSFEtrack=rbind(MSFEtrack,c(newclusters[j],MSFEsubaggTRUE(phi_list,theta_list,covarmat,GroupsPiv)))
      } 
      
      MSFE=MSFEtrack[MSFEtrack[,2]==min(MSFEtrack[,2]),2]
      newclust=MSFEtrack[MSFEtrack[,2]==min(MSFEtrack[,2]),1]
      GroupsPiv=Groups
      GroupsPiv[Groups==clust][stream]=newclust
    }
    else MSFE = MSFEcurr
    if(MSFE<MSFEcurr & stream<sum(Groups==clust)){
      stream=1
      Groups=GroupsPiv
      MSFEcurr=MSFE
    } 
    
    else if(MSFE<MSFEcurr & stream==sum(Groups==clust)){
      unfin.clust=unfin.clust[unfin.clust!=clust]
      Groups=GroupsPiv
      MSFEcurr=MSFE
      stream=1
    } 
    
    else if(MSFE>=MSFEcurr & stream<sum(Groups==clust)){
      stream=stream+1
    } 
    
    else if(MSFE>=MSFEcurr & stream==sum(Groups==clust)){
      unfin.clust=unfin.clust[unfin.clust!=clust]
      stream=1
    } 
    
  }
  
  retlist=list(Groups, MSFEcurr)
  retlist
  
}




pivotMSFE2=function(ProcessDataPiv,MSFEcurr, Groups,unfin.clust=1:max(Groups),ar.order=10,ma.order=10){ ## version 2 with loop.
  ##### select a random cluster and select a random stream in the cluster:
  ## clust represents the cluster from which the stream will be selected.
  ### Seems to converge to local minimum after 1 or 2 iterations.
  ## ar.order and ma.order specify the order of the ARMA models estimated for the subaggregated clusters
  
  Nobs=nrow(ProcessDataPiv)
  n.groups=max(Groups)
  stream=1
  
  while(length(unfin.clust)>0){
    
    if(stream == 1) clust=unfin.clust[1] ## this makes sure we only consider clusters not yet exhausted
    # changed from clust=sample(unfin.clust,1) on 11/13/18, after INFORMS
    newclusters=(1:n.groups)[-clust]
    MSFEtrack=c()
    for(j in 1:length(newclusters)){
      newclust=newclusters[j]
      GroupsPiv=Groups
      GroupsPiv[Groups==clust][stream]=newclust ## now clusters have been updated
      ClData=t(rowsum(t(ProcessDataPiv),GroupsPiv))
      
      res_clusters=c()
      order=c(ar.order,0,ma.order) ### below arma(ClData[,i],order=order)$residuals[(max(1,order)+1):Nobs] was changed to arima(ClData[,i],order=order,method="CSS")$residuals 12/13/18
      for(i in 1:ncol(ClData)) res_clusters=cbind(res_clusters,arima(ClData[,i],order=order,method="CSS")$residuals)
      MSFEtrack=rbind(MSFEtrack,c(newclusters[j],sum(cov(as.matrix(res_clusters)))))
      
    }
    
    MSFE=MSFEtrack[MSFEtrack[,2]==min(MSFEtrack[,2]),2]
    
    newclust=MSFEtrack[MSFEtrack[,2]==min(MSFEtrack[,2]),1]
    GroupsPiv=Groups
    GroupsPiv[Groups==clust][stream]=newclust
    
    
    if(MSFE<MSFEcurr & stream<sum(Groups==clust)){
      stream=1
      Groups=GroupsPiv
      MSFEcurr=MSFE
    } 
    
    else if(MSFE<MSFEcurr & stream==sum(Groups==clust)){
      unfin.clust=unfin.clust[unfin.clust!=clust]
      Groups=GroupsPiv
      MSFEcurr=MSFE
      stream=1
    } 
    
    else if(MSFE>=MSFEcurr & stream<sum(Groups==clust)){
      stream=stream+1
    } 
    
    else if(MSFE>=MSFEcurr & stream==sum(Groups==clust)){
      unfin.clust=unfin.clust[unfin.clust!=clust]
      stream=1
    } 
    
  }
  
  retlist=list(Groups, MSFEcurr)
  retlist
}





