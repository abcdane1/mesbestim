#set seed 
RNGkind("L'Ecuyer-CMRG")
set.seed(01051990)

data_aspire_fun_dfme<-function(m,csl,csu,sdvalout,sdvalsel,sdvalme,diffx=T,wdoc=F,
                                          sdvaldout,sdvaldsel,sdvaldme,
                                          dlb=2,dub=3){
  
  #cs<-sample(50:100,m,replace=T)
  cs<-sample(csl:csu,m,replace=T)
  #cs<-rep(csl,m)
  #treatment assignment
  trt<-rbinom(m,1,p=.5)
  
  #individual information
  x1<-rep(0,sum(cs))
  x2<-rep(0,sum(cs))
  x3<-rep(0,sum(cs))
  x4<-rep(0,sum(cs))
  
  y<-rep(0,sum(cs))
  ys<-rep(0,sum(cs))
  v<-rep(0,sum(cs))
  y0<-rep(0,sum(cs))
  y1<-rep(0,sum(cs))
  
  nds<-sample(dlb:dub,m,replace=T) #2 to 3
  if(dlb==dub){
    nds<-rep(dlb,m)
  }
  
  doc<-rep(0,sum(cs))
  
  d0<-1
  d1<-1
  for(i in 1:m){
    if(nds[i]==1){
      doc[d1:(d1+cs[i]-1)]<-d0
    }else{
      doc[d1:(d1+cs[i]-1)]<-sort(sample(d0:(d0+nds[i]-1),cs[i],replace=T))}
    d0<-d0+nds[i]
    d1<-d1+cs[i]
  }
  
  
  #induce within group correlation (like in glmm)
  re_function<-function(x,n){
    if(x==0){
      return(rep(0,n))
    }else{
      return(rnorm(n=n,mean=0,sd=x))
    }
  }
  
  #random effects of cluster
  bi<-re_function(sdvalout,m) #cluster outcome
  bis<-re_function(sdvalsel,m) #selection
  bime<-re_function(sdvalme,m) #measurement error
  
  #ensure entries are 0
  if(wdoc==F){
    sdvaldout<-0
    sdvaldsel<-0
    sdvaldme<-0}
  
  #doctor level random effects 
  bid<-re_function(sdvaldout,sum(nds))
  bids<-re_function(sdvaldsel,sum(nds))
  bidme<-re_function(sdvaldme,sum(nds))
  
  #data generation, adding cluster period
  c<-1
  for(i in 1:m){
    ##cluster information
    csi<-cs[i]
    clustidr<-c:(c+csi-1)
    #random effects
    bii<-bi[i]
    bisi<-bis[i]
    bimei<-bime[i]
    #treatment 
    trti<-trt[i]
    #doctor id
    bidi<-bid[doc[clustidr]]
    bidsi<-bids[doc[clustidr]]
    bidmei<-bidme[doc[clustidr]]
    
    #individual covariates
    mean1<-rep(1,csi)
    sigmacov<-diag(csi)
    valx1<-MASS::mvrnorm(n=1,mu=mean1,Sigma=sigmacov)
    x1[clustidr]<-valx1
    
    mean2<-rep(.5,csi)
    sigmacov2<-.5*diag(csi)+.05
    valx2<-MASS::mvrnorm(n=1,mu=mean2,Sigma=sigmacov2)
    x2[clustidr]<-valx2
    
    valx3<-rbinom(n=csi,1,p=.55)
    x3[clustidr]<-valx3
    
    valx4<-runif(1)
    x4[clustidr]<-valx4
    
    ##y1 and y0 generation
    #p0<-1/(1+exp(.75-.05*valx1-.1*valx2-.05*valx3+.05*valx4-bii-bidi))
    #p1<-1/(1+exp(.75-.75-.05*valx1-.1*valx2-.05*valx3+.05*valx4-bii-bidi))
    
    p0<-1/(1+exp(1-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
    p1<-1/(1+exp(1-.75-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
    
    y0val<-rbinom(n=csi,size=1,p0)
    y1val<-rbinom(n=csi,size=1,p1)
    
    y0[clustidr]<-y0val
    y1[clustidr]<-y1val
    
    ##me error generation
    #if doesn't depend on x
    if(diffx==F){
      # ps1<-1/(1+exp(1-2.5*y1val+.1-bimei))
      # ps0<-1/(1+exp(1-2.5*y0val-bimei))
      ps1<-1/(1+exp(1-2.5*y1val-bimei-bidmei))
      ps0<-1/(1+exp(1.25-1.5*y0val-bimei-bidmei))
    }
    
    #if does depend on x
    if(diffx==T){
      ps1<-1/(1+exp(.75-2.5*y1val+.25*valx1+.15*valx2+.25*valx3-.1*valx4-bimei-bidmei)) #re version
      ps0<-1/(1+exp(1.25-1.5*y0val-.25*valx1+.25*valx2+.15*valx3-.1*valx4-bimei-bidmei))
    }
    
    #measurement error probabilities
    ps<-ifelse(rep(trti,csi)==1,ps1,ps0)
    ysval<-rbinom(n=csi,size=1,p=ps)
    ys[clustidr]<-ysval
    
    ##selection and outcome models
    
    #selection probability
    pselect1<-1/(1+exp(.5-.15*y1val+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi))#- for now
    pselect0<-1/(1+exp(.25+.15*y0val+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi))
    
    # pselect1<-1/(1+exp(.5+.1*y1val+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi))#- for now
    # pselect0<-1/(1+exp(.25+.3*y0val+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi))
    
    pselect<-ifelse(rep(trti,csi)==1,pselect1,pselect0)
    #y values by treatment 
    yval<-ifelse(rep(trti,csi)==1,y1val,y0val)
    
    #selection indiciators
    selectind<-rbinom(n=csi,size=1,pselect)
    v[clustidr]<-selectind
    #observed ys
    y[clustidr]<-ifelse(selectind==1,yval,NA)
    
    c<-c+csi}
  df<-data.frame(y=y,ys=ys,Id=rep(1:m,cs),Doc=doc,a=rep(trt,cs),x1=x1,x2=x2,x3=x3,x4=x4,v=v,y0=y0,y1=y1)
  return(df)}

# data_aspire_fun_dfme<-function(m,csl,csu,sdvalout,sdvalsel,sdvalme,diffx=T,wdoc=F,
#                                sdvaldout,sdvaldsel,sdvaldme,
#                                dlb=2,dub=3){
#   
#   #cs<-sample(50:100,m,replace=T)
#   cs<-sample(csl:csu,m,replace=T)
#   #cs<-rep(csl,m)
#   #treatment assignment
#   trt<-rbinom(m,1,p=.5)
#   
#   #individual information
#   x1<-rep(0,sum(cs))
#   x2<-rep(0,sum(cs))
#   x3<-rep(0,sum(cs))
#   x4<-rep(0,sum(cs))
#   
#   y<-rep(0,sum(cs))
#   ys<-rep(0,sum(cs))
#   v<-rep(0,sum(cs))
#   y0<-rep(0,sum(cs))
#   y1<-rep(0,sum(cs))
#   
#   nds<-sample(dlb:dub,m,replace=T) #2 to 3
#   if(dlb==dub){
#     nds<-rep(dlb,m)
#   }
#   
#   doc<-rep(0,sum(cs))
#   
#   d0<-1
#   d1<-1
#   for(i in 1:m){
#     if(nds[i]==1){
#       doc[d1:(d1+cs[i]-1)]<-d0
#     }else{
#       doc[d1:(d1+cs[i]-1)]<-sort(sample(d0:(d0+nds[i]-1),cs[i],replace=T))}
#     d0<-d0+nds[i]
#     d1<-d1+cs[i]
#   }
#   
#   
#   #induce within group correlation (like in glmm)
#   re_function<-function(x,n){
#     if(x==0){
#       return(rep(0,n))
#     }else{
#       return(rnorm(n=n,mean=0,sd=x))
#     }
#   }
#   
#   #random effects of cluster
#   bi<-re_function(sdvalout,m) #cluster outcome
#   bis<-re_function(sdvalsel,m) #selection
#   bime<-re_function(sdvalme,m) #measurement error
#   
#   #ensure entries are 0
#   if(wdoc==F){
#     sdvaldout<-0
#     sdvaldsel<-0
#     sdvaldme<-0}
#   
#   #doctor level random effects 
#   bid<-re_function(sdvaldout,sum(nds))
#   bids<-re_function(sdvaldsel,sum(nds))
#   bidme<-re_function(sdvaldme,sum(nds))
#   
#   #data generation, adding cluster period
#   c<-1
#   for(i in 1:m){
#     ##cluster information
#     csi<-cs[i]
#     clustidr<-c:(c+csi-1)
#     #random effects
#     bii<-bi[i]
#     bisi<-bis[i]
#     bimei<-bime[i]
#     #treatment 
#     trti<-trt[i]
#     #doctor id
#     bidi<-bid[doc[clustidr]]
#     bidsi<-bids[doc[clustidr]]
#     bidmei<-bidme[doc[clustidr]]
#     
#     #individual covariates
#     mean1<-rep(1,csi)
#     sigmacov<-diag(csi)
#     valx1<-MASS::mvrnorm(n=1,mu=mean1,Sigma=sigmacov)
#     x1[clustidr]<-valx1
#     
#     mean2<-rep(.5,csi)
#     sigmacov2<-.5*diag(csi)+.05
#     valx2<-MASS::mvrnorm(n=1,mu=mean2,Sigma=sigmacov2)
#     x2[clustidr]<-valx2
#     
#     valx3<-rbinom(n=csi,1,p=.55)
#     x3[clustidr]<-valx3
#     
#     valx4<-runif(1)
#     x4[clustidr]<-valx4
#     
#     ##y1 and y0 generation
#     #p0<-1/(1+exp(.75-.05*valx1-.1*valx2-.05*valx3+.05*valx4-bii-bidi))
#     #p1<-1/(1+exp(.75-.75-.05*valx1-.1*valx2-.05*valx3+.05*valx4-bii-bidi))
#     
#     p0<-1/(1+exp(1-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
#     p1<-1/(1+exp(1-.75-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
#     
#     y0val<-rbinom(n=csi,size=1,p0)
#     y1val<-rbinom(n=csi,size=1,p1)
#     
#     y0[clustidr]<-y0val
#     y1[clustidr]<-y1val
#     
#     ##me error generation
#     #if doesn't depend on x
#     if(diffx==F){
#       # ps1<-1/(1+exp(1-2.5*y1val+.1-bimei))
#       # ps0<-1/(1+exp(1-2.5*y0val-bimei))
#       ps1<-1/(1+exp(1-2.5*y1val-bimei-bidmei))
#       ps0<-1/(1+exp(1.25-1.75*y0val-bimei-bidmei))
#     }
#     
#     #if does depend on x
#     if(diffx==T){
#       ps1<-1/(1+exp(.75-2.5*y1val+.25*valx1+.15*valx2+.25*valx3-.1*valx4-bimei-bidmei)) #re version
#       ps0<-1/(1+exp(1.25-1.75*y0val-.25*valx1+.25*valx2+.15*valx3-.1*valx4-bimei-bidmei))
#     }
#     
#     #measurement error probabilities
#     ps<-ifelse(rep(trti,csi)==1,ps1,ps0)
#     ysval<-rbinom(n=csi,size=1,p=ps)
#     ys[clustidr]<-ysval
#     
#     ##selection and outcome models
#     
#     #selection probability
#     pselect<-1/(1+exp(.5-.25*trti+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi)) #- for now
#     #y values by treatment 
#     yval<-ifelse(rep(trti,csi)==1,y1val,y0val)
#     
#     #selection indiciators
#     selectind<-rbinom(n=csi,size=1,pselect)
#     v[clustidr]<-selectind
#     #observed ys
#     y[clustidr]<-ifelse(selectind==1,yval,NA)
#     
#     c<-c+csi}
#   df<-data.frame(y=y,ys=ys,Id=rep(1:m,cs),Doc=doc,a=rep(trt,cs),x1=x1,x2=x2,x3=x3,x4=x4,v=v,y0=y0,y1=y1)
#   return(df)}


# 
# #empirical selection and mis-measurement proportions
feature_fun_plus<-function(df){
  #selection proportion
  select_prob<-mean(df$v)
  dfv<-df[which(df$v==1),]
  #overall mis-classification proportion
  overall_missclass<-mean(dfv$y!=dfv$ys)
  #treatment 1
  dfv1<-dfv[which(dfv$a==1),]
  #proportion of y=1 and y*=0
  trt1y1ys0_missclass<-mean(dfv1$y==1 & dfv1$ys==0)
  #proportion of y=0 and y*=1
  trt1y0ys1_missclass<-mean(dfv1$y==0 & dfv1$ys==1)
  #treatment 0
  dfv0<-dfv[which(dfv$a==0),]
  #proportion of y=1 and y*=0
  trt0y1ys0_missclass<-mean(dfv0$y==1 & dfv0$ys==0)
  #proportion of y=0 and y*=1
  trt0y0ys1_missclass<-mean(dfv0$y==0 & dfv0$ys==1)
  output<-c(select_prob,overall_missclass,trt1y1ys0_missclass,trt1y0ys1_missclass,trt0y1ys0_missclass,trt0y0ys1_missclass)
  names(output)<-c("SelectProp","OverallMissclass","Trt1Missclass_True1Meas0","Trt1Missclass_True0Meas1",
                   "Trt0Missclass_True1Meas0","Trt0Missclass_True0Meas1")
  return(output)
}

estimator_aspire_fun_v3<-function(df,xlabs,clabs,varprint=T,diffx=T,cor=T,ind=T){
  N<-nrow(df)
  dfv<-df[which(df$v==1),]
  
  new111<-data.frame(df[,c(xlabs,clabs)],a=1,y=1)
  new101<-data.frame(df[,c(xlabs,clabs)],a=1,y=0)
  new110<-data.frame(df[,c(xlabs,clabs)],a=0,y=1)
  new100<-data.frame(df[,c(xlabs,clabs)],a=0,y=0)
  
  if(diffx==T){
    formdiffx<-paste("ys","~",paste("a","y",paste(c(xlabs,clabs),
                                                  collapse="+"),paste(paste("a",c(xlabs,"y"),sep="*"),collapse="+"),sep="+"),sep="")
    formdiffx<-as.formula(formdiffx)
    mefit<-glm(formdiffx,family="binomial",data=dfv)
    
    dfvcov<-cbind(int=1,dfv[,c("a","y",xlabs,clabs)],dfv$a*dfv[,c(xlabs,"y")])
    df111cov<-cbind(int=1,new111[,c("a","y",xlabs,clabs)],new111[,c(xlabs,"y")])
    df101cov<-cbind(int=1,new101[,c("a","y",xlabs,clabs)],new101[,c(xlabs,"y")])
    df110cov<-cbind(int=1,new110[,c("a","y",xlabs,clabs)],0*new110[,c(xlabs,"y")])
    df100cov<-cbind(int=1,new100[,c("a","y",xlabs,clabs)],0*new100[,c(xlabs,"y")])
  }
  
  if(diffx==F){
    formsamex<-paste("ys","~",paste("a","y",paste(paste("a","y",sep="*"),collapse="+"),sep="+"),sep="")
    formsamex<-as.formula(formsamex)
    mefit<-glm(formsamex,family="binomial",data=dfv)
    
    dfvcov<-cbind(int=1,dfv[,c("a","y")],dfv$a*dfv[,"y"])
    df111cov<-cbind(int=1,new111[,c("a","y")],new111[,"y"])
    df101cov<-cbind(int=1,new101[,c("a","y")],new101[,"y"])
    df110cov<-cbind(int=1,new110[,c("a","y")],0*new110[,"y"])
    df100cov<-cbind(int=1,new100[,c("a","y")],0*new100[,"y"])
    
  }
  
  p111hat<-predict(mefit,newdata=new111,type="response")
  p101hat<-predict(mefit,newdata=new101,type="response")
  p110hat<-predict(mefit,newdata=new110,type="response")
  p100hat<-predict(mefit,newdata=new100,type="response")
  phat<-mean(df$a)
  
  u1hat<-1/N*sum((df$ys*df$a-phat*p101hat)/(phat*(p111hat-p101hat)))
  u0hat<-1/N*sum((df$ys*(1-df$a)-(1-phat)*p100hat)/((1-phat)*(p110hat-p100hat)))
  
  finalest<-u1hat-u0hat
  truth<-mean(df$y1-df$y0)
  
  if(varprint==F){
    return(c("Estimator"=finalest))
  }else{
    
    pmehat<-fitted.values(mefit)
    nparams<-ncol(dfvcov)
    
    # meatmat<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    # breadmat<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    meatmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    breadmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    m<-length(unique(df$Id))
    
    for(i in 1:m){
      clustid<-which(dfv$Id==i)
      dfviclustcov<-dfvcov[clustid,]
      pmehati<-pmehat[clustid]
      ysi<-dfv[clustid,]$ys
      estfuni<-t(t(dfviclustcov)%*%(ysi-pmehati)) #this is the estimating function by group, the phi_i, summing over j within the group
      
      clustidfull<-which(df$Id==i)
      dfi<-df[clustidfull,]
  
      
      mi5x<-sum(dfi$a-phat)
      #mi6x<-sum(dfi$ys*dfi$a-u1hat*phat*(p111hat[clustidfull]-p101hat[clustidfull])-phat*p101hat[clustidfull])
      mi6x<-sum((dfi$ys*dfi$a-phat*p101hat[clustidfull])/(phat*(p111hat[clustidfull]-p101hat[clustidfull]))-u1hat)
      #mi7x<-sum(dfi$ys*(1-dfi$a)-u0hat*(1-phat)*(p110hat[clustidfull]-p100hat[clustidfull])-(1-phat)*p100hat[clustidfull])
      mi7x<-sum((dfi$ys*(1-dfi$a)-(1-phat)*p100hat[clustidfull])/((1-phat)*(p110hat[clustidfull]-p100hat[clustidfull]))-u0hat)
      # meatmat<-meatmat+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
      meatmat1<-meatmat1+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
      
    }
    
    breadmat1[1:nparams,1:nparams]<--t(dfvcov*pmehat*(1-pmehat))%*%as.matrix(dfvcov)
    
    breadmat1[nparams+1,nparams+1]<--N
    
    breadmat1[nparams+2,1:nparams]<-apply((df101cov*phat^2*(-p101hat)*(1-p101hat)*(p111hat-p101hat)-(df$a*df$ys-phat*p101hat)*phat*(df111cov*p111hat*(1-p111hat)-df101cov*p101hat*(1-p101hat)))/
                                            (phat^2*(p111hat-p101hat)^2),2,sum)
    
    
    breadmat1[nparams+2,nparams+1]<-sum((-p101hat*(phat)-(df$a*df$ys-phat*p101hat))/(phat^2*(p111hat-p101hat)))
    
    breadmat1[nparams+2,nparams+2]<--N
    
    breadmat1[nparams+3,1:nparams]<-apply((df100cov*(1-phat)^2*(-p100hat)*(1-p100hat)*(p110hat-p100hat)-((1-df$a)*df$ys-(1-phat)*p100hat)*(1-phat)*(df110cov*p110hat*(1-p110hat)-df100cov*p100hat*(1-p100hat)))/
                                            ((1-phat)^2*(p110hat-p100hat)^2),2,sum)
  
    breadmat1[nparams+3,nparams+1]<-sum((p100hat*((1-phat))+((1-df$a)*df$ys-(1-phat)*p100hat))/((1-phat)^2*(p110hat-p100hat)))
    
    breadmat1[nparams+3,nparams+3]<--N
    
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat1%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)
    
    lb<-finalest+qnorm(.025)*sqrt(varest)
    ub<-finalest+qnorm(.975)*sqrt(varest)
    
    names_og<-c("Truth","Estimator","VarAsymCL","LBCL","UBCL")
    
    if(cor==F){
      output<-c(truth,finalest,varest,lb,ub)
      names(output)<-names_og
    }
    
    if(cor==T){
      varest1<-m/(m-7)*varest
      lb1Df<-finalest+qnorm(.025)*sqrt(varest1)
      ub1Df<-finalest+qnorm(.975)*sqrt(varest1)
      
      lb1T<-finalest+qt(.025,df=m-7)*sqrt(varest)
      ub1T<-finalest+qt(.975,df=m-7)*sqrt(varest)
      
      output<-c(truth,finalest,varest,lb,ub,varest1,lb1Df,ub1Df,varest,lb1T,ub1T)
      names(output)<-c(names_og,"VarAsymCLCorDf","LBCLCorDf","UBCLCorDf",
                       "VarAsymCLCorT","LBCLCorT","UBCLCorT")
    }
    
    if(ind==T){
      meatmat2<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
      #breadmat2<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
      for(i in 1:N){
        dfi<-df[i,]
        if(dfi$v==1){
          idv<-which(row.names(dfv)==i)
          pmehati<-pmehat[idv]
          dfvicov<-dfvcov[idv,]
          ysi<-dfv[idv,]$ys
          estfuni<-t(t(dfvicov)%*%(ysi-pmehati))}else{
            estfuni<-rep(0,nparams)}
        #this is the estimating function by group, the phi_i, summing over j within the group
        #mi5x<-sum(dfviint$a-phat)
        #mi6x<-sum(dfviint$ys*dfviint$a-u1hat*phat*(p111hat[clustid]-p101hat[clustid])-phat*p101hat[clustid])
        #mi7x<-sum(dfviint$ys*(1-dfviint$a)-u0hat*(1-phat)*(p110hat[clustid]-p100hat[clustid])-(1-phat)*p100hat[clustid])
        #meatmat1<-meatmat1+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
        
        mi5x<-sum(dfi$a-phat)
        #mi6x<-sum(dfi$ys*dfi$a-u1hat*phat*(p111hat[i]-p101hat[i])-phat*p101hat[i])
        #mi7x<-sum(dfi$ys*(1-dfi$a)-u0hat*(1-phat)*(p110hat[i]-p100hat[i])-(1-phat)*p100hat[i])
        mi6x<-sum((dfi$ys*dfi$a-phat*p101hat[i])/(phat*(p111hat[i]-p101hat[i]))-u1hat)
        mi7x<-sum((dfi$ys*(1-dfi$a)-(1-phat)*p100hat[i])/((1-phat)*(p110hat[i]-p100hat[i]))-u0hat)
        meatmat2<-meatmat2+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
      }
      
      varest2<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat2%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)
      lb2<-finalest+qnorm(.025)*sqrt(varest2)
      ub2<-finalest+qnorm(.975)*sqrt(varest2)
      output<-c(output,varest2,lb2,ub2)
      fin<-length(output)
      names(output)[(fin-2):fin]<-c("VarAsymInd","LBInd","UBInd")
    }
    return(output)}
}

#bootstrap
boot_aspire_v3<-function(iters=1000,xlabs,clabs,clust=T,diffx=T,dfog){
  estboot<-rep(0,iters)
  for(b in 1:iters){
    if(clust==T){
      idval<-length(unique(dfog$Id))
      clustidv<-sample(1:idval,idval,replace=T)
      dflist<-vector("list",length(idval))
      for(j in 1:idval){
        dflist[[j]]<-dfog[dfog$Id==clustidv[j],]
      }
      df<-do.call(rbind,dflist)
    }
    if(clust==F){
      rowind<-sample(1:nrow(dfog),nrow(dfog),replace=T)
      df<-dfog[rowind,]
    }
    
    estim_boot<-estimator_aspire_fun_v3(df,xlabs=xlabs,clabs=clabs,varprint=F,diffx=diffx) 
    #avoid degenerate fits
    if(!is.na(estim_boot) && abs(estim_boot)<=1){
      estboot[b]<-estim_boot
    }else{
      estboot[b]<-NA
    }}
  
  varboot<-var(estboot,na.rm=T) #for now 
  bounds<-quantile(estboot,c(.025,.975),na.rm=T)
  lb<-bounds[1]
  ub<-bounds[2]
  
  output<-c(varboot,lb,ub)
  if(clust==T){
    names(output)<-c("VarBootCl","LBBootCL","UBBootCL")
  }
  
  if(clust==F){
    names(output)<-c("VarBootInd","LBBootInd","UBBootInd")
  }
  
  return(output)}

#need to add truth fix mc error variance 

#failure check function
nmetrics<-42 #need to change for different scenarios
check_fail <- function(func, ...) {
  result <- tryCatch(
    {
      #call function with args
      func(...)
    },
    error = function(e) {
      #return NA in case of an error
      return(rep(NA,nmetrics))
    }
  )
  return(result)
}


# wrapper_aspire_funv3<-function(dfsim,xlabs,m,csl,csu,sdvalout,sdvalsel,sdvalme,diffx,diffx_mod,wdoc,
#                                sdvaldout,sdvaldsel,sdvaldme){
# wrapper_aspire_fun_v3<-function(df,xlabs,clabs,diffx,cor,ind,boot){
#   
#   feature_fun_plus(df=df)
#   points<-estimator_aspire_fun_v3(df=df,xlabs=xlabs,clabs=clabs,diffx=diffx,varprint=T,cor=cor,ind=ind)
#   
#   truth<-points[1]
#   point_estimator<-points[2]
#   
#   if(is.na(point_estimator) || abs(point_estimator)>1){
#     point_estimator<-NA #dividing by 0 by chance 
#   }
#   
#   varciasym<-points[-c(1:2)]
#   
#   #if(cor==T){
#   #  varciasymcor<-varciasym[4:6]
#   #  varciasym<-varciasym[1:3]
#   #}
#   
#   if(boot==T){
#     varbootcl<-boot_aspire_v3(iters=1000,xlabs=xlabs,clust=T,diffx=diffx,dfog=df)
#     varbootind<-boot_aspire_v3(iters=1000,xlabs=xlabs,clust=F,diffx=diffx,dfog=df)}
#   
#   if(boot==F){
#     varbootcl<-c("VarBootCl"=NA,"LBBootCL"=NA,"UBBootCL"=NA) 
#     varbootind<-c("VarBootInd"=NA,"LBBootInd"=NA,"UBBootInd"=NA) 
#   }
#   
#   bias<-point_estimator-truth
#   names(bias)<-"EmpBias"
#   
#   var_express<-sum(c(cor,ind))
#   #covnames<-c("CovCL","CovCLCor","CovInd")
#   covnames<-c("CovCL","CovCLCorDf","CovCLCorT","CovInd")
#   if(var_express==1 & cor==T){
#     covnames<-covnames[-4]
#   }
#   if(var_express==1 & ind==T){
#     covnames<-covnames[-c(2:3)]
#   }
#   if(var_express==0){
#     covnames<-covnames[1]
#   }
#   var_ind<-seq(1,1+3*(var_express+sum(cor)),by=3)
#   varestasym<-varciasym[var_ind]
#   covasym<-ifelse(varciasym[var_ind+1]<truth & varciasym[var_ind+2]>truth,1,0)
#   combined_var<-as.vector(rbind(varestasym,covasym))
#   names(combined_var)<-c(rbind(names(varestasym),covnames))
#   
#   #cl boot
#   varestclboot<-varbootcl[1]
#   covbootcl<-ifelse(varbootcl[2]<truth & varbootcl[3]>truth,1,0)
#   names(covbootcl)<-c("CovBootCl")
#   
#   #ind boot 
#   varestbootind<-varbootind[1]
#   covbootind<-ifelse(varbootind[2]<truth & varbootind[3]>truth,1,0)
#   names(covbootind)<-c("CovBootInd")
#   
#   comp_vals<-feature_fun_plus(df=df)
#   
#   results<-c(comp_vals,truth,point_estimator,bias,combined_var,varestclboot,covbootcl,varestbootind,covbootind)
#   
#   #asym cor 
#   # if(cor==T){
#   #   varestasymcor<-varciasymcor[1]
#   #   covasymcor<-ifelse(varciasymcor[2]<truth & varciasymcor[3]>truth,1,0)
#   #   results<-c(truth,estimator,bias,varestasym,covasym,varestasymcor,covasymcor,varestclboot,covbootcl,varestbootind,covbootind)
#   #   name_results<-c(name_results[1:5],"VarAsymCor","CovAsymCor",name_results[6:9])
#   # }
#   return(results)
# }

#fix this so that coverage is on average at end...see if there's a difference 
wrapper_aspire_fun_v3<-function(df,xlabs,clabs,diffx,cor,ind,boot){
  
  points<-estimator_aspire_fun_v3(df=df,xlabs=xlabs,clabs=clabs,diffx=diffx,varprint=T,cor=cor,ind=ind)
  
  truth<-points[1]
  point_estimator<-points[2]
  
  if(is.na(point_estimator) || abs(point_estimator)>1){
    point_estimator<-NA #dividing by 0 by chance 
  }
  
  varciasym<-points[-c(1:2)]
  
  #if(cor==T){
  #  varciasymcor<-varciasym[4:6]
  #  varciasym<-varciasym[1:3]
  #}
  
  if(boot==T){
    #cl boot
    varbootcl<-boot_aspire_v3(iters=1000,xlabs=xlabs,clabs=clabs,clust=T,diffx=diffx,dfog=df)
    #ind boot 
    varbootind<-boot_aspire_v3(iters=1000,xlabs=xlabs,clabs=clabs,clust=F,diffx=diffx,dfog=df)}
  
  if(boot==F){
    varbootcl<-c("VarBootCl"=NA,"LBBootCL"=NA,"UBBootCL"=NA) 
    varbootind<-c("VarBootInd"=NA,"LBBootInd"=NA,"UBBootInd"=NA) 
  }
  
  comp_vals<-feature_fun_plus(df=df)
  
  results<-c(comp_vals,truth,point_estimator,varciasym,varbootcl,varbootind)
  
  return(results)
}

#length results
wrapper_aspire_funerr_v3<-function(sim,m,csl,csu,sdvalout,sdvalsel,sdvalme,xlabs,clabs,diffx,wdoc,
                                   sdvaldout,sdvaldsel,sdvaldme,dlb,dub,cor,ind,boot){
  dfsim<-data_aspire_fun_dfme(m=m,csl=csl,csu=csu,sdvalout=sdvalout,sdvalsel=sdvalsel,sdvalme=sdvalme,diffx=diffx,
                                 wdoc=wdoc,sdvaldout=sdvaldout,sdvaldsel=sdvaldsel,sdvaldme=sdvaldme,dlb=dlb,dub=dub)
  
  #print(wrapper_aspire_fun_v3(df=dfsim,xlabs=xlabs,clabs,diffx=F,cor=cor,ind=ind,boot=boot))
  
  l1<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=T,cor=cor,ind=ind,boot=boot)
  l2<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=F,cor=cor,ind=ind,boot=boot)
  return(c(l1,l2))
  
  # if(diffx==F){
  #   l1<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=T,cor=cor,ind=ind,boot=boot)
  #   l2<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=F,cor=cor,ind=ind,boot=boot)
  #   return(c(l1,l2))
  #   # l1<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=F,cor=cor,ind=ind,boot=boot)
  #   # return(c(l1,l1))
  # }
  # 
  # if(diffx==T){
  #   l1<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=T,cor=cor,ind=ind,boot=boot)
  #   l2<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=F,cor=cor,ind=ind,boot=boot)
  #   return(c(l1,l2))
  # }
}

#wrapper_aspire_funerr_v3(sim=1,m=30,csl=100,csu=300,sdvalout=.2,sdvalsel=.2,sdvalme=.2,diffx=T,wdoc=F,
#                         sdvaldout=0,sdvaldsel=0,sdvaldme=0,dlb=2,dub=3,cor=T,ind=T,xlabs=c("x1","x2","x3"))

nc<-c(30,50)
#nc<-c(12,16,20)
#csl<-c(35,75,200)
#csu<-c(55,150,300)
csl<-c(100,500)
csu<-c(300,1000)
#all models

icc_vals<-c(0.01,0.05,0.1)

#wdoc F
varval_fun<-function(icc){(icc*pi^2/3)/(1-icc)}

#wdoc T
varval_fun_wdoc<-function(icc){
  wpc<-icc
  bpc<-icc/2
  varc<--bpc*(pi^2/3)/(wpc-1)
  return(varc)
}

table1scen<-data.frame(nc=rep(30,16),csl=rep(rep(c(100,500),each=2),4),csu=rep(rep(c(300,1000),each=2),4),wdoc=F,
              icc=c(rep(rep(c(.01,.1),each=4),2)),icc_me=0,diffx=c(rep(F,8),rep(T,8)),diffx_mod=rep(c(T,F),8))
              #diffx_mod=c(rep(F,8),rep(c(T,F),4)))
scenariomat<-table1scen
#c(rep(.01,4),rep(0,12))
#scenariomat<-rbind(expand.grid(nc=nc,icc=icc_vals,icc_me=c(0),diffx=F,csl=csl,wdoc=c(F,T)),
#                   expand.grid(nc=nc,icc=icc_vals,icc_me=c(0,0.01,0.05),diffx=T,csl=csl,wdoc=c(F,T)))
#scenariomat<-cbind(scenariomat[,-5],csl=scenariomat[,5],csu=ifelse(scenariomat[,5]==100,300,1000))
#remove scenariomat doc=T/icc.05
#wdocicc5id<-which(scenariomat$wdoc==T & scenariomat$icc==.05)
#scenariomat<-scenariomat[-wdocicc5id,]

sdvals_wodoc<-sqrt(varval_fun(scenariomat[,"icc"]))
sdvalsme_wodoc<-sqrt(varval_fun(scenariomat[,"icc_me"]))

sdvals_wdoc<-sqrt(varval_fun_wdoc(scenariomat[,"icc"]))
sdvalsme_wdoc<-sqrt(varval_fun_wdoc(scenariomat[,"icc_me"]))

sdvals<-ifelse(scenariomat$wdoc==T,sdvals_wdoc,sdvals_wodoc)
sdvalsme<-ifelse(scenariomat$wdoc==T,sdvalsme_wdoc,sdvalsme_wodoc)

sdvalsd<-ifelse(scenariomat$wdoc==T,sdvals_wdoc,0)
sdvalsdme<-ifelse(scenariomat$wdoc==T,sdvalsme_wdoc,0)

#scenariomat<-scenariomat[rep(1:nrow(scenariomat),each = 2),]
#scenariomat<-cbind(scenariomat,diffx_mod=rep(c(F,T),nrow(scenariomat)/2))
#rownames(scenariomat)<-1:nrow(scenariomat)

#36 is the number of metrics 
resultsmat<-matrix(rep(0,nrow(scenariomat)*36),nrow=nrow(scenariomat))

nsim<-5000

#bootvec<-rep(F,16) #for now 
bootvec<-ifelse(scenariomat$csl==100,T,F)

#for(i in 1:nrow(scenariomat)/2){
for(i in 1:(nrow(scenariomat)/2)){
  l<-parallel::mcmapply(wrapper_aspire_funerr_v3,sim=1:nsim,
                        MoreArgs=list(m=scenariomat[2*i-1,"nc"],csl=scenariomat[2*i-1,"csl"],
                                      csu=scenariomat[2*i-1,"csu"],sdvalout=sdvals[2*i-1],sdvalsel=sdvals[2*i-1],sdvalme=sdvalsme[2*i-1],diffx=scenariomat[2*i-1,"diffx"],
                                      wdoc=scenariomat[2*i-1,"wdoc"],sdvaldout=sdvalsd[2*i-1],sdvaldsel=sdvalsd[2*i-1],sdvaldme=sdvalsdme[2*i-1],
                                      dlb=2,dub=3,cor=T,ind=T,xlabs=c("x1","x2","x3"),clabs="x4",boot=bootvec[2*i-1]),
                        mc.cores=60)
  lw<-l[1:26,]
  lwo<-l[27:52,]
  
  finalw<-rowMeans(lw[-c(seq(10,25,by=3),seq(11,26,by=3)),],na.rm=T) #for now
  convfailw<-rowSums(is.na(lw[-c(seq(10,25,by=3),seq(11,26,by=3)),]),na.rm=T)
  true_atew<-finalw[7]
  biasw<-finalw[8]-true_atew
  mevarw<-var(lw[8,],na.rm=T) #for now, variance of l
  #coverage lb
  coveragew<-apply(lw[seq(10,25,by=3),]<=true_atew & lw[seq(11,26,by=3),]>=true_atew,1,mean)
  covnames<-c("CovCL","CovCLCorDf","CovCLCorT","CovInd","CovBootCL","CovBootInd")
  names(biasw)<-"EmpBias"
  names(mevarw)<-c("MCError")
  names(convfailw)<-paste("ConvError",names(finalw),sep="")
  varcovw<-c(rbind(finalw[9:14],coveragew))
  names(varcovw)<-c(rbind(names(finalw[9:14]),covnames))
  outputw<-c(finalw[1:8],biasw,mevarw,varcovw,convfailw)
  
  finalwo<-rowMeans(lwo[-c(seq(10,25,by=3),seq(11,26,by=3)),],na.rm=T) #for now
  convfailwo<-rowSums(is.na(lwo[-c(seq(10,25,by=3),seq(11,26,by=3)),]),na.rm=T)
  true_atewo<-finalwo[7]
  biaswo<-finalwo[8]-true_atewo
  mevarwo<-var(lwo[8,],na.rm=T) #for now, variance of l
  #coverage lb
  coveragewo<-apply(lwo[seq(10,25,by=3),]<=true_atewo & lwo[seq(11,26,by=3),]>=true_atewo,1,mean)
  covnames<-c("CovCL","CovCLCorDf","CovCLCorT","CovInd","CovBootCL","CovBootInd")
  names(biaswo)<-"EmpBias"
  names(mevarwo)<-c("MCError")
  names(convfailwo)<-paste("ConvError",names(finalwo),sep="")
  varcovwo<-c(rbind(finalwo[9:14],coveragewo))
  names(varcovwo)<-c(rbind(names(finalwo[9:14]),covnames))
  outputwo<-c(finalwo[1:8],biaswo,mevarwo,varcovwo,convfailwo)
  
  resultsmat[2*i-1,]<-outputw
  resultsmat[2*i,]<-outputwo
  
  colnames(resultsmat)<-names(outputw)
  
  if(i %% 1 == 0){
    filename<-paste("aspiresimsv2both",i,".csv",sep="")
    finalvals<-cbind(scenariomat,resultsmat)
    #print(finalvals)
    write.csv(finalvals,filename)
  }
}
