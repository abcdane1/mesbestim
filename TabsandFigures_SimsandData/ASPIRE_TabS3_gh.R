#set seed 
RNGkind("L'Ecuyer-CMRG")
set.seed(01051990)

data_aspire_fun_dfme<-function(m,csl,csu,sdvalout,sdvalsel,sdvalme,diffx=T,wdoc=F,
                               sdvaldout,sdvaldsel,sdvaldme,
                               dlb=2,dub=3){
  
  #cluster size 
  cs<-sample(csl:csu,m,replace=T)
  
  #treatment assignment
  trt<-rbinom(m,1,p=.5)
  
  #individual covariates
  x1<-rep(0,sum(cs))
  x2<-rep(0,sum(cs))
  x3<-rep(0,sum(cs))
  x4<-rep(0,sum(cs))
  
  #gold-standard obs outcome
  y<-rep(0,sum(cs))
  
  #silver-standard obs outcome
  ys<-rep(0,sum(cs))
  
  #selection variable
  v<-rep(0,sum(cs))
  
  #potential outcomes (for TRUE ATE)
  y0<-rep(0,sum(cs))
  y1<-rep(0,sum(cs))
  
  #clinician within cluster/site
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
  
  
  #induce within cluster correlation via random intercept (like in glmm)
  #if 0, return 0 vector 
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
  
  #clinician sds are 0 if no correlation
  if(wdoc==F){
    sdvaldout<-0
    sdvaldsel<-0
    sdvaldme<-0}
  
  #clinician level random effects 
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
    #clinician random effects
    bidi<-bid[doc[clustidr]]
    bidsi<-bids[doc[clustidr]]
    bidmei<-bidme[doc[clustidr]]
    
    #fill individual covariates
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
    p0<-1/(1+exp(1-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
    p1<-1/(1+exp(1-.75-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
    
    y0val<-rbinom(n=csi,size=1,p0)
    y1val<-rbinom(n=csi,size=1,p1)
    
    y0[clustidr]<-y0val
    y1[clustidr]<-y1val
    
    ##misclassification generation
    #if doesn't depend on x
    if(diffx==F){
      ps1<-1/(1+exp(1-2.5*y1val-bimei-bidmei))
      ps0<-1/(1+exp(1.25-1.5*y0val-bimei-bidmei))
    }
    
    #if does depend on x
    if(diffx==T){
      ps1<-1/(1+exp(.75-2.5*y1val+.25*valx1+.15*valx2+.25*valx3-.1*valx4-bimei-bidmei)) 
      ps0<-1/(1+exp(1.25-1.5*y0val-.25*valx1+.25*valx2+.15*valx3-.1*valx4-bimei-bidmei))
    }
    
    #measurement error probabilities
    ps<-ifelse(rep(trti,csi)==1,ps1,ps0)
    ysval<-rbinom(n=csi,size=1,p=ps)
    ys[clustidr]<-ysval
    
    ##selection models and observed data 
    
    #selection probability
    pselect1<-1/(1+exp(.5-.15*y1val+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi))
    pselect0<-1/(1+exp(.25+.15*y0val+.5*valx1+.5*valx2-.25*valx3+.25*valx4-bisi-bidsi))
    
    #observed selection probability under treatment 
    pselect<-ifelse(rep(trti,csi)==1,pselect1,pselect0)
    
    #potentially observed y values by treatment 
    yval<-ifelse(rep(trti,csi)==1,y1val,y0val)
    
    #selection indicators under treatment
    selectind<-rbinom(n=csi,size=1,pselect)
    v[clustidr]<-selectind
    
    #observed gold-standard outcomes under treatment 
    y[clustidr]<-ifelse(selectind==1,yval,NA)
    
    c<-c+csi}
  df<-data.frame(y=y,ys=ys,Id=rep(1:m,cs),Doc=doc,a=rep(trt,cs),x1=x1,x2=x2,x3=x3,x4=x4,v=v,y0=y0,y1=y1)
  return(df)}


##empirical selection and mis-measurement proportions
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

##SSW estimator
estimator_aspire_fun_v3<-function(df,xlabs,clabs,varprint=T,diffx=T,cor=T,ind=T){
  N<-nrow(df)
  #validation subset
  dfv<-df[which(df$v==1),]
  
  #regressor matrices for classification probabilities with set a and y
  new111<-data.frame(df[,c(xlabs,clabs)],a=1,y=1)
  new101<-data.frame(df[,c(xlabs,clabs)],a=1,y=0)
  new110<-data.frame(df[,c(xlabs,clabs)],a=0,y=1)
  new100<-data.frame(df[,c(xlabs,clabs)],a=0,y=0)
  
  #if classification model depends on covariates
  if(diffx==T){
    #classification model adjusts for outcome, treatment, xlabs, clabs, treatment*xlabs, and treatment*y
    formdiffx<-paste("ys","~",paste("a","y",paste(c(xlabs,clabs),
                                                  collapse="+"),paste(paste("a",c(xlabs,"y"),sep="*"),collapse="+"),sep="+"),sep="")
    formdiffx<-as.formula(formdiffx)
    mefit<-glm(formdiffx,family="binomial",data=dfv)
    
    #selection design matrix
    dfvcov<-cbind(int=1,dfv[,c("a","y",xlabs,clabs)],dfv$a*dfv[,c(xlabs,"y")])
    
    #overall design matrix with set a and y
    df111cov<-cbind(int=1,new111[,c("a","y",xlabs,clabs)],new111[,c(xlabs,"y")])
    df101cov<-cbind(int=1,new101[,c("a","y",xlabs,clabs)],new101[,c(xlabs,"y")])
    df110cov<-cbind(int=1,new110[,c("a","y",xlabs,clabs)],0*new110[,c(xlabs,"y")])
    df100cov<-cbind(int=1,new100[,c("a","y",xlabs,clabs)],0*new100[,c(xlabs,"y")])
  }
  
  #if classification model does not depend on covariates
  if(diffx==F){
    #saturated classification model adjusted for outcome, treatment, outcome*treatment
    formsamex<-paste("ys","~",paste("a","y",paste(paste("a","y",sep="*"),collapse="+"),sep="+"),sep="")
    formsamex<-as.formula(formsamex)
    mefit<-glm(formsamex,family="binomial",data=dfv)
    
    #selection design matrix
    dfvcov<-cbind(int=1,dfv[,c("a","y")],dfv$a*dfv[,"y"])
    
    #overall design matrix with set a and y
    df111cov<-cbind(int=1,new111[,c("a","y")],new111[,"y"])
    df101cov<-cbind(int=1,new101[,c("a","y")],new101[,"y"])
    df110cov<-cbind(int=1,new110[,c("a","y")],0*new110[,"y"])
    df100cov<-cbind(int=1,new100[,c("a","y")],0*new100[,"y"])
  }
  
  #predicted classification probabilities 
  p111hat<-predict(mefit,newdata=new111,type="response")
  p101hat<-predict(mefit,newdata=new101,type="response")
  p110hat<-predict(mefit,newdata=new110,type="response")
  p100hat<-predict(mefit,newdata=new100,type="response")
  
  #estimated proportion of individuals under treatment 1
  phat<-mean(df$a)
  
  #estimated values of mu1 and mu0
  u1hat<-1/N*sum((df$ys*df$a-phat*p101hat)/(phat*(p111hat-p101hat)))
  u0hat<-1/N*sum((df$ys*(1-df$a)-(1-phat)*p100hat)/((1-phat)*(p110hat-p100hat)))
  
  #SSW estimator 
  finalest<-u1hat-u0hat
  
  #difference in empirical means of potential outcomes
  truth<-mean(df$y1-df$y0)
  
  #just print SSW estimator
  if(varprint==F){
    return(c("Estimator"=finalest))
  }else{
    
    #for asymptotic variance estimation
    
    #fitted values 
    pmehat<-fitted.values(mefit)
    
    #number of parameters 
    nparams<-ncol(dfvcov)
    
    #meat matrix
    meatmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    #bread matrix 
    breadmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    m<-length(unique(df$Id))
    
    #unbiased estimating equations 
    for(i in 1:m){
      #gee on selection subset 
      clustid<-which(dfv$Id==i)
      dfviclustcov<-dfvcov[clustid,]
      pmehati<-pmehat[clustid]
      ysi<-dfv[clustid,]$ys
      estfuni<-t(t(dfviclustcov)%*%(ysi-pmehati)) #this is the estimating function by group, the phi_i, summing over j within the group
      
      #remaining components
      clustidfull<-which(df$Id==i)
      dfi<-df[clustidfull,]
      
      #mic for phat 
      mi5x<-sum(dfi$a-phat)
      #mic u1hat 
      mi6x<-sum((dfi$ys*dfi$a-phat*p101hat[clustidfull])/(phat*(p111hat[clustidfull]-p101hat[clustidfull]))-u1hat)
      #mic u0hat 
      mi7x<-sum((dfi$ys*(1-dfi$a)-(1-phat)*p100hat[clustidfull])/((1-phat)*(p110hat[clustidfull]-p100hat[clustidfull]))-u0hat)
      meatmat1<-meatmat1+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
    }
    
    ##filling components of bread matrix, sum of jacobian matrices 
    
    #theta theta
    breadmat1[1:nparams,1:nparams]<--t(dfvcov*pmehat*(1-pmehat))%*%as.matrix(dfvcov)
    
    #pi pi
    breadmat1[nparams+1,nparams+1]<--N
    
    #mu1 theta
    breadmat1[nparams+2,1:nparams]<-apply((df101cov*phat^2*(-p101hat)*(1-p101hat)*(p111hat-p101hat)-(df$a*df$ys-phat*p101hat)*phat*(df111cov*p111hat*(1-p111hat)-df101cov*p101hat*(1-p101hat)))/
                                            (phat^2*(p111hat-p101hat)^2),2,sum)
    #mu1 pi
    breadmat1[nparams+2,nparams+1]<-sum((-p101hat*(phat)-(df$a*df$ys-phat*p101hat))/(phat^2*(p111hat-p101hat)))
    
    #mu1 mu1
    breadmat1[nparams+2,nparams+2]<--N
    
    #mu0 theta
    breadmat1[nparams+3,1:nparams]<-apply((df100cov*(1-phat)^2*(-p100hat)*(1-p100hat)*(p110hat-p100hat)-((1-df$a)*df$ys-(1-phat)*p100hat)*(1-phat)*(df110cov*p110hat*(1-p110hat)-df100cov*p100hat*(1-p100hat)))/
                                            ((1-phat)^2*(p110hat-p100hat)^2),2,sum)
    #mu0 pi
    breadmat1[nparams+3,nparams+1]<-sum((p100hat*((1-phat))+((1-df$a)*df$ys-(1-phat)*p100hat))/((1-phat)^2*(p110hat-p100hat)))
    
    #mu0 mu0
    breadmat1[nparams+3,nparams+3]<--N
    
    #variance est
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat1%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)
    
    #wald interval
    lb<-finalest+qnorm(.025)*sqrt(varest)
    ub<-finalest+qnorm(.975)*sqrt(varest)
    
    names_og<-c("Truth","Estimator","VarAsymCL","LBCL","UBCL")
    
    #if no correction leave as is
    if(cor==F){
      output<-c(truth,finalest,varest,lb,ub)
      names(output)<-names_og
    }
    
    #if correction include two df corrections 
    if(cor==T){
      #df correction wald
      varest1<-m/(m-7)*varest
      lb1Df<-finalest+qnorm(.025)*sqrt(varest1)
      ub1Df<-finalest+qnorm(.975)*sqrt(varest1)
      
      #df correction t-interval (only one referenced)
      lb1T<-finalest+qt(.025,df=m-7)*sqrt(varest)
      ub1T<-finalest+qt(.975,df=m-7)*sqrt(varest)
      
      output<-c(truth,finalest,varest,lb,ub,varest1,lb1Df,ub1Df,varest,lb1T,ub1T)
      names(output)<-c(names_og,"VarAsymCLCorDf","LBCLCorDf","UBCLCorDf",
                       "VarAsymCLCorT","LBCLCorT","UBCLCorT")
    }
    
    #just individual version of above, indexing by ij for meat matrix
    if(ind==T){
      meatmat2<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
      for(i in 1:N){
        dfi<-df[i,]
        if(dfi$v==1){
          idv<-which(row.names(dfv)==i)
          pmehati<-pmehat[idv]
          dfvicov<-dfvcov[idv,]
          ysi<-dfv[idv,]$ys
          estfuni<-t(t(dfvicov)%*%(ysi-pmehati))}else{
            estfuni<-rep(0,nparams)} #this is the estimating function by person, the phi_ij, summing over ij
        
        mi5x<-sum(dfi$a-phat)
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

#bootstrap (for clusters 1-m)
boot_aspire_v3<-function(iters=1000,xlabs,clabs,clust=T,diffx=T,dfog){
  estboot<-rep(0,iters)
  #failure counter 
  failboot<-rep(0,iters)
  for(b in 1:iters){
    #if cluster boot, resample with replacement by cluster index
    if(clust==T){
      idval<-length(unique(dfog$Id))
      clustidv<-sample(1:idval,idval,replace=T)
      dflist<-vector("list",idval)
      for(j in 1:idval){
        dflist[[j]]<-dfog[dfog$Id==clustidv[j],]
      }
      df<-do.call(rbind,dflist)
    }
    #if not cluster boot, resample with replacement by individual index 
    if(clust==F){
      rowind<-sample(1:nrow(dfog),nrow(dfog),replace=T)
      df<-dfog[rowind,]
    }
    
    estim_boot<-estimator_aspire_fun_v3(df,xlabs=xlabs,clabs=clabs,varprint=F,diffx=diffx) 
    
    #avoid degenerate fits
    if(!is.na(estim_boot) && abs(estim_boot)<=1){
      estboot[b]<-estim_boot[1]
      failboot[b]<-0
    }else{
      estboot[b]<-NA
      failboot[b]<-1
    }}
  
  #estimator for variance 
  varboot<-var(estboot,na.rm=T) 
  bounds<-quantile(estboot,c(.025,.975),na.rm=T)
  lb<-bounds[1]
  ub<-bounds[2]
  
  #proportion of failures in bootstrap round
  failboot_prop<-mean(failboot)
  
  #outputs
  output<-c(varboot,lb,ub,failboot_prop)
  if(clust==T){
    names(output)<-c("VarBootCl","LBBootCL","UBBootCL","FailBootCL")
  }
  
  if(clust==F){
    names(output)<-c("VarBootInd","LBBootInd","UBBootInd","FailBootInd")
  }
  
  return(output)}

#failure check function
nmetrics<-28 #need to change for different scenarios
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


#wrapper function to allow for entries
wrapper_aspire_fun_v3<-function(df,xlabs,clabs,diffx,cor,ind,boot){
  
  #estimation vector
  points<-estimator_aspire_fun_v3(df=df,xlabs=xlabs,clabs=clabs,diffx=diffx,varprint=T,cor=cor,ind=ind)
  
  #difference in mean potential outcomes 
  truth<-points[1]
  
  #SSW estimator
  point_estimator<-points[2]
  
  #if either the point estimator is bigger than 1 or NA
  if(is.na(point_estimator) || abs(point_estimator)>1){
    point_estimator<-NA
  }
  
  #asymptotic variance vectors
  varciasym<-points[-c(1:2)]
  
  if(boot==T){
    #cl boot
    boot_cl<-boot_aspire_v3(iters=1000,xlabs=xlabs,clabs=clabs,clust=T,diffx=diffx,dfog=df)
    #ind boot 
    boot_ind<-boot_aspire_v3(iters=1000,xlabs=xlabs,clabs=clabs,clust=F,diffx=diffx,dfog=df)}
  
  if(boot==F){
    boot_cl<-c("VarBootCl"=NA,"LBBootCL"=NA,"UBBootCL"=NA,"FailBootCL"=NA) 
    boot_ind<-c("VarBootInd"=NA,"LBBootInd"=NA,"UBBootInd"=NA,"FailBootInd"=NA) 
  }
  
  comp_vals<-feature_fun_plus(df=df)
  
  results<-c(comp_vals,truth,point_estimator,varciasym,boot_cl,boot_ind)
  
  return(results)
}

#wrapper function to allow for parallelization across simulated data sets and failed outputs
wrapper_aspire_funerr_v3<-function(sim,m,csl,csu,sdvalout,sdvalsel,sdvalme,xlabs,clabs,diffx,wdoc,
                                   sdvaldout,sdvaldsel,sdvaldme,dlb,dub,cor,ind,boot){
  #data generation
  dfsim<-data_aspire_fun_dfme(m=m,csl=csl,csu=csu,sdvalout=sdvalout,sdvalsel=sdvalsel,sdvalme=sdvalme,diffx=diffx,
                              wdoc=wdoc,sdvaldout=sdvaldout,sdvaldsel=sdvaldsel,sdvaldme=sdvaldme,dlb=dlb,dub=dub)
  
  #classification model includes covs
  l1<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=T,cor=cor,ind=ind,boot=boot)
  
  #classification model does not include covs
  l2<-check_fail(func=wrapper_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,diffx=F,cor=cor,ind=ind,boot=boot)
  
  return(c(l1,l2))}

##icc to variance conversion

#if wdoc F, no clinician random intercept
varval_fun<-function(icc){(icc*pi^2/3)/(1-icc)}

#if wdoc T, clinician random intercept 
varval_fun_wdoc<-function(icc){
  wpc<-icc
  bpc<-icc/2
  varc<--bpc*(pi^2/3)/(wpc-1)
  return(varc)
}

##table 1 scenarios, error function added
table1scen<-data.frame(nc=rep(30,8),csl=rep(rep(c(100,500),each=2),2),csu=rep(rep(c(300,1000),each=2),2),wdoc=F,
                       icc=rep(.2,8),icc_me=0,diffx=c(rep(F,4),rep(T,4)),diffx_mod=rep(c(T,F),4))
scenariomat<-table1scen

#table 3 scenarios
# table1scen<-data.frame(nc=rep(30,8),csl=rep(rep(c(100,500),each=2),2),csu=rep(rep(c(300,1000),each=2),2),wdoc=F,
#                        icc=rep(.2,8),icc_me=rep(.2,8),diffx=c(rep(F,4),rep(T,4)),diffx_mod=rep(c(T,F),4))
# scenariomat<-table1scen

##computing standard deviations for data generation function
sdvals_wodoc<-sqrt(varval_fun(scenariomat[,"icc"]))
sdvalsme_wodoc<-sqrt(varval_fun(scenariomat[,"icc_me"]))

sdvals_wdoc<-sqrt(varval_fun_wdoc(scenariomat[,"icc"]))
sdvalsme_wdoc<-sqrt(varval_fun_wdoc(scenariomat[,"icc_me"]))

sdvals<-ifelse(scenariomat$wdoc==T,sdvals_wdoc,sdvals_wodoc)
sdvalsme<-ifelse(scenariomat$wdoc==T,sdvalsme_wdoc,sdvalsme_wodoc)

#at the clinician level, will be 0 if wdoc=F
sdvalsd<-ifelse(scenariomat$wdoc==T,sdvals_wdoc,0)
sdvalsdme<-ifelse(scenariomat$wdoc==T,sdvalsme_wdoc,0)

#38 is the number of metrics 
resultsmat<-matrix(rep(0,nrow(scenariomat)*38),nrow=nrow(scenariomat))

nsim<-5000

#bootvec<-rep(F,16)
bootvec<-ifelse(scenariomat$csl==100,T,F)

#obtaining output of simulations
for(i in 1:(nrow(scenariomat)/2)){
  l<-parallel::mcmapply(wrapper_aspire_funerr_v3,sim=1:nsim,
                        MoreArgs=list(m=scenariomat[2*i-1,"nc"],csl=scenariomat[2*i-1,"csl"],
                                      csu=scenariomat[2*i-1,"csu"],sdvalout=sdvals[2*i-1],sdvalsel=sdvals[2*i-1],sdvalme=sdvalsme[2*i-1],diffx=scenariomat[2*i-1,"diffx"],
                                      wdoc=scenariomat[2*i-1,"wdoc"],sdvaldout=sdvalsd[2*i-1],sdvaldsel=sdvalsd[2*i-1],sdvaldme=sdvalsdme[2*i-1],
                                      dlb=2,dub=3,cor=T,ind=T,xlabs=c("x1","x2","x3"),clabs="x4",boot=bootvec[2*i-1]),
                        mc.cores=60)
  #model 1 (classification with covariates) results
  lw<-l[1:28,]
  lw<-rbind(lw[-c(24,28),],lw[c(24,28),])
  
  #model 2 (classification without covariates) results
  lwo<-l[29:56,]
  lwo<-rbind(lwo[-c(24,28),],lwo[c(24,28),])
  
  #mean of estimates 
  finalw<-rowMeans(lw[-c(seq(10,25,by=3),seq(11,26,by=3),27,28),],na.rm=T) #for now
  
  #sum of convergence fails
  convfailw<-rowSums(is.na(lw[-c(seq(10,25,by=3),seq(11,26,by=3),27,28),]),na.rm=T)
  
  #true ate
  true_atew<-finalw[7]
  
  #bias
  biasw<-finalw[8]-true_atew
  
  #empirical variance 
  mevarw<-var(lw[8,],na.rm=T)
  
  #coverage 
  coveragew<-apply(lw[seq(10,25,by=3),]<=true_atew & lw[seq(11,26,by=3),]>=true_atew,1,mean)
  covnames<-c("CovCL","CovCLCorDf","CovCLCorT","CovInd","CovBootCL","CovBootInd")
  names(biasw)<-"EmpBias"
  names(mevarw)<-c("MCError")
  names(convfailw)<-paste("ConvError",names(finalw),sep="")
  varcovw<-c(rbind(finalw[9:14],coveragew))
  names(varcovw)<-c(rbind(names(finalw[9:14]),covnames))
  
  #failure rate
  bootfailclw<-mean(lw[27,])
  names(bootfailclw)<-"AvgBootCLFail"
  bootfailindw<-mean(lw[28,])
  names(bootfailindw)<-"AvgBootIndFail"
  
  #output
  outputw<-c(finalw[1:8],biasw,mevarw,varcovw,convfailw,bootfailclw,bootfailindw)
  
  #repeat of above without covariates (model 2)
  finalwo<-rowMeans(lwo[-c(seq(10,25,by=3),seq(11,26,by=3),27,28),],na.rm=T)
  convfailwo<-rowSums(is.na(lwo[-c(seq(10,25,by=3),seq(11,26,by=3),27,28),]),na.rm=T)
  true_atewo<-finalwo[7]
  biaswo<-finalwo[8]-true_atewo
  mevarwo<-var(lwo[8,],na.rm=T)
  coveragewo<-apply(lwo[seq(10,25,by=3),]<=true_atewo & lwo[seq(11,26,by=3),]>=true_atewo,1,mean)
  covnames<-c("CovCL","CovCLCorDf","CovCLCorT","CovInd","CovBootCL","CovBootInd")
  names(biaswo)<-"EmpBias"
  names(mevarwo)<-c("MCError")
  names(convfailwo)<-paste("ConvError",names(finalwo),sep="")
  varcovwo<-c(rbind(finalwo[9:14],coveragewo))
  names(varcovwo)<-c(rbind(names(finalwo[9:14]),covnames))
  bootfailclwo<-mean(lwo[27,])
  names(bootfailclwo)<-"AvgBootCLFail"
  bootfailindwo<-mean(lwo[28,])
  names(bootfailindwo)<-"AvgBootIndFail"
  outputwo<-c(finalwo[1:8],biaswo,mevarwo,varcovwo,convfailwo,bootfailclwo,bootfailindwo)
  
  resultsmat[2*i-1,]<-outputw
  resultsmat[2*i,]<-outputwo
  
  colnames(resultsmat)<-names(outputw)
  
  #save file for every ith iteration
  if(i %% 1 == 0){
    filename<-paste("aspiresimsv2_newicc",i,".csv",sep="")
    finalvals<-cbind(scenariomat,resultsmat)
    #print(finalvals)
    write.csv(finalvals,filename)
  }
}
