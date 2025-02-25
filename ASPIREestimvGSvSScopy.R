#set seed 
RNGkind("L'Ecuyer-CMRG")
set.seed(01051990)

data_aspire_fun_dfme_f2<-function(m,csl,csu,small_sel=T,small_valerr=T,sdvalout,sdvalsel,sdvalme,
                                  wdoc=F,sdvaldout,sdvaldsel,sdvaldme,
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
    p0<-1/(1+exp(1-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
    p1<-1/(1+exp(1-.75-.15*valx1-.2*valx2-.15*valx3+.15*valx4-bii-bidi))
    
    y0val<-rbinom(n=csi,size=1,p0)
    y1val<-rbinom(n=csi,size=1,p1)
    
    y0[clustidr]<-y0val
    y1[clustidr]<-y1val
    
    ##me error generation
    #if doesn't depend on x
    
    #if does depend on x
    
    
    if(small_valerr==T){
      #all <10%
      ps1<-1/(1+exp(2.75*(1-y1val)-3*y1val+.35*valx1+.25*valx2-.15*valx3+.1*valx4-bimei-bidmei)) #re version
      ps0<-1/(1+exp(2*(1-y0val)-2*y0val+.55*valx1+.35*valx2-.15*valx3+.1*valx4-bimei-bidmei))}else{
        
        #all >40%
        ps1<-1/(1+exp(.2*(1-y1val)-.75*y1val+.35*valx1+.25*valx2-.15*valx3-bimei-bidmei)) #re version
        ps0<-1/(1+exp(.25*(1-y0val)-.45*y0val+.5*valx1+.35*valx2-.15*valx3-bimei-bidmei))}
    
    #measurement error probabilities
    ps<-ifelse(rep(trti,csi)==1,ps1,ps0)
    ysval<-rbinom(n=csi,size=1,p=ps)
    ys[clustidr]<-ysval
    
    ##selection and outcome models
    
    if(small_sel==T){
      #20%
      pselect1<-1/(1+exp(.15+.15*y1val+.75*valx1+.75*valx2+.75*valx3-.15*valx4-bisi-bidsi)) #make it same for x
      pselect0<-1/(1+exp(.25-.15*y0val+.75*valx1+.75*valx2+.75*valx3-.15*valx4-bisi-bidsi))}else{
        
        #40%
        pselect1<-1/(1+exp(-.45+.15*y1val+.5*valx1+.5*valx2+.5*valx3-.1*valx4-bisi-bidsi)) #make it same for x
        pselect0<-1/(1+exp(-.7-.15*y0val+.5*valx1+.5*valx2+.5*valx3-.1*valx4-bisi-bidsi))}
    
    
    pselect<-ifelse(rep(trti,csi)==1,pselect1,pselect0)
    
    #selection probability
    #pselect<-1/(1+exp(-.1-.25*trti+valx1+valx2+.75*valx3-bisi-bidsi))
    
    #- for now
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

estimator_aspire_valid<-function(df,xlabs,varprint=T,cor=T){
  N<-nrow(df)
  m<-length(unique(df$Id))
  phat<-mean(df$a)
  p_sel_form<-paste("v","~",paste("a",paste(xlabs, collapse="+"),sep="+"),sep="")
  p_sel_form<-as.formula(p_sel_form)
  p_selfit<-glm(p_sel_form,family="binomial",data=df)
  fitted_v<-fitted.values(p_selfit)
  select_v_ind<-which(df$v==1)
  fitted_v_select<-fitted_v[select_v_ind]
  dfv<-df[select_v_ind,]
  u1hat<-sum((dfv$y*dfv$a)/(fitted_v_select*phat))/N
  u0hat<-sum((dfv$y*(1-dfv$a))/(fitted_v_select*(1-phat)))/N
  finalest<-u1hat-u0hat
  truth<-mean(df$y1-df$y0)
  
  if(varprint==F){
    return(c("Estimator"=finalest))
  }else{
    
    dfcov<-cbind(int=1,df[,c("a",xlabs)])
    dfvcov<-cbind(int=1,dfv[,c("a",xlabs)])
    
    nparams<-ncol(dfcov)
    meatmat<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    breadmat<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    for(i in 1:m){
      clustid<-which(df$Id==i)
      clustidv<-which(df$Id==i & df$v==1)
      dfi<-df[clustid,]
      dfclustcov<-dfcov[clustid,]
      fitted_vi<-fitted_v[clustid]
      dfvi<-dfi[which(dfi$v==1),]
      estfunvi<-t(t(dfclustcov)%*%(dfi$v-fitted_vi))
      mi3<-sum(dfi$a-phat)
      mi4<-sum(dfvi$y*dfvi$a/(fitted_v[clustidv]))-nrow(dfi)*phat*u1hat
      mi5<-sum(dfvi$y*(1-dfvi$a)/(fitted_v[clustidv]))-nrow(dfi)*(1-phat)*u0hat
      meatmat<-meatmat+c(estfunvi,mi3,mi4,mi5)%*%t(c(estfunvi,mi3,mi4,mi5))
    }
    
    #beta beta
    breadmat[1:nparams,1:nparams]<--t(dfcov*fitted_v*(1-fitted_v))%*%as.matrix(dfcov)
    
    #pi pi
    breadmat[nparams+1,nparams+1]<--N
    
    #u1 beta
    breadmat[nparams+2,1:nparams]<--apply(dfvcov*dfv$y*dfv$a*(1-fitted_v_select)/((fitted_v_select)),2,sum)
    
    #u1 pi
    breadmat[nparams+2,nparams+1]<--N*u1hat
    
    #u1 u1
    breadmat[nparams+2,nparams+2]<--N*phat
    
    #u0 beta
    breadmat[nparams+3,1:nparams]<--apply(dfvcov*dfv$y*(1-dfv$a)*(1-fitted_v_select)/((fitted_v_select)),2,sum)
    
    #u0 phat
    breadmat[nparams+3,nparams+1]<-N*u0hat
    
    #u0 u0
    breadmat[nparams+3,nparams+3]<--N*(1-phat)
    
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat)%*%meatmat%*%solve(t(breadmat))%*%c(rep(0,nparams+1),1,-1)
    
    lb<-finalest+qnorm(.025)*sqrt(varest)
    ub<-finalest+qnorm(.975)*sqrt(varest)
    
    names_og<-c("Truth","Estimator","VarAsymCL","LBCL","UBCL")
    
    if(cor==F){
      output<-c(truth,finalest,varest,lb,ub)
      names(output)<-names_og
    }
    
    if(cor==T){
      varest1<-m/(m-4)*varest
      lb1Df<-finalest+qnorm(.025)*sqrt(varest1)
      ub1Df<-finalest+qnorm(.975)*sqrt(varest1)
      
      lb1T<-finalest+qt(.025,df=m-4)*sqrt(varest)
      ub1T<-finalest+qt(.975,df=m-4)*sqrt(varest)
      
      output<-c(truth,finalest,varest,lb,ub,varest1,lb1Df,ub1Df,varest,lb1T,ub1T)
      names(output)<-c(names_og,"VarAsymCLCorDf","LBCLCorDf","UBCLCorDf",
                       "VarAsymCLCorT","LBCLCorT","UBCLCorT")
    }}
  return(output)}

nmetrics_outonly<-1
check_fail <- function(func, ...) {
  result <- tryCatch(
    {
      #call function with args
      func(...)
    },
    error = function(e) {
      #return NA in case of an error
      return(rep(NA,nmetrics_outonly))
    }
  )
  return(result)
}

wrapper_aspire_funerr_v3_f2<-function(sim,m,csl,csu,small_sel,small_valerr,
                                      sdvalout,sdvalsel,sdvalme,xlabs,clabs,wdoc,
                                      sdvaldout,sdvaldsel,sdvaldme,dlb,dub,cor,ind){
 
   dfsim<-data_aspire_fun_dfme_f2(m=m,csl=csl,csu=csu,small_sel=small_sel,small_valerr=small_valerr,
                                  sdvalout=sdvalout,sdvalsel=sdvalsel,sdvalme=sdvalme,wdoc=wdoc,
                                  sdvaldout=sdvaldout,sdvaldsel=sdvaldsel,sdvaldme=sdvaldme,dlb=2,dub=3)
   
   specs<-feature_fun_plus(dfsim)

   truth<-mean(dfsim$y1-dfsim$y0)
   l1<-check_fail(func=estimator_aspire_fun_v3,df=dfsim,xlabs=xlabs,clabs=clabs,varprint=F,diffx=T,cor=cor,ind=ind)
   l2<-check_fail(func=estimator_aspire_valid,df=dfsim,xlabs=c(xlabs,clabs),varprint=F,cor=cor)
   l3<-mean(dfsim$ys*dfsim$a)-mean(dfsim$ys*(1-dfsim$a))

   estimate_opts<-c(truth,l1,l2,l3)
   names(estimate_opts)<-c("TrueATE","Proposed_Estim","IPSW_Estim","SSOnly_Estim")
   return(c(specs,estimate_opts))
   }

dgm_f2_scens<-expand.grid(nc=30,csl=100,csu=300,icc=.01,icc_me=0,wdoc=F,
                          "small_sel"=c(T,F),"small_valerr"=c(T,F))
scenariomat<-dgm_f2_scens

#wdoc F
varval_fun<-function(icc){(icc*pi^2/3)/(1-icc)}

#wdoc T
varval_fun_wdoc<-function(icc){
  wpc<-icc
  bpc<-icc/2
  varc<--bpc*(pi^2/3)/(wpc-1)
  return(varc)
}

sdvals_wodoc<-sqrt(varval_fun(scenariomat[,"icc"]))
sdvalsme_wodoc<-sqrt(varval_fun(scenariomat[,"icc_me"]))

sdvals_wdoc<-sqrt(varval_fun_wdoc(scenariomat[,"icc"]))
sdvalsme_wdoc<-sqrt(varval_fun_wdoc(scenariomat[,"icc_me"]))

sdvals<-ifelse(scenariomat$wdoc==T,sdvals_wdoc,sdvals_wodoc)
sdvalsme<-ifelse(scenariomat$wdoc==T,sdvalsme_wdoc,sdvalsme_wodoc)

sdvalsd<-ifelse(scenariomat$wdoc==T,sdvals_wdoc,0)
sdvalsdme<-ifelse(scenariomat$wdoc==T,sdvalsme_wdoc,0)

nsim<-5000
res_list <- vector("list", length = 4)
names(res_list)<-c("Small_Sel,Small_ValErr","Large_Sel,Small_ValErr",
                   "Small_Sel,Large_ValErr", "Large_Sel,Large_ValErr")

for(i in 1:nrow(scenariomat)){
  
  l<-parallel::mcmapply(wrapper_aspire_funerr_v3_f2,sim=1:nsim,
                   MoreArgs=list(m=scenariomat[i,"nc"],csl=scenariomat[i,"csl"],csu=scenariomat[i,"csu"],
                                 small_sel=scenariomat[i,"small_sel"],small_valerr=scenariomat[i,"small_valerr"],
                                 sdvalout=sdvals[i],sdvalsel=sdvals[i],sdvalme=sdvalsme[i],
                                 wdoc=scenariomat[i,"wdoc"],sdvaldout=sdvalsd[i],sdvaldsel=sdvalsd[i],sdvaldme=sdvalsdme[i],
                                 dlb=2,dub=3,cor=T,ind=T,xlabs=c("x1","x2","x3"),clabs="x4"),mc.cores=60)
  res_list[[i]]<-l
  
  if(i %% 1 == 0){
    filename<-paste("aspiref2",i,".rds",sep="")
    #finalvals<-cbind(scenariomat,resultsmat)
    #print(res_list)
    saveRDS(res_list,filename)
  }
}
