#' Estimation of the ATE with measurement error and selection bias in the internal validation set.
#'
#' This function provides estimates of the ATE...
#' 
#' 
#' @param df Data frame containing all data required for estimation with user-specified names below.
#' @param diffx If True, fits a measurement error model with covariates.
#' @param xlabs Vector of covariate labels for measurement error model as characters to be interacted with treatment, usually individual level. Ignored if diffx=F.
#' @param clabs Vector of covariate labels as characters for measurement error model not to be interacted with treatment, usually at cluster level. Ignored if diffx=F.
#' @param varprint If varprint=T, returns variance of estimator. 
#' #' @param crob If crob=T, prints cluster-based variance estimates. If boot=T and crob=T, this is the non-parametric cluster bootstrap. 
#' If boot=F and crob=T, this is the cluster robust sandwich variance. If crob=F, prints iid based variances. If boot=T and crob=F, this is the non-parametric bootstrap. 
#' If boot=F and crob=F, this is the sandwich variance.
#' @param corcl If corcl=T, prints t-interval version of estimate with nc-7 df. If corcl=F, prints normal interval. If crob=F, we set corcl=F since this is intended as a method
#' for clustered data. Post-hoc corrections can be made using available corrections.
#' @param boot If boot=T, returns non-parametric bootstrap variance and percentile intervals. Otherwise, provides asymptotic variance estimates. 
#' @param iters Number of bootstrap iterations.
#' 
#' @return Point estimate for SACE and credible interval
#' 
#' @export

mesb_estim<-function(df,diffx=F,xlabs=NULL,clabs=NULL,varprint=T,crob=T,corcl=T,boot=F,iters=500){

if(crob==F & corcl==T){
  corcl<-F
  warning("No correction has been applied to individual sandwich variance.")
}
  
mesb_estim_sub<-function(df,diffx,xlabs,clabs,varprint,crob,corcl){
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
  
  if(varprint==F){
    return(c("Estimator"=finalest))
  }else{
    
    pmehat<-fitted.values(mefit)
    nparams<-ncol(dfvcov)
    
    if(crob==T){
    
    meatmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
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
    }else{
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
      }}
    
    breadmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
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
    
    if(crob==T){
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat1%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)}
    else{
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat2%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)
    }
    
    if(corcl==F){
      lb<-finalest+qnorm(.025)*sqrt(varest)
      ub<-finalest+qnorm(.975)*sqrt(varest)
      output<-c(finalest,varest,lb,ub)
      if(crob==T){
      names(output)<-c("Estimator","VarAsymCL","LBCL","UBCL")}
      else{
        names(output)<-c("Estimator","VarAsymInd","LBInd","UBInd")
      }
    }
    
    if(corcl==T){
      lb1T<-finalest+qt(.025,df=m-7)*sqrt(varest)
      ub1T<-finalest+qt(.975,df=m-7)*sqrt(varest)
      output<-c(finalest,varest,lb1T,ub1T)
      names(output)<-c("Estimator","VarAsymCL","LBCLCorT","UBCLCorT")
    }

    }
    return(output)}

if(varprint==F){
  output<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=varprint,crob=crob,corcl=corcl)
  return(output)
}

if(varprint==T & boot==F){
  var_output<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=varprint,crob=crob,corcl=corcl)
  return(var_output)
}


if(varprint==T & boot==T){
#bootstrap
boot_mesb_estim_sub<-function(dfog,iters,diffx,xlabs,clabs,crob){
  estboot<-rep(0,iters)
  for(b in 1:iters){
    if(crob==T){
      idval<-length(unique(dfog$Id))
      clustidv<-sample(1:idval,idval,replace=T)
      dflist<-vector("list",length(idval))
      for(j in 1:idval){
        dflist[[j]]<-dfog[dfog$Id==clustidv[j],]
      }
      df<-do.call(rbind,dflist)
    }
    if(crob==F){
      rowind<-sample(1:nrow(dfog),nrow(dfog),replace=T)
      df<-dfog[rowind,]
    }
    
    estim_boot<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=F) 
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
  
  if(crob==T){
    names(output)<-c("VarBootCl","LBBootCL","UBBootCL")
  }
  
  if(crob==F){
    names(output)<-c("VarBootInd","LBBootInd","UBBootInd")
  }
  
  return(output)}

finalest<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=F)
varvec<-boot_mesb_estim_sub(dfog=df,iters=iters,diffx=diffx,xlabs=xlabs,clabs=clabs,crob=crob)
output<-c(finalest,varvec)
return(output)

}}