#' Estimation of the average treatment effect with misclassified outcomes and non-random validation subsets using cluster-randomized data
#'
#' This function estimates the average treatment effect (ATE) for a binary outcome using silver-standard outcome measures, which are available for all individuals but subject to misclassification, 
#' together with a non-randomly selected validation subset of gold-standard outcome measures. By default, the function implements inference for cluster-randomized data (as in Isenberg et al. 2025), 
#' but it can also be specified for the i.i.d. setting. The classification model, fit on the internal validation subset, 
#' regresses gold-standard outcomes on silver-standard outcomes, treatment, and if user-specified, covariates. If covariates are specified, the classification model is assumed to be a logistic GEE with an independence working correlation. 
#' If no covariates are specified, the ATE is completely non-parametrically estimated (saturated GEE). Variance options for cluster-randomized data allow for an asymptotic cluster-robust sandwich variance estimate with or without small-sample correction (t-interval, with #clusters-7 df)
#' or the non-parametric cluster bootstrap. Corresponding options for individually-randomized data are provided as well, but no small sample corrections are included for asymptotic variance estimates.
#' 
#' 
#' @param df Data frame containing all columns required for estimation.
#' @param trt Name of treatment as character (required). 
#' @param ssvar Name of silver-standard outcome as character (required). 
#' @param gsvar Name of gold-standard outcome as character (required).
#' @param vald Name of validation variable as character (required).
#' @param cl Name of cluster id as character (required for clustered data). 
#' @param xlabs Vector of covariate labels as characters for the classification model that are to be interacted with treatment. Categorical variables must be dummy coded. Interactions that are also interacted with treatment must appear as additional columns in your data frame whose names are included in `xlabs`.
#' @param clabs Vector of covariate labels as characters for classification model that are not to be interacted with treatment. While these
#' can be any covariates, it is recommended in very small samples that they should include observed cluster variables; hence, the naming convention. Categorical variables must be dummy coded.
#' Additional interactions not interacted with treatment must appear as additional columns in your data frame whose names are included in `clabs`.
#' @param varprint If `varprint=T`, returns variance of estimator. 
#' @param crob If ``crob=T``, returns cluster-based variance estimates. If `boot=T` and `crob=T`, this is the non-parametric cluster bootstrap. 
#' If `boot=F` and `crob=T`, this is the cluster-robust sandwich variance. If `crob=F`, returns iid-based variances. If `boot=T` and `crob=F`, this is the non-parametric bootstrap. 
#' If `boot=F` and `crob=F`, this is the sandwich variance.
#' @param corcl If `corcl=T`, returns t-interval version of estimate with #clusters-7 degrees of freedom. If `corcl=F`, prints normal interval. If `crob=F`, we set `corcl=F` since this is intended as a method
#' for clustered data, and post-hoc corrections can be made by the user.
#' @param boot If `boot=T`, returns non-parametric bootstrap variance and percentile intervals. Otherwise, provides asymptotic variance estimates. 
#' @param iters Number of bootstrap iterations. Ignored if `boot=F`.
#' @param alpha Value strictly between 0 and 1 for significance level (i.e, 1-`alpha/2` confidence interval). Critical values and percentile intervals are based on \{`alpha`/2, 1-`alpha`/2\}.
#' 
#' @return Point and interval estimates for the ATE.  
#' 
#' @references
#' \enumerate{
#'    \item{Isenberg, D., Mitra, N., Marcus, S.C., Beidas, R.S., and Linn, K.A, 2025. Estimating the average treatment effect in cluster-randomized trials with misclassified outcomes and non-random validation subsets. \emph{In Progress}}
#'    \item{Shu, D., and Yi, G.Y., 2019. Causal inference with measurement error in outcomes: Bias analysis and estimation methods. \emph{Statistical Methods in Medical Research, 28(7)}, pp.2049-2068.}
#'    \item{Shen, J., Isenberg, D., Linn, K.A., and Hubbard, R.A., 2025. Integrating Misclassified EHR Outcomes With Validated Outcomes From a Non‐Probability Sample. \emph{Statistics in Medicine, 44(15-17)}, p.e70127.}
#'    \item{Stefanski L.A. and Boos D.D., 2002. The calculus of M-estimation. \emph{The American Statistician, 56(1)}, pp.29–38.}
#'    \item{Field C.A. and Welsh A.H., 2007. Bootstrapping clustered data. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology, 69(3)}, pp.369–390.}
#' }
#' 
#' @export

#ATE estimator function 
mesb_estim<-function(df,trt,ssvar,gsvar,vald,cl,
                     xlabs=NULL,clabs=NULL,varprint=T,crob=T,corcl=T,boot=F,iters=500,alpha=.05){

#rename variables for simplicity 
colnames(df)[colnames(df) == trt] <- "a"
colnames(df)[colnames(df) == ssvar] <- "ys"
colnames(df)[colnames(df) == gsvar] <- "y"
colnames(df)[colnames(df) == vald] <- "v"
colnames(df)[colnames(df) == cl] <- "Id"

#if not cluster robust, no correction
if(crob==F & corcl==T & boot==F){
  corcl<-F
  warning("No correction has been applied to individual sandwich variance.")
}


if (is.null(xlabs) && is.null(clabs)) {
  diffx<-F
}else{
  diffx<-T
}
  
#SSW estimator function
mesb_estim_sub<-function(df,diffx=diffx,xlabs,clabs,varprint,crob,corcl){
  N<-nrow(df)
  
  #validation subset
  dfv<-df[which(df$v==1),]
  
  #regressor matrices for classification probabilities with set a and y
  new111<-data.frame(df[,c(xlabs,clabs),drop=F],a=1,y=1)
  new101<-data.frame(df[,c(xlabs,clabs),drop=F],a=1,y=0)
  new110<-data.frame(df[,c(xlabs,clabs),drop=F],a=0,y=1)
  new100<-data.frame(df[,c(xlabs,clabs),drop=F],a=0,y=0)
  
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
    df100cov<-cbind(int=1,new100[,c("a","y",xlabs,clabs)],0*new100[,c(xlabs,"y")])}
  
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
  phat<-mean(df$a)
  
  #estimated values of mu1 and mu0
  u1hat<-1/N*sum((df$ys*df$a-phat*p101hat)/(phat*(p111hat-p101hat)))
  u0hat<-1/N*sum((df$ys*(1-df$a)-(1-phat)*p100hat)/((1-phat)*(p110hat-p100hat)))
  
  #SSW estimator 
  finalest<-u1hat-u0hat
  
  #just point SSW estimator is printed
  if(varprint==F){
    return(c("Estimator"=finalest))
  }else{
    
    ##for asymptotic variance estimation
    
    #fitted values
    pmehat<-fitted.values(mefit)
    
    #number of parameters 
    nparams<-ncol(dfvcov)
    
    #cluster robust sandwich variance 
    if(crob==T){
    
    #meat matrix 
    meatmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    #bread matrix 
    #number of clusters
    m<-length(unique(df$Id))
    
    #unbiased estimating equations 
    for(i in unique(df$Id)){
      
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
      
      #mic for u1hat 
      mi6x<-sum((dfi$ys*dfi$a-phat*p101hat[clustidfull])/(phat*(p111hat[clustidfull]-p101hat[clustidfull]))-u1hat)
      
      #mic for u0hat
      mi7x<-sum((dfi$ys*(1-dfi$a)-(1-phat)*p100hat[clustidfull])/((1-phat)*(p110hat[clustidfull]-p100hat[clustidfull]))-u0hat)
      meatmat1<-meatmat1+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
      
    }
    }else{
      #individual sandwich variance
      #if individual index above by ij
      meatmat2<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
      for(i in row.names(df)){
        dfi<-df[i,]
        #logistic regression on validation
        if(dfi$v==1){
          idv<-which(row.names(dfv)==i)
          pmehati<-pmehat[idv]
          dfvicov<-dfvcov[idv,]
          ysi<-dfv[idv,]$ys
          estfuni<-t(t(dfvicov)%*%(ysi-pmehati))}else{
            estfuni<-rep(0,nparams)} #this is the estimating function by person, the phi_ij, summing over ij

        #remaining components
        mi5x<-dfi$a-phat
        mi6x<-(dfi$ys*dfi$a-phat*p101hat[i])/(phat*(p111hat[i]-p101hat[i]))-u1hat
        mi7x<-(dfi$ys*(1-dfi$a)-(1-phat)*p100hat[i])/((1-phat)*(p110hat[i]-p100hat[i]))-u0hat
        meatmat2<-meatmat2+c(estfuni,mi5x,mi6x,mi7x)%*%t(c(estfuni,mi5x,mi6x,mi7x))
      }}
    
    ##filling components of bread matrix, sum of jacobian matrices
    breadmat1<-matrix(rep(0,(nparams+3)^2),ncol=nparams+3)
    
    #theta theta
    breadmat1[1:nparams,1:nparams]<--t(dfvcov*pmehat*(1-pmehat))%*%as.matrix(dfvcov)
    
    #pi pi
    breadmat1[nparams+1,nparams+1]<--N
    
    #mu1 theta
    breadmat1[nparams+2,1:nparams]<-apply((df101cov*phat^2*(-p101hat)*(1-p101hat)*(p111hat-p101hat)-(df$a*df$ys-phat*p101hat)*phat*(df111cov*p111hat*(1-p111hat)-df101cov*p101hat*(1-p101hat)))/
                                            (phat^2*(p111hat-p101hat)^2),2,sum)
    
    #mu1 pi
    breadmat1[nparams+2,nparams+1]<-sum(-df$a*df$ys/(phat^2*(p111hat-p101hat)))
    
    #mu1 mu1
    breadmat1[nparams+2,nparams+2]<--N
    
    #mu0 theta
    breadmat1[nparams+3,1:nparams]<-apply((df100cov*(1-phat)^2*(-p100hat)*(1-p100hat)*(p110hat-p100hat)-((1-df$a)*df$ys-(1-phat)*p100hat)*(1-phat)*(df110cov*p110hat*(1-p110hat)-df100cov*p100hat*(1-p100hat)))/
                                            ((1-phat)^2*(p110hat-p100hat)^2),2,sum)
    
    #mu0 pi
    breadmat1[nparams+3,nparams+1]<-sum(((1-df$a)*df$ys)/((1-phat)^2*(p110hat-p100hat)))
    
    #mu0 mu0
    breadmat1[nparams+3,nparams+3]<--N
    
    #delta method
    if(crob==T){
    #for cluster robust
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat1%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)}
    else{
    #for individual sand
    varest<-t(c(rep(0,nparams+1),1,-1))%*%solve(breadmat1)%*%meatmat2%*%solve(t(breadmat1))%*%c(rep(0,nparams+1),1,-1)
    }
    
    #if no small sample correction
    if(corcl==F){
      lb<-finalest+qnorm(alpha/2)*sqrt(varest)
      ub<-finalest+qnorm(1-alpha/2)*sqrt(varest)
      output<-c(finalest,varest,lb,ub)
      if(crob==T){
      names(output)<-c("Estimator","VarAsymCL","LBCL","UBCL")}
      else{
        names(output)<-c("Estimator","VarAsymInd","LBInd","UBInd")
      }
    }
    
    #if small sample correction
    if(corcl==T){
      lb1T<-finalest+qt(alpha/2,df=m-7)*sqrt(varest)
      ub1T<-finalest+qt(1-alpha/2,df=m-7)*sqrt(varest)
      output<-c(finalest,varest,lb1T,ub1T)
      names(output)<-c("Estimator","VarAsymCL","LBCLCorT","UBCLCorT")
    }

    }
    return(output)}

#return outputs
#warning if ATE does not respect bounds
if(varprint==F){
  output<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=varprint,crob=crob,corcl=corcl)
  if(abs(output[1])>1){
    warning("Estimator outside valid bounds of -1 to 1 for ATE.")
  }
  
  return(output)
}

#if asymptotic variance 
if(varprint==T & boot==F){
  var_output<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=varprint,crob=crob,corcl=corcl)
  if(abs(var_output[1])>1){
    warning("Estimator outside valid bounds of -1 to 1 for ATE.")
  }
  
  return(var_output)
}

#if bootstrap methods  
if(varprint==T & boot==T){

#non-par bootstrap function
boot_mesb_estim_sub<-function(dfog,iters,diffx=diffx,xlabs,clabs,crob){
  estboot<-rep(0,iters)
  failboot<-rep(0,iters)
  for(b in 1:iters){
    #resample with replacement from clusters
    if(crob==T){
      idval<-length(unique(dfog$Id))
      clustidv<-sample(unique(dfog$Id),idval,replace=T)
      dflist<-vector("list",idval)
      for(j in 1:idval){
        dflist[[j]]<-dfog[dfog$Id==clustidv[j],]
      }
      df<-do.call(rbind,dflist)
    }
    #resamples with replacement from individuals
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
      failboot[b]<-1
    }}
  
  #bootstrap estimate of variance
  varboot<-var(estboot,na.rm=T)
  
  #percentile interval
  bounds<-quantile(estboot,c(alpha/2,1-alpha/2),na.rm=T)
  lb<-bounds[1]
  ub<-bounds[2]
  output<-c(varboot,lb,ub,mean(failboot))
  
  #cluster boot
  if(crob==T){
    names(output)<-c("VarBootCl","LBBootCL","UBBootCL","Percent Invalid")
  }
  
  #individual boot
  if(crob==F){
    names(output)<-c("VarBootInd","LBBootInd","UBBootInd","Percent Invalid")
  }
  
  return(output)}

finalest<-mesb_estim_sub(df=df,diffx=diffx,xlabs=xlabs,clabs=clabs,varprint=F)
varvec<-boot_mesb_estim_sub(dfog=df,iters=iters,diffx=diffx,xlabs=xlabs,clabs=clabs,crob=crob)
output<-c(finalest,varvec)
if(abs(output[1])>1){
  warning("Estimator outside valid bounds of -1 to 1 for ATE.")
}

return(output)

}}