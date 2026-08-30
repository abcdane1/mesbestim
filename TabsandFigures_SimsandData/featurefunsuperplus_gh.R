#function to obtain additional features regarding the selection subset and classification

#not reported in manuscript (expect to obtain proportion of Y=1)

##empirical selection and mis-measurement proportions
feature_fun_superplus<-function(df){
  #selection proportion
  select_prob<-mean(df$v)
  dfv<-df[which(df$v==1),]
  #overall mis-classification proportion
  overall_missclass<-mean(dfv$y!=dfv$ys)
  #prevalence
  prev<-mean(dfv$y)
  prevys<-mean(dfv$ys)
  #proportion of y=1 and y*=0
  y1ys0_missclass<-mean(dfv$y==1 & dfv$ys==0)
  #proportion of y=0 and y*=1
  y0ys1_missclass<-mean(dfv$y==0 & dfv$ys==1)
  #y*=1|y=0
  ys1givy0<-y0ys1_missclass/(1-prev)
  #y*=0|y=1
  ys0givy1<-y1ys0_missclass/prev
  #y=0|y*=1
  y0givys1<-y0ys1_missclass/prevys
  #y=1|y*=0
  y1givys0<-y1ys0_missclass/(1-prevys)
  
  #trt1select
  trt1select_prop<-mean(dfv$a)
  
  #treatment 1
  dfv1<-dfv[which(dfv$a==1),]
  #proportion of y=1 and y*=0
  trt1y1ys0_missclass<-mean(dfv1$y==1 & dfv1$ys==0)
  #proportion of y=0 and y*=1
  trt1y0ys1_missclass<-mean(dfv1$y==0 & dfv1$ys==1)
  #prevalence 
  prev1<-mean(dfv1$y)
  prevys1<-mean(dfv1$ys)
  
  #y*=1|y=0
  trt1ys1givy0<-trt1y0ys1_missclass/(1-prev1)
  #y*=0|y=1
  trt1ys0givy1<-trt1y1ys0_missclass/prev1
  #y=0|y*=1
  trt1y0givys1<-trt1y0ys1_missclass/prevys1
  #y=1|y*=0
  trt1y1givys0<-trt1y1ys0_missclass/(1-prevys1)
  
  
  #treatment 0
  dfv0<-dfv[which(dfv$a==0),]
  #proportion of y=1 and y*=0
  trt0y1ys0_missclass<-mean(dfv0$y==1 & dfv0$ys==0)
  #proportion of y=0 and y*=1
  trt0y0ys1_missclass<-mean(dfv0$y==0 & dfv0$ys==1)
  #prevalence
  prev0<-mean(dfv0$y)
  prevys0<-mean(dfv0$ys)
  
  #y*=1|y=0
  trt0ys1givy0<-trt0y0ys1_missclass/(1-prev0)
  #y*=0|y=1
  trt0ys0givy1<-trt0y1ys0_missclass/prev0
  #y=0|y*=1
  trt0y0givys1<-trt0y0ys1_missclass/prevys0
  #y=1|y*=0
  trt0y1givys0<-trt0y1ys0_missclass/(1-prevys0)
  
  output<-c(select_prob,trt1select_prop,overall_missclass,y1ys0_missclass,y0ys1_missclass,
            trt1y1ys0_missclass,trt1y0ys1_missclass,trt0y1ys0_missclass,trt0y0ys1_missclass,prev,prevys,prev1,prevys1,prev0,prevys0,
            ys1givy0,ys0givy1,y0givys1,y1givys0,trt1ys1givy0,trt1ys0givy1,trt1y0givys1,trt1y1givys0,
            trt0ys1givy0,trt0ys0givy1,trt0y0givys1,trt0y1givys0
  )
  
  names(output)<-c("SelectProp","Trt1PropGivSelect",
                   "OverallMissclass","Missclass_True1Meas0","Missclass_True0Meas1",
                   "Trt1Missclass_True1Meas0","Trt1Missclass_True0Meas1",
                   "Trt0Missclass_True1Meas0","Trt0Missclass_True0Meas1",
                   "Sel_Prev","Sel_Prevys","Sel_PrevTrt1","Sel_PrevysTrt1","Sel_PrevTrt0","Sel_PrevysTrt0",
                   "ys1givy0","ys0givy1","y0givys1","y1givys0","trt1ys1givy0","trt1ys0givy1","trt1y0givys1","trt1y1givys0",
                   "trt0ys1givy0","trt0ys0givy1","trt0y0givys1","trt0y1givys0")
  return(output)
}