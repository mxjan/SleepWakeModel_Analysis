GetTime_TS<-function(expr_train,time_train,expr_pred,time_pred,doplot=F,a=0.5,s=exp(-1.42),fittedOsciModel=NULL,cols="black"){
  
  # Train
  TSorig <- trainTimeStamp(
    expr=t(expr_train), 
    subjIDs=rep("1",nrow(expr_train)) ,
    times=time_train,
    trainFrac=1, 
    recalib=FALSE, 
    a=a, s=s, 
    plot=F 
  )
  
  #plot(TSorig$cv.fit)
  
  # new-points
  newx<-as.matrix(expr_pred)
  
  plot(time_pred,
       XY2dectime(predict(TSorig$cv.fit,newx,s=s)[,,1]),xaxt="n",yaxt="n",ylim=c(0,24),xlab="Sampled Time",ylab="Predicted ZT",col=cols,pch=19)
  rect(xleft = 48,xright = 54,ybottom = -5,ytop = 28,col = rgb(1,0,0,.3))
  #rect(xleft = 0,xright = 120,ybottom = 12,ytop = 24,col = rgb(0,0,0,.3))
  axis(1,unique(time_pred),paste(unique(time_pred),"\n(ZT",unique(time_pred) %% 24,")",sep=""),padj = 0.25)
  axis(2,seq(0,24,by=6),seq(0,24,by=6) %% 24,las=2)
  abline(h=seq(0,24,by=6),col="grey75",lty=2)
  
  if (! is.null(fittedOsciModel)){
    newx<-t(fittedOsciModel)[,colnames(expr_train)]
    time<-as.numeric(colnames(fittedOsciModel))
    preds<-XY2dectime(predict(TSorig$cv.fit,newx,s=s)[,,1])#+rep(cumsum(rep(24,11))-48,each=240)
    
    # 
    dd<-diff(preds)
    dd[dd< -12]<- 24 + dd[dd< -12] 
    preds2<-cumsum(dd)
    # lines(time,rep(preds[time >= 24 & time<48],11),col="black",lty=2)
    lines(time[-1],preds2,col="black")
    lines(time[-1],preds2-24,col="black")
    lines(time[-1],preds2-48,col="black")
    lines(time[-1],preds2-72,col="black")
    lines(time[-1],preds2-96,col="black")
    
    # bsl
    lines(time[-1]+24,preds2,col="black",lty=2)
    lines(time[-1]+48,preds2,col="black",lty=2)
    lines(time[-1]+72,preds2,col="black",lty=2)
    lines(time[-1]+96,preds2,col="black",lty=2)
    
  }
  newx<-as.matrix(expr_pred)
  return(list(TSorig=TSorig,preds=list(Time=time_pred,preds=XY2dectime(predict(TSorig$cv.fit,newx,s=s)[,,1])),
              PredsLineTime=time[-1],PredsLine=preds2))
  
}


PlotTS<-function(TSres,tissue="Cortex"){
  
  dd<-cbind.data.frame(Time=TSres$preds$Time,Preds=TSres$preds$preds)
  gg<-ggplot(aes(x=Time,y=Preds),data=dd)+annotate("rect",xmin=48,xmax=54,ymin=-Inf,ymax=Inf,fill=rgb(1,0,0,.25),color="black")
  gg<-gg+scale_x_continuous(breaks=seq(24,102,by=12),labels=c("24\n(ZT0)","36\n(ZT12)","48\n(ZT0)",
                                                              "60\n(ZT12)","72\n(ZT0)","84\n(ZT12)","96\n(ZT0)"))+
    scale_y_continuous(breaks=seq(0,24,by=6),labels=c(0,6,12,18,0))
  gg<-gg+geom_point(color=colorCode[[tissue]])
  
  gg<-gg + annotate("path",x=TSres$PredsLineTime,TSres$PredsLine)
  gg<-gg + annotate("path",x=TSres$PredsLineTime,TSres$PredsLine-24)
  gg<-gg + annotate("path",x=TSres$PredsLineTime,TSres$PredsLine-48)
  gg<-gg + annotate("path",x=TSres$PredsLineTime,TSres$PredsLine-72)
  gg<-gg + annotate("path",x=TSres$PredsLineTime,TSres$PredsLine-96)

  gg<-gg + annotate("path",x=TSres$PredsLineTime+24,TSres$PredsLine,linetype="22")
  gg<-gg + annotate("path",x=TSres$PredsLineTime+48,TSres$PredsLine,linetype="22")
  gg<-gg + annotate("path",x=TSres$PredsLineTime+72,TSres$PredsLine,linetype="22")
  gg<-gg + annotate("path",x=TSres$PredsLineTime+96,TSres$PredsLine,linetype="22")
  
  gg<-gg+coord_cartesian(xlim = c(24,102),ylim = c(0,24))
  
  gg<-gg+theme_bw()+theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
  
  gg<-gg+ylab("Corresponding expression\n in baseline as ZT")
  
  gg
  
}

GetPenalization_LeaveOneOut<-function(expr_train,time_train){
  
  # results
  res<-matrix(nrow=0,ncol=3);colnames(res)<-c("a","s","rss")
  
  # sequence of tested alpha
  aseq<-seq(0,1,by=.1)#seq(0.4,0.6,by=.05)
  # tested lambda
  sseq<-exp(seq(-3,-1,by=.5))
  # Samples
  samples<-1:nrow(expr_train)
  
  # Run
  for (a in aseq){
    for (s in sseq){
      
      RSSloo<-c() # store residuals
      for (i in samples){
        
        TSorig <- trainTimeStamp(
          expr=t(expr_train[samples[-i],]), 
          subjIDs=rep("1",nrow(expr_train[samples[-i],])) ,
          times=time_train[samples[-i]],
          trainFrac=1, 
          recalib=FALSE, 
          a=a, s=s,
          plot=F,grouped=FALSE
        )
        preds<- XY2dectime(t(as.matrix(predict(TSorig$cv.fit,as.matrix(expr_train)[samples[i],,drop=F],s=s)[,,1])))
        RSS<-sum(apply(abs(cbind(preds - time_train[i] %% 24,preds - time_train[i] %% 24 -24)),1,min)^2)
        RSSloo<-c(RSSloo,RSS)
      }
      
      print(paste(a,s,sum(RSSloo)))
      res<-rbind(res,c(a,s,sum(RSSloo)))
      
    }
  }
  
  print(res[which.min(res[,3]),])
  return(list("a"=res[which.min(res[,3]),1],"s"=res[which.min(res[,3]),2]))
  
}

TScoef <- function(timestamp,s=timestamp$cv.fit$lambda.min){
  out <- as.matrix(do.call(cbind,coef(timestamp$cv.fit, s=s)))
  out <- out[rowSums(out!=0)>0,][-1,]
  colnames(out) <- c("sunX","sunY")
  return(out)
}

Train4FCV<-function(expr_train,time_train,seed=1){
  
  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  mm<-matrix(nrow=0,ncol=3)
  
  for (a in seq(0,1,by=.1)){
    for (s in exp(seq(-8,-2,by=.5))){
      RSSrands<-c()
      for (j in seq(1,100)){
        set.seed(j)
        sampl<-chunk2(sample(1:nrow(expr_train)),4)
        sampl_rss<-sampl
        for (i in names(sampl)){
          train<-unlist(sampl[names(sampl)[! names(sampl) %in% i]])
          pred<-sampl[[names(sampl)[names(sampl) %in% i]]]
          # Train
          TSorig <- trainTimeStamp(
            expr=t(expr_train[train,]), 
            subjIDs=rep("1",nrow(expr_train[train,])) ,
            times=time_train[train],
            trainFrac=1, 
            recalib=FALSE, 
            a=a, s=s, 
            plot=F 
          )
          preds<- XY2dectime(predict(TSorig$cv.fit,as.matrix(expr_train)[pred,],s=s)[,,1])
          RSS<-sum(apply(abs(cbind(preds - time_train[pred] %% 24,preds - time_train[pred] %% 24 -24)),1,min)^2)
          sampl_rss[[i]]<-RSS
        }
        
        RSSrands<-c(RSSrands,sum(unlist(sampl_rss)))
        
      }
      
      print(paste(a,s,mean(RSSrands)))
      mm<-rbind(mm,c(a,s,mean(RSSrands)))
    }
  }
  print(mm[which.min(mm[,3]),])
  return(list("a"=mm[which.min(mm[,3]),1],"s"=mm[which.min(mm[,3]),2]))
}
