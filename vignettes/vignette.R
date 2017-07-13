## ----fig1, echo=FALSE, message=FALSE-------------------------------------
library(RMST)
D <- rmst2.sample.data()
plot.rmst2=function(x, xlab="", ylab="", col="black", col.RMST="gray80", col.RMTL="orange",density=80, angle=85,...){
  
  if(is.null(x$unadjusted.result)) stop("Please run rmst2 without covariates")
  
  if(!is.null(x$unadjusted.result)){
    
    ft1=x$RMST.arm1$fit
    tau=x$tau
    
    par(mfrow=c(1,1))
    
    #=== arm 1 ===
    fit=ft1
    
    tmp.xx=c(0, fit$time); tmp.yy=c(1, fit$surv) ;
    idx=tmp.xx<=tau
    y.tau = min(tmp.yy[idx])
    xx=c(tmp.xx[idx],   tau)
    yy=c(tmp.yy[idx], y.tau)  
    x.step=sort(c(0, xx, xx))
    y.step=rev(sort(c(1,1,yy, yy[-length(yy)])))
    
    #--------
    plot(fit, mark.time=F, conf.int=F, lwd=2, main="", xlab=xlab, ylab=ylab, col=col)
    
    for (i in seq(1, length(x.step), by=2)){  
      polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(0, 0, y.step[i+1], y.step[i]), col= col.RMST, density=density, angle=angle, lwd=5)
    }

    
    x.step=sort(c(0, tmp.xx, tmp.xx))
    y.step=rev(sort(c(1,1,tmp.yy, tmp.yy[-length(tmp.yy)])))
    lines(x.step, y.step, col=col, lwd=3) 
    # text(5,0.4, paste(round(rmst$rmst[1], digits=2),"years"), cex=0.9)
    

    
  }
  
  
}

km <- rmst2(D$time, D$status, D$arm,tau=5)
plot(km)

## ---- message=FALSE------------------------------------------------------
# Set the number of digits to be shown in outputs to be 3
options(digits = 3)

## ------------------------------------------------------------------------
fit <- rmst(Surv(time,status) ~ arm, data=D, trunc=5, alpha=0.05, method="KM")
fit[[1]]

## ------------------------------------------------------------------------
fit[[2]]

## ---- message=FALSE------------------------------------------------------
fit <- rmst(Surv(time,status) ~ arm, data=D, trunc=5, alpha=0.05, method="pseudo")
fit[[1]]

## ---- message=FALSE------------------------------------------------------
fit[[2]]

## ---- message=FALSE------------------------------------------------------
fit <- rmst(Surv(time,status) ~ arm, data=D, trunc=5, alpha=0.05, method="RP", knots=2)
fit[[1]]

## ---- message=FALSE------------------------------------------------------
fit[[2]]

