context("testRMST")

# Simulate data to test
set.seed(12345)
samp <- function(u, b, shape) { -(exp(b)^(-1))^(1 / shape) * (log(1 - u))^(1 / shape) }
wei  <- function(x) (exp(b1)*shape*x^(shape-1))*exp(-exp(b1)*x^shape)
weib <- function(x) (1/wei(0))*wei(x) # scale the weibull curve so that it's 1 at time 0

median <- 8; shape <- 1; HR <- 0.7

b1     <- log((-1 / (median^shape)) * log(0.5))
b2     <- -log(1 / HR)

surv <-  c(vapply(runif(500), samp, b = b1, shape = shape, FUN.VALUE = numeric(1)),
           vapply(runif(500), samp, b = b1 + b2, shape = shape, FUN.VALUE = numeric(1)))

arm   <- c(rep(0,500),rep(1,500))
event <- rep(1,length(surv))
dat   <- data.frame("surv" = surv,"event" = event,"arm" = arm)

test_that("UseCaseKM",{
  out1 <- rmst2(dat$surv,dat$event,dat$arm, tau=10)
  out2 <- rmst(Surv(surv, event) ~ arm, data=dat, trunc=10, alpha=0.05, method="KM")

  expect_equal(unname(out1$RMST.arm0$rmst[1]),as.numeric(as.character(out2[[1]][1,1]))) # RMST in arm 0
  expect_equal(unname(out1$RMST.arm1$rmst[1]),as.numeric(as.character(out2[[1]][2,1]))) # RMST in arm 1
  expect_equal(round((out1$unadjusted.result[1,3]-out1$unadjusted.result[1,1])/1.96,4),
               round(as.numeric(out2$diff[,"se"]),4))                                  # SE of difference

  out6 <- rmst(Surv(surv, event) ~ arm, data=dat[dat$arm==0,], trunc=10, alpha=0.05, method="KM")
  expect_equal(as.numeric(as.character(out6[[1]][1,1])),7.015411) # RMST works in single arm
})

test_that("UseCasePseudo",{
  pseudo <- pseudomean(dat$surv,dat$event,10)
  a      <- cbind(dat,"pseudo" = pseudo,"id" = 1:nrow(dat))
  fit    <- geese(pseudo ~ arm,
                  data = a, id=id, jack = TRUE, family=gaussian,
                  corstr="independence", scale.fix=FALSE)
  sfit   <- summary(fit)

  out3 <- rmst(Surv(surv, event) ~ arm, data=dat, trunc=10, alpha=0.05, method="pseudo")

  expect_equal(unname(fit$beta[1]),as.numeric(out3[[1]][1,1]))             # RMST in arm 0
  expect_equal(unname(fit$beta[1]+fit$beta[2]),as.numeric(out3[[1]][2,1])) # RMST in arm 1
  expect_equal(sfit$mean[2,2],as.numeric(out3[[2]][2]))                    # SE of difference

  out7 <- rmst(Surv(surv, event) ~ arm, data=dat[dat$arm==0,], trunc=10, alpha=0.05, method="pseudo")
  expect_equal(as.numeric(as.character(out7[[1]][1,1])),7.015411) # RMST works in single arm
})

test_that("UseCaseRP",{
  out4 <- flexsurvspline(Surv(surv, event) ~ arm, data=dat, k=2)
  surv1 <- function(x,arm) unlist(summary(out4, newdata=data.frame("arm"=arm),t=x)[[1]]["est"])

  out5 <- rmst(Surv(surv, event) ~ arm, data=dat, trunc=10, alpha=0.05, method="RP", knots=2)

  var       <- c(out4$res.t[out4$dlist$pars,"est"], out4$res[out4$covpars,"est"])
  modelvcov <- vcov(out4)
  slength   <- length(out4$res.t[out4$dlist$pars,"est"])

  surv2 <- function(x,var,dframe){
    unlist(modflexsummary(var,slength,out4,t=x,ci=FALSE,newdata=dframe)[[1]]["est"])
  }

  diffunc <- function(var, end) {abs(integrate(surv2,0,end,  dframe=data.frame(arm=0),var=var)$value -
                                       integrate(surv2,0,end,dframe=data.frame(arm=1),var=var)$value )}

  SEdif <- function(var, end){
    gd <- grad(diffunc, var, end = end)
    sqrt(gd %*% modelvcov %*% gd)
  }

  expect_equal(integrate(surv1,0,10,arm=0)$value,as.numeric(out5[[1]][1,1])) # RMST in arm 0
  expect_equal(integrate(surv1,0,10,arm=1)$value,as.numeric(out5[[1]][2,1])) # RMST in arm 1
  expect_equal(as.numeric(SEdif(var,10)), as.numeric(out5[[2]][2]))         # SE of difference

  out8 <- rmst(Surv(surv, event) ~ arm, data=dat[dat$arm==0,], trunc=10, alpha=0.05, method="RP")
  expect_equal(round(as.numeric(as.character(out8[[1]][1,"RMST"])),5),7.01454) # RMST works in single arm
  out9 <- rmst(Surv(surv, event) ~ 1, data=dat, trunc=10, alpha=0.05, method="RP")
  expect_equal(round(as.numeric(as.character(out9[[1]][1,"RMST"])),6),7.298897) # RMST works in single arm
})

test_that("RPWeibull",{
  model <- flexsurvreg(Surv(surv, event) ~ arm, data=dat,dist="weibull")
  surv <- function(x) unlist(summary(model,newdata=data.frame("arm"=0),t=x,ci=F)[[1]]["est"])

  expect_equal(
  rmst(Surv(surv, event) ~ arm, data=dat, trunc=5, alpha=0.05, method="RP", knots=0)[[1]][1,1],
  integrate(surv,0,5)$value,
  tolerance=1e-6
  )
})

test_that("KMequalpseudo",{
  expect_equal(
    rmst(Surv(surv, event) ~ arm, data=dat, trunc=10, alpha=0.05, method="KM")[[2]][1],
    as.numeric(rmst(Surv(surv, event) ~ arm, data=dat, trunc=10, alpha=0.05, method="pseudo")[[2]][1])
  )
})

test_that("ModifiedSummary",{
  testmodel <- flexsurvspline(Surv(surv,event) ~ arm, data=dat, k=2)

  var       <- c(testmodel$res.t[testmodel$dlist$pars,"est"], testmodel$res[testmodel$covpars,"est"])
  modelvcov <- vcov(testmodel)
  slength   <- length(testmodel$res.t[testmodel$dlist$pars,"est"])

  expect_equal(
    modflexsummary(var,slength,testmodel,t=10,ci=FALSE,newdata=data.frame("arm"=1))[[1]]["est"],
    summary(testmodel,t=10,ci=FALSE,newdata=data.frame("arm"=1))[[1]]["est"]
  )
})


