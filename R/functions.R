##' Modified flexsurv summary function for use with the delta method
##'
##' This a modified flexsurv summary function which takes the coefficients of the model, and the number of ancillary parameters as the first arguments of function to facilitate compatibility with numDeriv. It is intended primarily for internal useage. Please use the summary.flexsurv function unless computation of a standard error is required via the delta method.
##' @param var Coefficients from a flexsurvspline model
##' @param sl The number of ancillary (gamma) covariates
##' @param object Output from flexsurvspline, representing a fitted survival model object.
##' @param newdata Data frame containing covariate values to produce fitted values for. Or a list that can be coerced to such a data frame. There must be a column for every covariate in the model formula, and one row for every combination of covariates the fitted values are wanted for. These are in the same format as the original data, with factors as a single variable, not 0/1 contrasts. \cr\cr If this is omitted, if there are any continuous covariates, then a single summary is provided with all covariates set to their mean values in the data - for categorical covariates, the means of the 0/1 indicator variables are taken. If there are only factor covariates in the model, then all distinct groups are used by default.
##' @param X Alternative way of defining covariate values to produce fitted values for. Since version 0.4, newdata is an easier way that doesn't require the user to create factor contrasts, but X has been kept for backwards compatibility. \cr\cr Columns of X represent different covariates, and rows represent multiple combinations of covariate values. For example matrix(c(1,2),nrow=2) if there is only one covariate in the model, and we want survival for covariate values of 1 and 2. A vector can also be supplied if just one combination of covariates is needed. \cr\cr For "factor" (categorical) covariates, the values of the contrasts representing factor levels (as returned by the contrasts function) should be used. For example, for a covariate agegroup specified as an unordered factor with levels 20-29, 30-39, 40-49, 50-59, and baseline level 20-29, there are three contrasts. To return summaries for groups 20-29 and 40-49, supply X = rbind(c(0,0,0), c(0,1,0)), since all contrasts are zero for the baseline level, and the second contrast is "turned on" for the third level 40-49.
##' @param type "survival" for survival probabilities. \cr "cumhaz" for cumulative hazards. \cr "hazard" for hazards. \cr Ignored if "fn" is specified.
##' @param fn Custom function of the parameters to summarise against time. This has optional first two arguments t representing time, and start representing left-truncation points, and any remaining arguments must be parameters of the distribution. It should return a vector of the same length as t.
##' @param t Times to calculate fitted values for. By default, these are the sorted unique observation (including censoring) times in the data - for left-truncated datasets these are the "stop" times.
##' @param start Optional left-truncation time or times. The returned survival, hazard or cumulative hazard will be conditioned on survival up to this time. \cr\cr A vector of the same length as t can be supplied to allow different truncation times for each prediction time, though this doesn't make sense in the usual case where this function is used to calculate a predicted trajectory for a single individual. This is why the default start time was changed for version 0.4 of flexsurv - this was previously a vector of the start times observed in the data.
##' @param ci Set to FALSE to omit confidence intervals.
##' @param B Number of simulations from the normal asymptotic distribution of the estimates used to calculate confidence intervals. Decrease for greater speed at the expense of accuracy, or set B=0 to turn off calculation of CIs.
##' @param cl Width of symmetric confidence intervals, relative to 1.
##' @param tidy If TRUE, then the results are returned as a tidy data frame instead of a list. This can help with using the ggplot2 package to compare summaries for different covariate values.
##' @param ... Further arguments passed to or from other methods. Currently unused.
##' @return If tidy=FALSE, a list with one component for each unique covariate value (if there are only categorical covariates) or one component (if there are no covariates or any continuous covariates). Each of these components is a matrix with one row for each time in t, giving the estimated survival (or cumulative hazard, or hazard) and 95\% confidence limits. These list components are named with the covariate names and values which define them. \cr\cr If tidy=TRUE, a data frame is returned instead. This is formed by stacking the above list components, with additional columns to identify the covariate values that each block corresponds to. \cr\cr If there are multiple summaries, an additional list component named X contains a matrix with the exact values of contrasts (dummy covariates) defining each summary. \cr\cr The plot.flexsurvreg function can be used to quickly plot these model-based summaries against empirical summaries such as Kaplan-Meier curves, to diagnose model fit. \cr\cr Confidence intervals are obtained by random sampling from the asymptotic normal distribution of the maximum likelihood estimates (see, e.g. Mandel (2013)).
##' @export
modflexsummary <- function(var, sl, object, newdata=NULL, X=NULL, type="survival", fn=NULL,
                           t=NULL, start=0, ci=TRUE, B=1000, cl=0.95, tidy=FALSE,
                           ...)
{
                                  x <- object
  dat <- x$data
  Xraw <- model.frame(x)[,unique(attr(model.frame(x),"covnames.orig")),drop=FALSE]
  isfac <- sapply(Xraw,is.factor)
  type <- match.arg(type, c("survival","cumhaz","hazard"))
  if (is.null(newdata)){
    if (is.vector(X)) X <- matrix(X, nrow=1)
    if (x$ncovs > 0 && is.null(X)) {
      ## if any continuous covariates, calculate fitted survival for "average" covariate value
      if (!all(isfac))
        X <- matrix(colMeans(model.matrix(x)) ,nrow=1)
      ## else calculate for all different factor groupings
      else {
        X <- unique(model.matrix(x))
        ## build names like "COVA=value1,COVB=value2"
        nam <- as.matrix(unique(Xraw))
        for (i in 1:ncol(nam)) nam[,i] <- paste(colnames(nam)[i], nam[,i], sep="=")
        rownames(X) <- apply(nam, 1, paste, collapse=",")
      }
    }
    else if (is.null(X)) X <- as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1))
    else if (!is.matrix(X) || (is.matrix(X) && ncol(X) != x$ncoveffs)) {
      plural <- if (x$ncoveffs > 1) "s" else ""
      stop("expected X to be a matrix with ", x$ncoveffs, " column", plural, " or a vector with ", x$ncoveffs, " element", plural)
    }
  } else
    
    X <- form.model.matrix(object, as.data.frame(newdata))
  if (is.null(t))
    t <- sort(unique(dat$Y[,"stop"]))
  if (length(start)==1)
    start <- rep(start, length(t))
  else if (length(start) != length(t))
    stop("length of \"start\" is ",length(start),". Should be 1, or length of \"t\" which is ",length(t))

  if (is.null(fn)) {
    
    fn <- summary.fns(x, type)
  }
  
  fn <- expand.summfn.args(fn)
  fncall <- list(t,start)
  beta <- if (x$ncovs==0) 0 else var[(sl+1):length(var)]
  if (ncol(X) != length(beta)){
    ## don't think we should ever reach here - error should be caught in newdata or X
    isare <- if(length(beta)==1) "is" else "are"
    plural <- if(ncol(X)==1) "" else "s"
    pluralc <- if(length(beta)==1) "" else "s"
    stop("Supplied X has ", ncol(X), " column",plural," but there ",isare," ",
         length(beta), " covariate effect", pluralc)
  }
  dlist <- x$dlist
  ret <- vector(nrow(X), mode="list")
  if(!is.null(newdata)){
    nd <- attr(X, "newdata")
    covnames <- apply(as.data.frame(nd), 1, function(x)paste0(names(nd), "=", x, collapse=", "))
  }
  else covnames <- rownames(X)
  names(ret) <- covnames
  for (i in 1:nrow(X)) {
    
    basepars.mat <- add.covs(x, var[1:sl], beta, X[i,,drop=FALSE], transform=FALSE)
    basepars <- as.list(as.data.frame(basepars.mat))
    fncall[dlist$pars] <- basepars
    for (j in seq_along(x$aux)){
      fncall[[names(x$aux)[j]]] <- x$aux[[j]]
    }
    y <- do.call(fn, fncall)
    if (ci){
      
      res.ci <- cisumm.flexsurvreg(x, t, start, X[i,,drop=FALSE], fn=fn, B=B, cl=cl)
      ly <- res.ci[,1]
      uy <-  res.ci[,2]
    }
    ret[[i]] <- data.frame(time=t, est=y, row.names=NULL)
    if (ci) { ret[[i]]$lcl <- ly; ret[[i]]$ucl <- uy}
  }
  if (x$ncovs>0) attr(ret,"X") <- X
  if (tidy) {
    ret <- do.call("rbind", ret)
    covdf <- unique(Xraw)[rep(seq_len(nrow(unique(Xraw))), each=length(t)), , drop=FALSE]
    rownames(ret) <- NULL
    ret <- cbind(ret, covdf)
  }
  class(ret) <- c("summary.flexsurvreg",class(ret))
  ret
}

#########################################################################################################
#########################################################################################################

##' Calculate the RMST in each arm and the between-group difference
##'
##' This function calculates an estimate of the RMST using the specified method
##' @return The Restricted Mean Survival Time computated
##' @param formula A formula object, with the response on the left of a ~ operator, and the a single term for arm on the right. The response must be a survival object as returned by the Surv function.
##' @param data A data.frame in which to interpret the variables named in the formula.
##' @param trunc The time point to truncated the estimate of the RMST at.
##' @param alpha The default is 0.05. (1-alpha) confidence intervals are reported.
##' @param method The method to estimate the RMST. Note that 2 knots will be used by default when using the Royston-Parmar method unless otherwise specified.
##' @param ... Additional arguments to be passed to the method
##'
##' @details
##' The RMST is estimate which can be found by calculating the area under the survival curve. This function facilitates the implementation the three most discussed approaches in current literature. For further details on how the RMST is calculated for each method see the document for the appropriate method, i.e. \code{?rmstKM}, \code{?rmstPseudo}, or \code{?rmstRP}.
##' @examples
##' # Load in an example data set
##' library(survRM2)
##' D <- rmst2.sample.data()
##'
##' # Calculate the RMST using pseudovalues
##' rmst(Surv(time, status) ~ arm, data=D, trunc=5, alpha=0.05, method="pseudo")
##' @export
rmst <- function(formula, data, trunc, alpha, method=c("KM", "pseudo", "RP"), ...) {
  Call <- match.call()

  result <- switch(match.arg(method),
                   KM=rmstKM(formula, data, trunc, alpha, ...),
                   pseudo=rmstPseudo(formula, data, trunc, alpha, ...),
                   RP=rmstRP(formula, data, trunc, alpha, ...),
                   stop("Please enter a valid method."))

  result$call <- Call
  result
}

##' @method plot rmst
##' @export
plot.rmst <- function(x, ...) stop("No plots available for this method")

#########################################################################################################
#########################################################################################################

# Function required for rmstKM

rmstfunc <- function(dat, tau, alpha){

  indat <- dat

  # Taken from rmst2 package and embedded here for version control
  onearm <- function(time, status, tau, alpha, arm){
    ft         <- survival::survfit(Surv(time, status) ~ 1)
    idx        <- ft$time <= tau
    wk.time    <- sort(c(ft$time[idx], tau))
    wk.surv    <- ft$surv[idx]
    wk.n.risk  <- ft$n.risk[idx]
    wk.n.event <- ft$n.event[idx]
    time.diff  <- diff(c(0, wk.time))
    areas      <- time.diff * c(1, wk.surv)
    rmst       <- sum(areas)
    wk.var     <- ifelse((wk.n.risk - wk.n.event) == 0, 0, wk.n.event/(wk.n.risk * (wk.n.risk - wk.n.event)))
    wk.var     <- c(wk.var, 0)
    rmst.var   <- sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
    rmst.se    <- sqrt(rmst.var)
    out        <- matrix(0, 1, 5)
    out[1, ]   <- c(paste("RMST ",arm,":",sep=""),
                    rmst, rmst.se, rmst - qnorm(1 - alpha/2) * rmst.se,
                    rmst + qnorm(1 - alpha/2) * rmst.se)
    return(list("RMST"=out,"var"=rmst.var))
  }

  multiarm <- function(val){
    onearm(indat$time  [indat$arm == val],
           indat$status[indat$arm == val],
           tau   = tau,
           alpha = alpha,
           arm   = val)
  }

  uniqarm <- as.character(unique(sort(indat$arm)))
  allarms <- lapply(uniqarm,multiarm)

  if (length(uniqarm) > 1) {

    combi <- expand.grid(seq(1,length(uniqarm)),seq(1,length(uniqarm)))
    combi <- t(apply(combi, 1, sort))
    combi <- combi[combi[,1] != combi[,2],]
    combi <- matrix(combi[!duplicated(combi),],ncol=2)

    difference <- function(comb){
      dat1 <- allarms[[comb[1]]]
      dat2 <- allarms[[comb[2]]]
      rmst.diff.10     <- as.numeric(dat1$RMST[2]) - as.numeric(dat2$RMST[2])
      rmst.diff.10.se  <- sqrt(dat1$var + dat2$var)
      rmst.diff.10.low <- rmst.diff.10 - qnorm(1 - alpha/2) * rmst.diff.10.se
      rmst.diff.10.upp <- rmst.diff.10 + qnorm(1 - alpha/2) * rmst.diff.10.se
      rmst.diff.pval   <- pnorm(-abs(rmst.diff.10)/rmst.diff.10.se) * 2
      rmst.diff.result <- c(paste("RMST Dif (",uniqarm[comb[1]]," - ",uniqarm[comb[2]],"):",sep=""),
                            rmst.diff.10,rmst.diff.10.se, rmst.diff.10.low,
                            rmst.diff.10.upp, rmst.diff.pval)
      out              <- matrix(0, 1, 6)
      out[1, ]         <- rmst.diff.result
      return(out)
    }

    out22          <- t(apply(combi, 1, difference))
    out2           <- matrix(out22[,-1],ncol=5)
    class(out2)    <- "numeric"
    rownames(out2) <- out22[,1]
    colnames(out2) <- c("Est.",
                        "se",
                        paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                        paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""),
                        "p")

  } else {out2 <- NULL}


  out11 <- lapply(paste("allarms[[",seq(1,length(uniqarm)),"]]$RMST",sep=""),function(x) eval(parse(text=x)) )
  out11 <- as.matrix(do.call(rbind, out11))
  out1  <- matrix(out11[,-1],ncol=4)
  class(out1)    <- "numeric"
  rownames(out1) <- out11[,1]
  colnames(out1) <- c("Est.",
                      "se",
                      paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                      paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""))

  return(list("RMST"=out1,"diff"=out2))
}

#########################################################################################################
#########################################################################################################

# Function required for rmstPseudo

pseudofunc <- function(indat, tau, alpha){

  indat$pseudo <- pseudo::pseudomean(time  = indat$time,
                             event = indat$status,
                             tmax  = tau)
  indat$id     <- 1:dim(indat)[1]
  indat$arm    <- factor(indat$arm)

  if (length(unique(indat$arm))==1){
    rmstestimate <- mean(indat$pseudo)
    rmstse <- sd(indat$pseudo)/sqrt(length(indat$pseudo))
    lci <- rmstestimate - qnorm(1 - alpha/2) * rmstse
    uci <- rmstestimate + qnorm(1 - alpha/2) * rmstse
    arms <- data.frame("estimate"=rmstestimate,"san.se"=rmstse,"lci"=lci,"uci"=uci)
    arms <- arms[,c("estimate","san.se","lci","uci")]
    colnames(arms)<- c("Est.",
                       "se",
                       paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                       paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""))
    row.names(arms) <- paste("RMST ",row.names(arms),":",sep="")

    diff <- NULL
  } else {


  uniqarm <- as.character(unique(sort(indat$arm)))
  allarms <- lapply(uniqarm, functest, data=indat)

  combi <- expand.grid(seq(1,length(uniqarm)),seq(1,length(uniqarm)))
  combi <- t(apply(combi, 1, sort))
  combi <- combi[combi[,1] != combi[,2],]
  combi <- matrix(combi[!duplicated(combi),],ncol=2)

  difffunc <- function(input) {
    dif <- allarms[[input[1]]][input[2],]
    rownames(dif) <- c(paste("RMST Dif (",uniqarm[input[1]]," - ",uniqarm[input[2]],")",sep=""))
    return(dif)
  }

  arms <- do.call(rbind, lapply(allarms, head, 1))
  arms$lci <- arms$estimate - qnorm(1 - alpha/2) * arms$san.se
  arms$uci <- arms$estimate + qnorm(1 - alpha/2) * arms$san.se
  arms <- arms[,c("estimate","san.se","lci","uci")]
  colnames(arms)<- c("Est.",
                     "se",
                     paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                     paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""))
  row.names(arms) <- paste("RMST ",row.names(arms),":",sep="")

  diff <- do.call("rbind", apply(combi, 1, difffunc))
  diff[,1] <- -diff[,1]
  diff$lci <- diff$estimate - qnorm(1 - alpha/2) * diff$san.se
  diff$uci <- diff$estimate + qnorm(1 - alpha/2) * diff$san.se
  diff <- diff[,c("estimate","san.se","lci","uci","p")]
  colnames(diff)<- c("Est.",
                     "se",
                     paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                     paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""),
                     "p")
  }

  return(list("RMST"=arms,"diff"=diff))
}

functest <- function(ar,data){
    indat <- base::within(data, arm <- relevel(arm, ref = ar))

    fit1 <- geepack::geese(pseudo ~ arm,
                  data = indat,
                  id=indat$id,
                  jack = TRUE,
                  family=gaussian,
                  corstr="independence",
                  scale.fix=FALSE)

    sf11 <- summary(fit1)
    rownames(sf11$mean)[1] <- c(paste(ar,sep=""))
    return(sf11$mean)
}

#########################################################################################################
#########################################################################################################

# Function required for rmstRP

# Single arm
surv <- function(x, var, dframe, slength, RPfit){
  unlist(modflexsummary(var,
                        slength,
                        RPfit,
                        t=x,
                        ci=FALSE)[[1]]["est"])
}

RMST <- function(var, end, slength, RPfit){
  integrate(surv,0,end,var=var,slength=slength,RPfit=RPfit)$value
}

SE <- function(var, end, slength, RPfit, modelvcov){
  gd <- numDeriv::grad(RMST, var, end = end, slength=slength,RPfit=RPfit)
  sqrt(gd %*% modelvcov %*% gd)
}

# Multiple arms
surv_mult <- function(x,var,dframe,slength,RPfit){
  unlist(modflexsummary(var,
                        slength,
                        RPfit,
                        t=x,
                        ci=FALSE,
                        newdata=dframe)[[1]]["est"])
}

RMST_mult <- function(var, end, dframe, slength, RPfit){
  integrate(surv_mult,0,end,dframe=dframe,var=var,slength=slength,RPfit=RPfit)$value
}

SE_mult <- function(var, end, dframe, slength, RPfit, modelvcov){
  gd <- numDeriv::grad(RMST_mult, var, end = end, dframe = dframe, slength=slength, RPfit=RPfit)
  sqrt(gd %*% modelvcov %*% gd)
}

diffunc_mult <- function(var, end, dframe, slength, RPfit) {
  abs(integrate(surv_mult,0,end,dframe=dframe[[1]],var=var, slength=slength, RPfit=RPfit)$value -
      integrate(surv_mult,0,end,dframe=dframe[[2]],var=var, slength=slength, RPfit=RPfit)$value )}

SEdif_mult <- function(var, end, dframe, slength, RPfit, modelvcov){
  gd <- numDeriv::grad(diffunc_mult, var, end = end, dframe = dframe, slength=slength, RPfit=RPfit)
  sqrt(gd %*% modelvcov %*% gd)
}

#########################################################################################################
#########################################################################################################
# Imported function from flexsurv

form.model.matrix <- function (object, newdata) 
{
  mfo <- model.frame(object)
  covnames <- attr(mfo, "covnames")
  missing.covs <- unique(covnames[!covnames %in% names(newdata)])
  if (length(missing.covs) > 0) {
    missing.covs <- sprintf("\"%s\"", missing.covs)
    plural <- if (length(missing.covs) > 1) 
      "s"
    else ""
    stop(sprintf("Value%s of covariate%s ", plural, plural), 
         paste(missing.covs, collapse = ", "), " not supplied in \"newdata\"")
  }
  tt <- attr(mfo, "terms")
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, xlev = .getXlevels(tt, 
                                                       mfo))
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, mf)
  forms <- object$all.formulae
  mml <- vector(mode = "list", length = length(object$dlist$pars))
  names(mml) <- names(forms)
  forms[[1]] <- delete.response(terms(forms[[1]]))
  for (i in names(forms)) {
    mml[[i]] <- model.matrix(forms[[i]], mf)
  }
  X <- compress.model.matrices(mml)
  attr(X, "newdata") <- mf
  X
}

summary.fns <- function (x, type) 
{
  switch(type, survival = function(t, start, ...) {
    ret <- (1 - x$dfns$p(t, ...))/(1 - x$dfns$p(start, ...))
    ret[t < start] <- 1
    ret
  }, hazard = function(t, start, ...) {
    ret <- x$dfns$h(t, ...) * (1 - x$dfns$p(start, ...))
    ret[t < start] <- 0
    ret
  }, cumhaz = function(t, start, ...) {
    ret <- x$dfns$H(t, ...) - x$dfns$H(start, ...)
    ret[t < start] <- 0
    ret
  })
}

expand.summfn.args <- function (summfn) 
{
  summfn2 <- summfn
  args <- c(alist(t = , start = ), formals(summfn))
  formals(summfn2) <- args[!duplicated(names(args))]
  body(summfn2) <- body(summfn)
  summfn2
}

add.covs <- function (x, pars, beta, X, transform = FALSE) 
{
  nres <- nrow(X)
  if (!is.matrix(pars)) 
    pars <- matrix(pars, nrow = nres, ncol = length(pars), 
                   byrow = TRUE)
  if (!is.matrix(beta)) 
    beta <- matrix(beta, nrow = 1)
  for (j in seq(along = x$dlist$pars)) {
    covinds <- x$mx[[x$dlist$pars[j]]]
    if (length(covinds) > 0) {
      pars[, j] <- pars[, j] + beta[, covinds] %*% t(X[, 
                                                       covinds, drop = FALSE])
    }
    if (!transform) 
      pars[, j] <- x$dlist$inv.transforms[[j]](pars[, j])
  }
  colnames(pars) <- x$dlist$pars
  pars
}

cisumm.flexsurvreg <- function (x, t, start, X, fn, B = 1000, cl = 0.95) 
{
  if (any(is.na(x$res[, 2])) || (B == 0)) 
    ret <- array(NA, dim = c(length(t), 2))
  else {
    ret <- normbootfn.flexsurvreg(x = x, t = t, start = start, 
                                  X = X, fn = fn, B = B)
    ret <- apply(ret, c(1, 3), function(x) quantile(x, c((1 - 
                                                            cl)/2, 1 - (1 - cl)/2), na.rm = TRUE))
    ret <- t(ret[, 1, ])
  }
  ret
}

form.model.matrix <- function (object, newdata) 
{
  mfo <- model.frame(object)
  covnames <- attr(mfo, "covnames")
  missing.covs <- unique(covnames[!covnames %in% names(newdata)])
  if (length(missing.covs) > 0) {
    missing.covs <- sprintf("\"%s\"", missing.covs)
    plural <- if (length(missing.covs) > 1) 
      "s"
    else ""
    stop(sprintf("Value%s of covariate%s ", plural, plural), 
         paste(missing.covs, collapse = ", "), " not supplied in \"newdata\"")
  }
  tt <- attr(mfo, "terms")
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, xlev = .getXlevels(tt, 
                                                       mfo))
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, mf)
  forms <- object$all.formulae
  mml <- vector(mode = "list", length = length(object$dlist$pars))
  names(mml) <- names(forms)
  forms[[1]] <- delete.response(terms(forms[[1]]))
  for (i in names(forms)) {
    mml[[i]] <- model.matrix(forms[[i]], mf)
  }
  X <- compress.model.matrices(mml)
  attr(X, "newdata") <- mf  
  X
}

compress.model.matrices <- function (mml) 
{
  cbind.drop.intercept <- function(...) do.call("cbind", lapply(list(...), 
                                                                function(x) x[, -1, drop = FALSE]))
  X <- do.call("cbind.drop.intercept", mml)
  loc.cnames <- colnames(mml[[1]])[-1]
  anc.cnames <- unlist(mapply(function(x, y) sprintf("%s(%s)", 
                                                     x, y), names(mml[-1]), lapply(mml[-1], function(x) colnames(x)[-1])))
  cnames <- c(loc.cnames, anc.cnames)
  colnames(X) <- cnames
  X
}
