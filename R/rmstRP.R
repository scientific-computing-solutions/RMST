##' Calculate the RMST in each arm and the between-group difference (Royston-Parmar model)
##'
##' This function calculates an estimate of the RMST by fitting a Royston-Parmar model to the survival data with arm as a covariate. The model is fit using the function \code{flexsurvspline} from the package \code{flexsurv}. For further information see the \code{flexsurv} documentation. The standard errors are calculated using the Delta method.
##' @return The Restricted Mean Survival Time computated by fitting Royston-Parmar model to the survival data
##' @param formula a formula object, with the response on the left of a ~ operator, and the a single term for arm on the right. The response must be a survival object as returned by the Surv function.
##' @param data a data.frame in which to interpret the variables named in the formula.
##' @param trunc the time point to truncated the estimate of the RMST at.
##' @param alpha The default is 0.05. (1-alpha) confidence intervals are reported.
##' @param knots the number of internal knots to fit the Royston-Parmar model with.
##' @examples
##' # Load in an example data set
##' D <- rmst2.sample.data()
##'
##' # Calculate the RMST based on the Royston-Parmar model with two knots
##' rmstRP(Surv(time, status) ~ arm, data=D, trunc=5, alpha=0.05, knots=2)
##' @export
rmstRP <- function (formula, data, trunc, alpha = 0.05, knots = 2) {
  Call <- match.call()
  temp <- Call[c(1, match(c("formula", "data"), names(Call), nomatch=0))]
  temp[[1]] <- as.name("model.frame")

  m <- eval.parent(temp)

  response <- model.extract(m, "response")

  ## check that response is an appropriately-sized thing
  if (!survival::is.Surv(response)) {
    stop("Response must be a survival object")
  }

  if (attr(response, "type") != "right") {
    stop("Require right-censored data")
  }

  Terms <- terms(formula, c("strata", "cluster"))
  ord <- attr(Terms, "order")
  if (length(ord) & any(ord != 1))
    stop("Interaction terms are not valid for this function")

  ll <- attr(Terms, "term.labels")
  if (length(ll) > 1) {
    stop("The survival formula must only contain arm as a covariate")
  }

  if (!is.numeric(trunc) | trunc < 0){
    stop("The trunction time must be a positive number")
  }

  if (identical(ll, character(0))){
    arm <- gl(1,1)
  } else {
  arm <- survival::strata(m[ll])
  armlab <- levels(arm)
  }

  if (nlevels(arm) > 2) {
    stop("The arm variable can't contain more than two levels for the RP method currently")
  }

  my.levels <- unique(m[ll])

  if (nlevels(arm)==1){
    indat <- data.frame("time" = as.numeric(m[,1])[1:dim(m)[1]],
                        "status" = as.numeric(m[,1])[-1:-dim(m)[1]],
                        "arm"  = 1)

    RPfit <- flexsurv::flexsurvspline(Surv(time, status) ~ 1, data=indat, k=2)

    var       <- c(RPfit$res.t[RPfit$dlist$pars,"est"], RPfit$res[RPfit$covpars,"est"])
    modelvcov <- vcov(RPfit)
    slength   <- length(RPfit$res.t[RPfit$dlist$pars,"est"])
    
    RM.est <- RMST(var,trunc,slength,RPfit)
    se.est <- SE(var,trunc,slength,RPfit,modelvcov)

    RMout <- data.frame(" "        = "",
                        "RMST"     = RM.est,
                        "SE"       = se.est,
                        "Lower CI" = RM.est - qnorm(1-alpha/2)*se.est,
                        "Upper CI" = RM.est + qnorm(1-alpha/2)*se.est)

    names(RMout) <- c("","RMST","SE","Lower CI","Upper CI")

    RMSTdif <- NULL
  } else {

  RPfit <- flexsurv::flexsurvspline(formula, data=data, k=knots)

  var       <- c(RPfit$res.t[RPfit$dlist$pars,"est"], RPfit$res[RPfit$covpars,"est"])
  modelvcov <- vcov(RPfit)
  slength   <- length(RPfit$res.t[RPfit$dlist$pars,"est"])

  ndata1 <- data.frame("test" = my.levels[,1][1])
  ndata2 <- data.frame("test" = my.levels[,1][2])
  colnames(ndata1) <- c(names(m[ll])[1])
  colnames(ndata2) <- c(names(m[ll])[1])
  df <- list(ndata1,ndata2)

  RM.est <- lapply(df, RMST_mult, var=var, end=trunc, slength=slength, RPfit=RPfit)
  se.est <- lapply(df, SE_mult,   var=var, end=trunc, slength=slength, RPfit=RPfit, modelvcov=modelvcov)

  RMout <- data.frame(" "        = c(armlab[1],armlab[2]),
                      "RMST"     = unlist(RM.est),
                      "SE"       = unlist(se.est),
                      "Lower CI" = unlist(RM.est) - qnorm(1-alpha/2)*unlist(se.est),
                      "Upper CI" = unlist(RM.est) + qnorm(1-alpha/2)*unlist(se.est))

  names(RMout) <- c("","RMST","SE","Lower CI","Upper CI")

  dif    <- diffunc_mult(var, trunc, df, slength=slength, RPfit=RPfit)
  se.dif <- SEdif_mult(var, trunc, df, slength=slength, RPfit=RPfit, modelvcov=modelvcov)

  RMSTdif <- data.frame("Difference" = dif,
                        "SE" = se.dif,
                        "Lower CI" = dif - qnorm(1-alpha/2)*se.dif,
                        "Upper CI" = dif + qnorm(1-alpha/2)*se.dif,
                        "p-value"   = 1 - pchisq((dif/se.dif)^2, df = 1))
}

  result <- list("RMST"=RMout,"diff"=RMSTdif, call=Call)
  class(result) <- c("rmstRP", "rmst")
  return(result)
}
