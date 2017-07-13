##' Calculate the RMST in each arm and the between-group difference (Pseudovalues)
##'
##' This function calculates an estimate of the RMST by generating a pseudovalue for each patient, and then fitting a GEE model to these pseudovalues with arm as a term in the model. The pseudovalues are generated using the function \code{pseudomean} from the package \code{pseudo} and the GEE is fit using the function \code{geese} from the package \code{geepack}. For further information see the \code{pseudo} or \code{geepack} documentation.
##' @return The Restricted Mean Survival Time computated by fitting a GEE to pseudovalues generated for each patient
##' @param formula a formula object, with the response on the left of a ~ operator, and the a single term for arm on the right. The response must be a survival object as returned by the Surv function.
##' @param data a data.frame in which to interpret the variables named in the formula.
##' @param trunc the time point to truncated the estimate of the RMST at.
##' @param alpha The default is 0.05. (1-alpha) confidence intervals are reported.
##' @examples
##' # Load in an example data set
##' D <- rmst2.sample.data()
##'
##' # Calculate the RMST based on the pseudovalue approach
##' rmstPseudo(Surv(time, status) ~ arm, data=D, trunc=5, alpha=0.05)
##' @export
rmstPseudo <- function (formula, data, trunc, alpha = 0.05) {
  Call <- match.call()
  temp <- Call[c(1, match(c("formula", "data"), names(Call), nomatch=0))]
  temp[[1]] <- as.name("model.frame")

  m        <- eval.parent(temp)
  response <- model.extract(m, "response")
  Terms    <- terms(formula, c("strata", "cluster"))
  ord      <- attr(Terms, "order")
  ll       <- attr(Terms, "term.labels")

  if (!survival::is.Surv(response)) {
    stop("Response must be a survival object")
  }

  if (attr(response, "type") != "right") {
    stop("Require right-censored data")
  }

  if (length(ord) & any(ord != 1))
    stop("Interaction terms are not valid for this function")

  if (length(ll) > 1) {
    stop("The survival formula must only contain arm as a covariate")
  }

  if (!is.numeric(trunc) | trunc < 0){
    stop("The trunction time must be a positive number")
  }

  if (identical(ll,character(0))){
    m[,2] <- 1
  }

  indat <- data.frame("time" = as.numeric(m[,1])[1:dim(m)[1]],
                      "status" = as.numeric(m[,1])[-1:-dim(m)[1]],
                      "arm"  = m[,2])

  if (any(tapply(indat$time, indat$arm, max) < trunc)) {
    stop(paste("The truncation time must be shorter than the minimum of the largest observed time in each group: ",
               sprintf("%.3f",min(tapply(indat$time, indat$arm, max))),
               sep=""))
  }

  out <- pseudofunc(indat, tau=trunc, alpha=alpha)

  result <- list("RMST"=out$RMST, "diff"=out$diff, call=Call)
  class(result) <- c("rmstPseudo", "rmst")
  return(result)
}
