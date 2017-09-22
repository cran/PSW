
#' @title Propensity score weighting with doubly robust estimation
#' @description \code{psw.dr} is the function to estimate the doubly robust estimator for mean difference
#' (mean outcome difference for \code{"gaussian"} family and risk difference for \code{"binomial"} family).
#' The doubly robust estimator is consistent when either the propensity score model or the outcome model is correctly specified.
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, \code{"OVERLAP"}, and \code{"TRAPEZOIDAL"}.
#' @param form.outcome outcome model.
#' @param family outcome family, either \code{"gaussian"} or \code{"binomial"}. \code{family="gaussian"} by default.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight. The estimand is
#' closer to the average treatment effect (ATE) with larger value of \code{K}. \code{K=4} by default.
#' @details \code{psw.dr} is used to estimate the doubly robust estimator, \eqn{\hat{\Delta}_{DR}},
#' and make inference using the sandwich variance that adjusts for the sampling variability in the estimated propensity score.
#' @return A list of weighting method, fitted propensity score model, estimated propenstity scores, estimated propensity score weights,
#' doubly robust estimator and standard deviation estimator.
#' \item{weight}{weighting method.}
#' \item{ps.model}{object returned by fitting the propensity score model using \code{glm} with \code{"binomial"} family.}
#' \item{ps.hat}{estimated propensity score.}
#' \item{W}{estimated propensity score weight.}
#' \item{est.dr}{doubly robust estimator for mean difference when \code{family = "gaussian"}.}
#' \item{std.dr}{standard deviation for \code{est.dr}.}
#' \item{est.risk.dr}{doubly robust estimator for risk difference when \code{family = "binomial"}.}
#' \item{std.risk.dr}{standard deviation for \code{est.risk.dr}.}
#' @export
#' @examples
#' # Load the test data set
#' data(test_data);
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#' # Outcome model
#' form.out <- "Y ~ X1 + X2 + X3 + X4";
#' tmp <- psw.dr( data = test_data, weight = "ATE", form.ps = form.ps,
#' form.outcome = form.out, family="gaussian" );
#'
psw.dr <- function( data, form.ps, weight, form.outcome, family="gaussian", K=4 ) {

  # Outcome family can only be "gaussian" or "binomial"
  if ( !( family %in% c( "gaussian", "binomial" ) ) ) {
    stop( "Family should be gaussian or binomial." );
  }

  # check if data is a data frame
  if ( !( is.data.frame( data ) ) ) {
    stop( "Input data must be a data frame." );
  }

  # fit propensity score model
  out.ps <- ps.model( dat = data, form.ps = as.formula( form.ps ) );
  trt.var <- as.character( terms( out.ps$fm )[[ 2 ]] );
  Xname <- names( coef( out.ps$fm ) )[-1];  # covariate name in propensity score model

  ps.hat <- out.ps$ps.hat;  # estimated ps
  beta.hat <- as.numeric(coef(out.ps$fm));
  omega <- sapply( ps.hat, calc.omega, weight = weight, delta = 0.002, K = K);
  Q <- data[ , trt.var ] * ps.hat + ( 1 - data[ , trt.var ] ) * ( 1 - ps.hat );  # denominator of generic weight;
  W <- omega/Q;  # generic weight;
  if ( weight == "MW" ) {
    # adjust the MW weight ( for MW, the omega function is multiplied by 2 )
    W <- W / 2;
  }
  if ( weight == "OVERLAP" ) {
    # adjust the OVERLAP weight ( for MW, the omega function is multiplied by 4 )
    W <- W / 4;
  }

  res <- list( weight = weight,
               ps.model = out.ps$fm,
               ps.hat = ps.hat,
               W = W
             );

  out.outcome <- outcome.model( dat=data, form=as.formula(form.outcome), trt.var=trt.var, family=family );  # fit outcome model;
  out.var <- as.character( terms( out.outcome$fm1 )[[2]] );
  res.dr <- psw.dr.core( dat = data, beta.hat = beta.hat, omega = omega, Q = Q, out.ps = out.ps, out.outcome = out.outcome,
                         trt.var = trt.var, out.var = out.var, weight = weight, family = family, delta = 0.002, K = K );
  if ( family == "gaussian" ) {
    res$est.dr <- res.dr$est;
    res$std.dr <- res.dr$std;
  }
  if ( family == "binomial" ) {
    # risk difference
    res$est.risk.dr <- res.dr$est.risk;
    res$std.risk.dr <- res.dr$std.risk;
  }

  return( res );
}
