
#' @title Propensity score weighting with augmented estimation
#' @description \code{psw.aug} is the function to estimate the augmented estimator for mean difference
#' (mean outcome difference for \code{"gaussian"} family and risk difference for \code{"binomial"} family).
#' The augmented estimator is consistent for the estimand defined by the corresponding propensity score model.
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, \code{"OVERLAP"}, and \code{"TRAPEZOIDAL"}.
#' @param form.outcome outcome model.
#' @param family outcome family, either \code{"gaussian"} or \code{"binomial"}. \code{family="gaussian"} by default.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight. The estimand is
#' closer to the average treatment effect (ATE) with larger value of \code{K}. \code{K=4} by default.
#' @details \code{psw.aug} is used to estimate the augmented estimator, \eqn{\hat{\Delta}_{aug}},
#' and make inference using the sandwich variance that adjusts for the sampling variability in the estimated propensity score.
#' @return A list of weighting method, fitted propensity score model, estimated propenstity scores, estimated propensity score weights,
#' augmented estimator and associated standard error.
#' \item{weight}{weighting method.}
#' \item{ps.model}{object returned by fitting the propensity score model using \code{glm} with \code{"binomial"} family.}
#' \item{ps.hat}{estimated propensity score.}
#' \item{W}{estimated propensity score weight.}
#' \item{est.aug}{augmented estimator for mean difference when \code{family = "gaussian"}.}
#' \item{std.aug}{standard error for \code{est.aug}.}
#' \item{est.risk.aug}{augmented estimator for risk difference when \code{family = "binomial"}.}
#' \item{std.risk.aug}{standard error for \code{est.risk.aug}.}
#' @export
#' @examples
#' # Load the test data set
#' data(test_data);
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#' # Outcome model
#' form.out <- "Y ~ X1 + X2 + X3 + X4";
#' tmp <- psw.aug( data = test_data, form.ps = form.ps, weight = "ATE",
#' form.outcome = form.out, family="gaussian" );
#'
psw.aug <- function( data, form.ps, weight, form.outcome, family="gaussian", K=4 ) {

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
  res.aug <- psw.aug.core( dat = data, beta.hat = beta.hat, omega = omega, Q = Q, out.ps = out.ps, out.outcome = out.outcome,
                         trt.var = trt.var, out.var = out.var, weight = weight, family = family,  K = K, delta = 0.002 );
  if ( family == "gaussian" ) {
    res$est.aug <- res.aug$est;
    res$std.aug <- res.aug$std;
  }
  if ( family == "binomial" ) {
    # risk difference
    res$est.risk.aug <- res.aug$est.risk;
    res$std.risk.aug <- res.aug$std.risk;
  }

  return( res );
}
