
#' @title Propensity score weighting estimator
#' @description \code{psw.wt} is used to estimate the weighted treatment effect estimator (without double robustness).
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, \code{"OVERLAP"}, and \code{"TRAPEZOIDAL"}.
#' @param out.var outcome variable.
#' @param family outcome family, either \code{"gaussian"} or \code{"binomial"}. \code{family="gaussian"} by default.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight. The estimand is
#' closer to the average treatment effect (ATE) with larger value of \code{K}. \code{K=4} by default.
#' @details \code{psw.wt} is used to estimate the weighted estimator, \eqn{\hat{\Delta}}, and make inference using the sandwich variance estimator
#' that takes into account the sampling variability in the estimated propensity score.
#' @return A list of weighting method, fitted propensity score model, estimated propenstity scores, estimated propensity score weights,
#' weighted estimator and standard deviation estimator
#' \item{weight}{weighting method.}
#' \item{ps.model}{object returned by fitting the propensity score model using \code{glm} with \code{"binomial"} family.}
#' \item{ps.hat}{estimated propensity score.}
#' \item{W}{estimated propensity score weight.}
#' \item{est.wt}{weighted estimator for mean difference when \code{family = "gaussian"}.}
#' \item{std.wt}{standard deviation for \code{est.wt}.}
#' \item{est.risk.wt}{weighted estimator for risk difference when \code{family = "binomial"}.}
#' \item{std.risk.wt}{standard deviation for \code{est.risk.wt}.}
#' \item{est.rr.wt}{weighted estimator for relative risk when \code{family = "binomial"}.}
#' \item{std.rr.wt}{standard deviation for \code{est.rr.wt}.}
#' \item{est.or.wt}{weighted estimator for odds ratio when \code{family = "binomial"}.}
#' \item{std.or.wt}{standard deviation for \code{est.or.wt}.}
#' \item{est.lor.wt}{weighted estimator for log odds ratio when \code{family = "binomial"}.}
#' \item{std.lor.wt}{standard deviation for \code{est.lor.wt}.}
#' @seealso \link{psw}
#' @import stats
#' @export
#' @examples
#' # Load the test data set
#' data(test_data);
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#' tmp <- psw.wt( data = test_data, weight = "ATE", form.ps = form.ps,
#' out.var = "Y", family = "gaussian" );
#'
psw.wt <- function( data, form.ps, weight, out.var, family="gaussian", K=4 ) {

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

  res.wt <- psw.wt.core( dat = data, beta.hat = beta.hat, omega=omega, Q = Q, trt.var = trt.var,
                         out.var = out.var, family = family, Xname = Xname, weight = weight, delta = 0.002, K = K );
  if ( family == "gaussian" ) {
    res$est.wt <- res.wt$est;
    res$std.wt <- res.wt$std;
  }
  if ( family == "binomial" ) {
    # risk difference
    res$est.risk.wt <- res.wt$est.risk;
    res$std.risk.wt <- res.wt$std.risk;

    # relative risk
    res$est.rr.wt <- res.wt$est.rr;
    res$std.rr.wt <- res.wt$std.rr;

    # odds ratio
    res$est.or.wt <- res.wt$est.or;
    res$std.or.wt <- res.wt$std.or;

    # log odds ratio
    res$est.lor.wt <- res.wt$est.lor;
    res$std.lor.wt <- res.wt$std.lor;
  }
  return( res );
}
