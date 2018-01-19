
#' @title Propensity score model specification test
#' @description \code{psw.spec.test} is used to test the sufficiency of propensity score model in balancing covariates between groups.
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, and \code{"TRAPEZOIDAL"}.
#' Note that specification test is not available for \code{weight = "OVERLAP"}.
#' @param V.name a vector of covariates on which the specification test is performed.
#' @param trans.type a vector of the same length as \code{V.name} that specifies the transformation type for each element in \code{V.name}.
#' Available transformations are \code{"identity"}, \code{"log"}, \code{"logit"}, \code{"sqrt"}, and \code{"Fisher"}.
#' No transformation is perfomred with \code{"identity"}.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight. \code{K=4} by default.
#' @return A list of model specification test results.
#' \item{weight}{weighting method.}
#' \item{ps.model}{object returned by fitting the propensity score model using \code{glm} with \code{"binomial"} family.}
#' \item{ps.hat}{estimated propensity score.}
#' \item{W}{estimated propensity score weight.}
#' \item{V.name}{covariates in the specification test.}
#' \item{g.B1.hat}{a vector of transformed weighted average for covariates in the treated group when \code{spec.test=T}.}
#' \item{g.B0.hat}{a vector of transformed weighted average for covariates in the control group when \code{spec.test=T}.}
#' \item{B.hat}{difference between \code{eta.B1.hat} and \code{eta.B0.hat} when \code{spec.test=T}.}
#' \item{var.B.hat}{covariance matrix for \code{B.hat} when \code{spec.test=T}.}
#' \item{test.stat}{test statistic for the specification test, which follows the \eqn{\chi^2_{df}} distribution under the null, available when \code{spec.test=T}.}
#' \item{df}{degree of freedom for the specification test, \code{df=rank(var.B.hat)}, available when \code{spec.test=T}.}
#' \item{pvalue}{pvalue of the specification test when \code{spec.test=T}.}
#' @details In the data set, treatment indicator should be numerically specified such that a value of \code{1} indicates the treated
#' and a value of \code{0} indicates the control. The null hypothesis is that the propensity score model is correctly specified; the
#' alternative is that the propensity score model is misspecified. Therefore, this test is a goodness-of-fit test of propensity score model,
#' with the test statistic being a metric of covariate balance.
#' #'
#'   Rejection of the specification test implies current propensity score model is inadquate
#' for balancing covariates between groups.
#' @seealso \link{psw}, \link{psw.balance}
#' @export
#' @examples
#' # Load the test data set
#' data(test_data);
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#' # A vector of covariates
#' V.name <- c( "X1", "X2", "X3", "X4" );
#' # A vector of transformation types for covariates in V.name.
#' trans.type <- c( "identity", "identity", "logit", "logit" );
#' tmp <- psw.spec.test( data = test_data, form.ps = form.ps,
#' weight = "MW", V.name = V.name, trans.type = trans.type );
#'
psw.spec.test <- function( data, form.ps, weight, V.name, trans.type, K=4 ) {

  # check if data is a data frame
  if ( !( is.data.frame( data ) ) ) {
    stop( "Input data must be a data frame." );
  }
  if ( weight == "OVERLAP" ) {
    stop( "Specification test is not available for OVERLAP weight." );
  }

  n <- nrow(data);  # number of sujects
  out.ps <- ps.model( dat = data, form.ps = form.ps );
  X.name <- names( coef( out.ps$fm ) )[-1]
  trt.var <- as.character( terms(out.ps$fm)[[2]] );  # treatment indicator

  ps.hat <- out.ps$ps.hat;
  Q <- data[ , trt.var ] * ps.hat + ( 1 - data[ , trt.var ] ) * ( 1 - ps.hat );
  omega <- sapply( ps.hat, calc.omega, weight = weight, delta = 0.002, K = K );
  W <- omega/Q;  # generic weight
  beta.hat <- as.numeric( coef(out.ps$fm) );  # ps model coefficients
  if ( weight == "MW" ) {
    # adjust the MW weight ( for MW, the omega function is multiplied by 2 )
    W <- W / 2;
  }

  X.mat <- as.matrix( cbind( rep( 1, n ), data[ , X.name, drop=F ] ) ) ;
  V.mat <- as.matrix( data[ , V.name, drop=F ] ) ;
  trt <- as.numeric( data[ , trt.var ] ) ;  # treatment value

  ans <- psw.spec.test.core( X.mat = X.mat, V.mat = V.mat, V.name =V.name, trt = trt, beta.hat = beta.hat, omega = omega,
                             Q = Q, trans.type = trans.type, weight = weight, delta = 0.002, K = K );

  res <- list( weight = weight,
               ps.model = out.ps$fm,
               ps.hat = ps.hat,
               W = W,
               V.name = ans$V.name,
               g.B1.hat = ans$g.B1.hat,
               g.B0.hat = ans$g.B0.hat,
               B.hat = ans$B.hat,
               var.B.hat = ans$var.B.hat,
               test.stat = ans$test.stat,
               df = ans$df,
               pvalue = ans$pvalue );

  return( res );
}
