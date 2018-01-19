
#' @title Balance checking using standardized mean difference
#' @description \code{psw.balance} is used to compute the standardized mean difference (in percentage) for balance diagnosis.
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, \code{"OVERLAP"}, and \code{"TRAPEZOIDAL"}.
#' @param V.name a vector of covariates on which standardized mean difference is computed. If \code{V.name = NULL}, the covariates in propensity score model are used.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight.
#' @return A list of weighting method, fitted propensity score model, estimated propenstity scores, estimated propensity score weights,
#' standardized mean difference before and after weighting adjustment.
#' \item{weight}{weighting method.}
#' \item{ps.model}{object returned by fitting the propensity score model using \code{glm} with \code{"binomial"} family.}
#' \item{ps.hat}{estimated propensity score.}
#' \item{W}{estimated propensity score weight.}
#' \item{std.diff.before}{A data frame of weighed mean, variance, and standardized mean difference for covariates in \code{V.name} by treatment groups before weighting.
#' \code{V.name} is the row name and \code{"treated.mean"}, \code{"treated.var"}, \code{"control.mean"}, \code{"control.var"}, \code{"std.diff.pct"} are column names.}
#' \item{std.diff.after}{A data frame of weighed mean, variance, and standardized mean difference for covariates in \code{V.name} by treatment groups after weighting.}
#' @seealso \link{psw}, \link{psw.spec.test}
#' @export
#' @examples
#' # Load the test data set
#' data(test_data);
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#' # A vector of covariates
#' V.name <- c( "X1", "X2", "X3", "X4" );
#' tmp <- psw.balance( data = test_data, weight = "MW", form.ps = form.ps,
#' V.name = V.name );
#'
psw.balance <- function(data, form.ps, weight, V.name = NULL, K = 4 ) {

  # check if data is a data frame
  if ( !( is.data.frame( data ) ) ) {
    stop( "Input data must be a data frame." );
  }

  n <- nrow(data);  # number of sujects
  out.ps <- ps.model( dat = data, form.ps = form.ps );
  trt.var <- as.character( terms( out.ps$fm )[[ 2 ]] );
  Xname <- names( coef( out.ps$fm ) )[-1];  # covariate name in propensity score model

  # covariates for balance checking
  if ( is.null( V.name ) ) {
    V.name <- Xname;
  }

  ps.hat <- out.ps$ps.hat;  # estimated ps
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

  V.mat <- as.matrix( data[ , V.name, drop=F ] ) ;
  trt <- as.numeric( data[ , trt.var ] ) ;
  tmp <- psw.balance.core( Xmat = V.mat, Xname = V.name,
                           weight = weight, W = W, Z = as.numeric( data[ , trt.var ] ) );

  res <- list( weight = weight,
               ps.model = out.ps$fm,
               ps.hat = ps.hat,
               W = W,
               std.diff.before = tmp$std.diff.before,
               std.diff.after = tmp$std.diff.after
             );

  return( res );
}
