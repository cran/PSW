#' @title Mirror histogram
#' @description \code{psw.mirror.hist} is used to plot the mirror histogram that visualizes the propensity score distributions in both treatment groups.
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, \code{"OVERLAP"}, and \code{"TRAPEZOIDAL"}.
#' @param add.weight add propensity score weights to the mirror histogram, \code{add.weight=FALSE} by default and it is not available for \code{weight="ATE", "ATT"} or \code{"ATC"}.
#' @param nclass number of breaks in the mirror histogram.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight.
#' @details See \code{psw}.
#' @return \code{NULL}.
#' @seealso \link{psw}
#' @import graphics
#' @export
#' @examples
#' # Load the test data set
#' data(test_data);
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#' tmp <- psw.mirror.hist( data = test_data, weight = "MW", form.ps = form.ps,
#' add.weight = TRUE );
#'
psw.mirror.hist <- function( data, form.ps, weight, add.weight = FALSE, nclass=50, K=4 ) {

  # check if data is a data frame
  if ( !( is.data.frame( data ) ) ) {
    stop( "Input data must be a data frame." );
  }

  # disable add.weight option for "ATE", "ATT", and "ATC"
  if ( weight %in% c( "ATE", "ATT", "ATC" ) ) {
    add.weight <- FALSE;
  }

  # fit propensity score model
  out.ps <- ps.model( dat = data, form.ps = as.formula( form.ps ) );
  trt.var <- as.character( terms( out.ps$fm )[[ 2 ]] );
  ps.hat <- out.ps$ps.hat;  # estimated ps
  Xname <- names( coef( out.ps$fm ) )[-1];  # covariate name in propensity score model
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

  mirror.hist.core( ps.above = ps.hat[ data[ , trt.var] == 0 ],
                    ps.below = ps.hat[ data[ , trt.var] == 1 ],
                    wt.above = W[ data[ , trt.var] == 0 ],
                    wt.below = W[ data[ , trt.var] == 1 ],
                    add.weight = add.weight,
                    label.above = paste( trt.var, "=0", sep="" ),
                    label.below = paste( trt.var, "=1", sep="" ),
                    nclass = nclass );
}
