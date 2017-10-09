
#' @title Propensity score weighting
#' @description \code{psw} is the main function to perfrom propensity score weighting analysis for (1) visualization of the propensity score distribution in both treatment groups,
#' (2) covariate balance diagnosis, (3) propensity score model specification test, (4) treatment effect estimation and inference, and (5) doubly robust estimation when applicable.
#' @param data data frame to be used.
#' @param form.ps propensity score model.
#' @param weight weighting method to be used. Available methods are \code{"ATE"},  \code{"ATT"}, \code{"ATC"}, \code{"MW"}, \code{"OVERLAP"}, and \code{"TRAPEZOIDAL"}.
#' @param std.diff calculate standardized mean difference as a percentage, \code{std.diff=FALSE} by default.
#' @param mirror.hist mirror histogram showing the propensity score distributions in both treatment groups, \code{mirror.hist=FALSE} by default.
#' @param add.weight add propensity score weights to the mirror histogram, \code{add.weight=FALSE} by default and it is not available for \code{weight="ATE", "ATT"} or \code{"ATC"}.
#' @param nclass number of breaks in the mirror histogram.
#' @param wt estimate the weighted estimator, \code{wt=FALSE} by default.
#' @param out.var outcome variable, needed when \code{wt=TRUE}.
#' @param family outcome family, either \code{"gaussian"} or \code{"binomial"}, \code{family="gaussian"} by default.
#' @param dr estimate the doubly robust estimator, \code{dr=FALSE} by default.
#' @param form.outcome outcome model, needed when \code{dr=TRUE}.
#' @param spec.test propensity score model specification test, \code{spec.test=FALSE} by default.
#' Note that specification test is not available for \code{weight="OVERLAP"}.
#' @param V.name a vector of covariates on which the specification test is performed, needed when \code{spec.test=TRUE}.
#' @param trans.type a vector of the same length as \code{V.name} that specifies the transformation type for each element in \code{V.name}.
#' Available transformations are \code{"identity"}, \code{"log"}, \code{"logit"}, \code{"sqrt"}, \code{"Fisher"}.
#' Needed when \code{spec.test=T}, and no transformation is perfomred with \code{"identity"}. See Details.
#' @param K value of \eqn{K} in \eqn{\omega(e_i) = min(1, K min(e_i, 1-e_i)) } for \code{"TRAPEZOIDAL"} weight. The estimand is
#' closer to the average treatment effect (ATE) with larger value of \code{K}. \code{K=4} by default.
#'
#' @details In package \code{PSW}, treatment indicator (left handside of \code{form.ps}) should be dummy coded
#' such that a value of 1 indicates the treated and a value of 0 indicates the control. All categorical covariates need to be dummy coded too.
#' If the outcome belongs to the \code{"gaussian"} family, causal estimation is based on the mean differnce between treatment groups. If the
#' outcome belongs to the \code{"binomial"} family, causal estimation is based on risk difference, risk ratio, odds ratio or log odds ratio.
#' The Delta method is used for variance estimation when applicable.
#'
#' Let \eqn{Z_i} be the treatment indicator of subject \eqn{i}, \eqn{e_i} be the corresponding propensity score. Then
#' propensity score weight, \eqn{W_i}, is defined as \deqn{W_i = \frac{\omega(e_i)}{Z_i e_i + (1-Z_i)(1-e_i)},} where \eqn{\omega(e_i)} is a function depending
#' on \eqn{e_i}. For \code{"ATE"}, \eqn{\omega(e_i) = 1}, which leads to estimating the average treatment effect. For \code{"ATT"}, \eqn{\omega(e_i) = e_i},
#' which leads to estimating average treatment effect for the treated. For \code{"ATC"}, \eqn{\omega(e_i) = 1 - e_i}, which leads to estimating average treatment effect
#' for the controls. For \code{"MW"}, \eqn{\omega(e_i) = min( e_i, 1 - e_i )}. For \code{"OVERLAP"}, \eqn{\omega(e_i) = e_i(1 - e_i)}. For \code{"TRAPEZOIDAL"},
#' \eqn{\omega(e_i) = min( 1, K min( e_i, 1 - e_i ) )}. This type of weights are studied by Hirano, Imbens and Ridder (2003) and Li et al (2016).
#' The \eqn{\omega(e_i)} function is specified by the \code{weight} argument.
#'
#' The matching weight (\code{"MW"}) was proposed by Li and Greene (2013). The overlap weight (\code{"OVERLAP"}) was propsed by Li et al (2016).
#' These methods down weight subjects with propensity score close to 0 or 1. and hence improve the stability of computation.
#'
#' A mirror histogram is produced to visualize the propensity score distributions in both treatment groups. In the mirror histogram, above the horizontal line
#' is the histogram of the propensiy scores of the control group, below is that of the treated group. The vertical axis of the histogram is the frequency. When
#' \code{mirror.hist.weight=TRUE}, the height of the shaded bar is the summation of the weights of subjects within the corresponding propensity score stratum.
#' For weighting methods of \code{"ATE"}, \code{"ATT"}, \code{"ATC"}, \code{add.weight} is not recommended for visual display because weights may be larger than 1.
#'
#' Standardized mean difference for a covariate is defiend as
#'  \deqn{ \frac{100 (\bar{x}_1 - \bar{x}_0)}{\sqrt{\frac{s_1^2 + s_0^2}{2} } },}
#' where \eqn{\bar{x}_1} and \eqn{s_1^2} are weighted mean and standard deviation for the treated group, and \eqn{ \bar{x_0}} and \eqn{s_0^2}
#' are defined similarly for the control group. A plot showing the standardized mean difference before and after weighting adjustement will be generated to
#' facilitate covariate balance diagnosis. It is sometimes recommended that the absolute values of standardized mean differences of all covariates should be less
#' than \code{10\%} in order to claim covariate balance.
#'
#' For the proensity score model specification test (Li and Greene, 2013), the quantity of interest is
#' \deqn{ \hat{B} = \boldsymbol{g} \left\{ \frac{ \sum_{i=1}^n W_i Z_i \boldsymbol{V}_i}{\sum_{i=1}^n W_i Z_i}\right\} - \boldsymbol{g} \left\{ \frac{ \sum_{i=1}^n W_i (1-Z_i) \boldsymbol{V}_i}{\sum_{i=1}^n W_i (1-Z_i)}\right\}, }
#' where \eqn{\boldsymbol{V}_i} is a vector of covariates whose balance are examined, and \eqn{\boldsymbol{g}(.)} is a vector of monotone smooth transformations for the input.
#' Transformation type is specified by argument \code{trans.type}, and available transformation types are \code{"identity"}, \code{"log"}, \code{"logit"}, \code{"sqrt"}, \code{"Fisher"}.
#' These transformations are recommended to improve the finite sample performance of the specification test. Log transformation (\code{"log"}) and square root transformation (\code{"sqrt"})
#' are recommended for skewed data, Logit transformation (\code{"logit"}) for binary data, and Fisher z-transformation (\code{"Fisher"}) for bounded data between \eqn{(-1, 1)}.
#' The current version of model specification test is not available for \code{weight="OVERLAP"} because it results in zero standardized difference.
#'
#' For estimation of mean difference (\code{"gaussian"} family) or risk difference (\code{"binomial"} family), the weighted estimator is
#' \deqn{ \hat{\Delta} = \frac{\sum_{i=1}^n W_i Z_i Y_i}{\sum_{i=1}^n W_i Z_i} - \frac{\sum_{i=1}^n W_i (1-Z_i) Y_i}{\sum_{i=1}^n W_i (1-Z_i)}, }
#' and the doubly robust estimator is
#'\deqn{\hat{\Delta}_{DR} = \frac{ \sum_{i=1}^n \omega(e_i) \{ m_{1i} - m_{0i} \}}{ \sum_{i=1}^n \omega(e_i) } + \frac{ \sum_{i=1}^n W_i Z_i \{ Y_i - m_{1i} \}}{ \sum_{i=1}^n W_i Z_i } - \frac{ \sum_{i=1}^n W_i (1-Z_i) \{ Y_i - m_{0i} \}}{ \sum_{i=1}^n W_i (1-Z_i)},}
#'where \eqn{m_{1i} = E[Y_i | \boldsymbol{X_i}, Z_i=1]} is the conditional expectation of outcome when treated given covariates \eqn{\boldsymbol{X}_i},
#'and \eqn{m_{0i} = E[Y_i | \boldsymbol{X_i}, Z_i=0]} is the conditional expectation of outcome when control given covariates \eqn{\boldsymbol{X}_i}.
#'When the outcome belongs to the \code{"binomial"} family, the marginal probability is used to estimate risk ratio, odds ratio and log odds ratio.
#'Sandwich variance estimation is used to adjust for the uncertainty in the estimated propensity score (Li and Greene, 2013).
#'
#' @return \code{psw} returns a list of elements depending on the supplied arguments.
#' \item{weight}{weighting method.}
#' \item{ps.model}{object returned by fitting the propensity score model using \code{glm} with \code{"binomial"} family.}
#' \item{ps.hat}{estimated propensity score.}
#' \item{W}{estimated propensity score weight.}
#' \item{std.diff.before}{A data frame of weighed mean, variance, and standardized mean difference for covariates in \code{V.name} (see below) by treatment groups before weighting.
#' \code{V.name} is the row name and \code{"treated.mean"}, \code{"treated.var"}, \code{"control.mean"}, \code{"control.var"}, \code{"std.diff.pct"} are column names.}
#' \item{std.diff.after}{A data frame of weighed mean, variance, and standardized mean difference for covariates in \code{V.name} by treatment groups after weighting.}
#' \item{est.wt}{weighted estimator for mean difference when \code{wt=T} and \code{family = "gaussian"}.}
#' \item{std.wt}{standard error for \code{est.wt}.}
#' \item{est.dr}{doubly robust estimator for mean difference when \code{dr=T} and \code{family = "gaussian"}.}
#' \item{std.dr}{standard error for \code{est.dr}.}
#' \item{est.risk.wt}{weighted estimator for risk difference when \code{wt=T} and \code{family = "binomial"}.}
#' \item{std.risk.wt}{standard error for \code{est.risk.wt}.}
#' \item{est.risk.dr}{doubly robust estimator for risk difference when \code{dr=T} and \code{family = "binomial"}.}
#' \item{std.risk.dr}{standard error for \code{est.risk.dr}.}
#' \item{est.rr.wt}{weighted estimator for relative risk when \code{wt=T} and \code{family = "binomial"}.}
#' \item{std.rr.wt}{standard error for \code{est.rr.wt}.}
#' \item{est.or.wt}{weighted estimator for odds ratio when \code{wt=T} and \code{family = "binomial"}.}
#' \item{std.or.wt}{standard error for \code{est.or.wt}.}
#' \item{est.lor.wt}{weighted estimator for log odds ratio when \code{wt=T} and \code{family = "binomial"}.}
#' \item{std.lor.wt}{standard error for \code{est.lor.wt}.}
#' \item{V.name}{covariates in the specification test or balance diagnosis.}
#' \item{g.B1.hat}{a vector of transformed weighted average for covariates in the treated group when \code{spec.test=T}.}
#' \item{g.B0.hat}{a vector of transformed weighted average for covariates in the control group when \code{spec.test=T}.}
#' \item{B.hat}{difference between \code{eta.B1.hat} and \code{eta.B0.hat} when \code{spec.test=T}.}
#' \item{var.B.hat}{covariance matrix for \code{B.hat} when \code{spec.test=T}.}
#' \item{test.stat}{test statistic for the specification test, which follows the \eqn{\chi^2_{df}} distribution under the null, available when \code{spec.test=T}.}
#' \item{df}{degree of freedom for the specification test, \code{df=rank(var.B.hat)}, available when \code{spec.test=T}.}
#' \item{pvalue}{pvalue of the specification test when \code{spec.test=T}.}
#'
#' @examples
#' # Load the test data set
#' data(test_data);
#'
#' # Propensity score model
#' form.ps <- "Z ~ X1 + X2 + X3 + X4";
#'
#' #1. Standardized differnce with "ATE"
#' tmp1 <- psw( data = test_data, form.ps = form.ps, weight = "ATE" );
#'
#' #2. Mirror histogram and add weights to it with "MW".
#' tmp2 <- psw( data = test_data, form.ps = form.ps, weight = "MW",
#' add.weight = T );
#'
#' #3. Estimate average treatment effect with "ATE"
#' tmp3 <- psw( data = test_data, form.ps = form.ps, weight = "ATE", wt = TRUE,
#' out.var = "Y", family = "gaussian" );
#'
#' #4. Doubly robust estimator with "OVERLAP"
#' # outcome model
#' form.out <- "Y ~ X1 + X2 + X3 + X4";
#' tmp4 <- psw( data = test_data, form.ps = form.ps, weight = "OVERLAP", dr = TRUE,
#' form.outcome = form.out, family = "gaussian" );
#'
#' #5. Propensity score model specification test with "MW".
#' # A vector of covariates
#' V.name <- c( "X1", "X2", "X3", "X4" );
#' # A vector of transformation types for covariates in V.name.
#' trans.type <- c( "identity", "identity", "logit", "logit" );
#' tmp5 <- psw( data = test_data, form.ps = form.ps, weight = "MW", spec.test = TRUE,
#' V.name = V.name, trans.type = trans.type );
#'
#' @references Hirano K, Imbens GW and Ridder G. "Efficient estimation of average treatment effects using the estimated propensity score." Econometrica 2003; 71(4): 1161-1189.
#' @references Li F, Morgan KL and Zaslavsky AM. "Balancing covariates via propensity score weighting." J Am Stat Assoc 2016; DOI:10.1080/01621459.2016.1260466.
#' @references Li L and Greene T. "A weighting analogue to pair matching in propensity score analysis." Int J Biostat 2013; 9(2):215-234.
#'
#' @seealso \link{psw.balance}, \link{psw.spec.test}, \link{psw.wt}, \link{psw.dr}, \link{psw.mirror.hist}.
#' @import stats
#' @export
psw <- function( data, form.ps, weight, std.diff = FALSE, mirror.hist = FALSE, add.weight = FALSE, nclass = 50,
                 wt = FALSE, out.var = NULL, family = "gaussian",
                 dr = FALSE, form.outcome = NULL,
                 spec.test = F, V.name = NULL, trans.type = NULL,
                 K = 4 ) {
  # Main function to peform propoensity weighting analysis
  #
  # Args
  #   data: data set
  #   form.ps: formula for propensity score model
  #   form.outcome: formula for outcome model
  #   out.var: outcome variable
  #   weight: propensity score weighting method
  #   mirror.hist: plot mirror histogram for propensity score distribution, mirror.hist = F by default
  #   nclass: number of breaks in mirror histogram
  #   std.diff: calculate standardized mean difference? TRUE by default
  #   wt: estimate weighted estimator? TRUE by defaulst
  #   dr: estimate doubly robust estimator? FALSE by default
  #   spec.test: peform propensity score model specification test? spec.test = F by default
  #   V.name: names of covariates that are used for propensity score model specification test
  #   trans.type: tranformation type for propensity score model specification test
  #   K: coefficient of trapezoidal weight, K is the slope of trapezoidal edge, K = 4 by default
  #
  # Return
  #   A list composed of point estimation (est), standard error (std), supplied proprnsity model, propensity score coefficients,
  #   fitted propensity score, propensity score weights.

  # Outcome family can only be "gaussian" or "binomial"
  if ( !( family %in% c( "gaussian", "binomial" ) ) ) {
    stop( "Family should be gaussian or binomial." );
  }

  # check if data is a data frame
  if ( !( is.data.frame( data ) ) ) {
    stop( "Input data must be a data frame." );
  }

  n <- nrow(data);  # number of sujects
  out.ps <- ps.model( dat = data, form.ps = form.ps );

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

  ## standardized mean difference and plot
  if ( std.diff ) {

    #res.balance <- psw.balance.core( Xmat = as.matrix( data[ , Xname, drop=F ] ), Xname = Xname,
    #                                 weight = weight, W = W, Z = as.numeric( data[ , trt.var ] ) );
    res.balance <- psw.balance.core( Xmat = as.matrix( data[ , V.name, drop=F ] ), Xname = V.name,
                                     weight = weight, W = W, Z = as.numeric( data[ , trt.var ] ) );
    res$std.diff.before <- res.balance$std.diff.before;
    res$std.diff.after <- res.balance$std.diff.after;
  }

  ## mirror histogram
  if ( mirror.hist ) {
    if ( weight %in% c( "ATE", "ATT", "ATC" ) ) {
      add.weight <- FALSE;
    }
    mirror.hist.core( ps.above = ps.hat[ data[ , trt.var] == 0 ],
                      ps.below = ps.hat[ data[ , trt.var] == 1 ],
                      wt.above = W[ data[ , trt.var] == 0 ],
                      wt.below = W[ data[ , trt.var] == 1 ],
                      add.weight = add.weight,
                      label.above = paste( trt.var, "=0", sep=""),
                      label.below = paste( trt.var, "=1", sep=""),
                      nclass = nclass );
  }

  ## weighted estimator
  if ( wt ) {
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
  }

  ## doubly robust estimator
  if ( dr ) {
    # fit outcome model
    out.outcome <- outcome.model( dat = data, form = as.formula(form.outcome), trt.var = trt.var, family=family );
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

      # risk ratio
      res$est.rr.dr <- res.dr$est.rr;
      res$std.rr.dr <- res.dr$std.rr;

      # odds ratio
      res$est.or.dr <- res.dr$est.or;
      res$std.or.dr <- res.dr$std.or;

      # log odds ratio
      res$est.lor.dr <- res.dr$est.lor;
      res$std.lor.dr <- res.dr$std.lor;
    }
  }

  ## propensity score model specification test
  if ( spec.test ) {

    if ( weight == "OVERLAP" ) {
      stop( "Specification test is not available for OVERLAP weight." );
    }

    X.mat <- as.matrix( cbind( rep( 1, n ), data[ , Xname, drop=F ] ) ) ;
    V.mat <- as.matrix( data[ , V.name, drop=F ] ) ;

    res.spec.test <- psw.spec.test.core( X.mat = X.mat, V.mat = V.mat, V.name =V.name, trt = data[ , trt.var ], beta.hat = beta.hat,
                                         omega = omega,Q = Q, trans.type = trans.type, weight = weight, delta = 0.002, K = K );
    res$V.name <- res.spec.test$V.name;
    res$g.B1.hat <- res.spec.test$g.B1.hat;
    res$g.B0.hat <- res.spec.test$g.B0.hat;
    res$B.hat <- res.spec.test$B.hat;
    res$var.B.hat <- res.spec.test$var.B.hat;
    res$test.stat <- res.spec.test$test.stat;
    res$df <- res.spec.test$df;
    res$pvalue <- res.spec.test$pvalue;
  }
  return( res );
}
