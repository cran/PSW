######################################################################################
##
## Internal functions, not intended for users
##
######################################################################################
#library( Matrix );  # bdiag()

ps.model <- function(dat, form.ps) {
  # Fit propensity score model
  #
  # Args:
  #   dat: the data from which estimated ps to be calculated
  #   form.ps: propensity score model
  #
  # Return:
  #   Estimated propensity score and object returned by glm fitting

  fm <- glm( formula = as.formula( form.ps ), data = dat, family = binomial(link = "logit"));
  ps.hat <- as.numeric(predict(fm, newdata = dat, type = "response"));
  ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999);  # ps.hat cannot be exactly 0 or 1

  return( list( ps.hat = ps.hat, fm = fm ) );
}

psw.balance.core <- function( Xmat, Xname, weight, W, Z ) {
  # check and plot covariate balance before and after weighting
  #
  # Args:
  #   Xmat: matrix of covariate value
  #   Xname: name of covariates whose balance is to be checked
  #   weight: propensity score weighting method
  #   W: propensity score weight
  #   Z: treatment group
  #
  # Returns:
  #   Two data frames of covariate balance before and after applying weighting and a standardized difference plot

  out1 <- NULL ;
  for ( j in 1 : length( Xname ) ) {
    this.name <- Xname[j];
    this.var <- as.numeric( Xmat[ , this.name ] );
    tmp <- calc.std.diff( var.value = this.var, wt = rep(1, length = nrow( Xmat ) ), Z = Z );
    out1 <- rbind(out1, tmp );
  }
  rownames(out1) <- Xname;
  colnames(out1) <- c("treated.mean","treated.sd","control.mean","control.sd","std.diff.pct");
  std.diff.before <- out1;

  out1 <- NULL;
  for ( j in 1 : length( Xname ) ) {
    this.name <- Xname[j];
    this.var <- as.numeric( Xmat[ , this.name ] );
    tmp <- calc.std.diff( var.value = this.var, wt = W, Z = Z);
    out1 <- rbind(out1, tmp);
  }
  rownames(out1) <- Xname;
  colnames(out1) <- c("treated.mean","treated.sd","control.mean","control.sd","std.diff.pct");
  std.diff.after <- out1;

  # plot difference before and after weighting
  diff.plot( diff.before = std.diff.before[ , "std.diff.pct" ], diff.after = std.diff.after[ , "std.diff.pct" ], name = Xname, weight = weight );

  return( list( std.diff.before = std.diff.before,
                std.diff.after = std.diff.after ) );
}

psw.wt.core <- function( dat, beta.hat, omega, Q, trt.var, out.var, family, Xname, weight, delta=0.002, K=4 ) {
  # Core function to perform propensity score weighting for weighted estimator
  #
  # Args
  #   dat: data set
  #   beta.hat: propensity score model coefficients
  #   omega: omega() function
  #   Q: denominator in the formula of calculating propensity score weight
  #   trt.var: treatment indicator, must be numerically 0 or 1
  #   out.var: outcome variable
  #   family: outcome family, "gaussian" or "binomial"
  #   Xname: covariate name in propensity score model
  #   weight: weighting method
  #   delta: closeness to non-differential point in omega function, delta=0.002 by default
  #   K: coefficient of trapezoidal weight, K is the slope of left trapezoidal edge, K=4 by default
  #
  # Return
  #   A list composed of point estimation (est), standard deviation (std).

  n <- nrow(dat);  # number of sujects

  ## point estimation
  mu1.hat <- sum( ( omega / Q ) * dat[ , trt.var ] * dat[ , out.var ] ) / sum( ( omega / Q ) * dat[ , trt.var ] );
  mu0.hat <- sum( ( omega / Q ) * ( 1 - dat[ , trt.var ]) * dat[ , out.var ] ) / sum( ( omega / Q ) * ( 1 - dat[ , trt.var ] ) );

  if ( family == "gaussian" ) {
    est <- mu1.hat - mu0.hat;
  }
  if ( family == "binomial" ) {
    est.risk <- mu1.hat - mu0.hat;
    est.rr <- mu1.hat / mu0.hat;
    est.or <- ( mu1.hat / ( 1 - mu1.hat ) ) / ( mu0.hat / ( 1 - mu0.hat ) );
    est.lor <- log( est.or );
  }


  ## standard deviation estimation
  Amat <- Bmat <- 0;  # An and Bn matrix for variance calculation
  for (i in 1 : n) {
    Xi <- as.numeric( c( 1, dat[ i, Xname ] ) );
    Zi <- dat[ i, trt.var ];
    Yi <- dat[ i, out.var ];

    ei <- calc.ps.Xbeta( Xmat = Xi, beta = beta.hat );
    ei.deriv1 <- calc.ps.deriv1( Xmat = Xi, beta = beta.hat );
    ei.deriv2 <- calc.ps.deriv2( Xi = Xi, beta = beta.hat );
    omega.ei <- omega[ i ];
    omegaei.deriv <- omega.derive.ei( ps = ei, weight = weight, delta = delta, K=K );
    Qi <- Q[i];
    Qi.deriv <- 2*Zi - 1;

    phi <- c( Zi * ( Yi - mu1.hat ) * omega.ei / Qi,
              ( 1 - Zi ) * ( Yi - mu0.hat ) * omega.ei / Qi,
              ( Zi -ei ) / ( ei * ( 1 - ei ) ) * ei.deriv1 );

    Bmat <- Bmat + outer( phi, phi );  # Bn matrix

    # first row of phi's first derivative w.r.t theta
    first.row <- c( - Zi * omega.ei / Qi, 0, Zi * ( Yi - mu1.hat ) * ei.deriv1 * ( Qi * omegaei.deriv - omega.ei * Qi.deriv ) / Qi^2 );

    # second row of phi's first derivative w.r.t theta
    second.row <- c( 0, - ( 1 - Zi ) * omega.ei / Qi, ( 1 - Zi ) * ( Yi - mu0.hat ) * ei.deriv1 * ( Qi * omegaei.deriv - omega.ei * Qi.deriv ) / Qi^2 );

    # third row of phi's first derivative w.r.t theta
    tmp0 <- matrix( 0, nrow = length( beta.hat ), ncol = 2 );
    tmp1 <- - ei * ( 1 - ei ) * colVec( Xi ) %*% Xi;
    third.row <- cbind( tmp0, tmp1 );

    phi.deriv <- rbind( first.row, second.row, third.row );
    Amat <- Amat + phi.deriv;
  }

  Amat <- Amat/n;
  Bmat <- Bmat/n;
  Amat.inv <- solve( Amat );
  var.mat <- ( Amat.inv %*% Bmat %*% t( Amat.inv ) ) / n;
  #tmp <- c(1, -1, rep(0, length(beta.hat)));
  var.mat <- var.mat[ c( 1 : 2 ), c( 1 : 2 ) ];  # var-covar matrix for point estimators of invididual groups

  # continuous outcome
  if ( family == "gaussian" ) {
    tmp1 <- c( 1, -1 );
    var.est <- rowVec( tmp1 ) %*% var.mat %*% colVec( tmp1 );
    std <- sqrt( as.numeric( var.est ) );
    ans <- list( est = est, std = std );
  }

  # binary outcome
  if ( family == "binomial" ) {
    # risk difference
    tmp.risk <- c( 1, -1 );
    var.risk <- rowVec( tmp.risk ) %*% var.mat %*% colVec( tmp.risk );
    std.risk <- sqrt( as.numeric( var.risk ) );

    # risk ratio
    tmp.rr <- c( 1 / mu1.hat, - mu1.hat / ( mu0.hat^2 ) );
    var.rr <- rowVec( tmp.rr ) %*% var.mat %*% colVec( tmp.rr );
    std.rr <- sqrt( as.numeric( var.rr ) );

    # odds ratio
    tmp1 <- ( 1 - mu0.hat ) / ( mu0.hat * ( 1 - mu1.hat )^2 );
    tmp0 <- - mu1.hat / ( ( 1 - mu1.hat ) * mu0.hat^2 );
    tmp.or <- c( tmp1, tmp0 );
    var.or <- rowVec( tmp.or ) %*% var.mat %*% colVec( tmp.or );
    std.or <- sqrt( as.numeric( var.or ) );

    # log odds ratio
    tmp1 <- 1/( mu1.hat * ( 1 - mu1.hat ) );
    tmp0 <- - 1/( mu0.hat * ( 1 - mu0.hat ) );
    tmp.lor <- c( tmp1, tmp0 );
    var.lor <- rowVec( tmp.lor ) %*% var.mat %*% colVec( tmp.lor );
    std.lor <- sqrt( as.numeric( var.lor ) );

    ans <- list( est.risk = est.risk, std.risk = std.risk,
                 est.rr = est.rr, std.rr = std.rr,
                 est.or = est.or, std.or = std.or,
                 est.lor = est.lor, std.lor = std.lor );
  }
  return( ans );
}

psw.dr.core <- function(dat, beta.hat, omega, Q, out.ps, out.outcome, trt.var, out.var, weight, family, delta=0.002, K=4){
  # Core function for double robust propensity score weighting analysis
  #
  # Args:
  #   dat: date to be used.
  #   out.ps: the object returned by fitting the propensity score model with \code{glm()}.
  #   out.outcome: the object returned by fitting the outcome model with \code{lm()}.
  #   trt.var: dichotomous treatment indicator, \code{1} indicates treated and \code{0} indicates control.
  #   out.var: outcome variable.
  #   weight: propensity score weighting method.
  #   family: outcome family
  #   delta: closeness to non-differential point in omega function, delta=0.002 by default
  #   K: coefficient of trapezoidal weight, K is the slope of left trapezoidal edge, K=4 by default
  # Return
  #   A list of point estimation and standard deviation estimation for the double roubst estimator

  n <- nrow( dat );
  W <- omega/Q;  # generic weight;

  tmp1 <- W*dat[ , trt.var ];
  mu1.hat <- sum( tmp1*(dat[ , out.var ] - out.outcome$Y.hat.m1) )/sum( tmp1 );
  mu2.hat <- sum( omega*out.outcome$Y.hat.m1 )/sum( omega );

  tmp0 <- W*( 1 - dat[ , trt.var ] );
  mu3.hat <- sum( tmp0*(dat[ , out.var ] - out.outcome$Y.hat.m0) )/sum( tmp0 );
  mu4.hat <- sum( omega*out.outcome$Y.hat.m0 )/sum( omega );

  if ( family == "gaussian" ) {
    est <- mu1.hat + mu2.hat - mu3.hat - mu4.hat;
  }
  if ( family == "binomial" ) {
    p1.hat <- mu1.hat + mu2.hat;
    p0.hat <- mu3.hat + mu4.hat;
    est.risk <- p1.hat - p0.hat;
    est.rr <- p1.hat/p0.hat;
    est.or <- (p1.hat/(1-p1.hat)) / (p0.hat/(1-p0.hat));
    est.lor <- log( est.or );
  }

  alpha1.hat <- as.numeric( coef(out.outcome$fm1) );
  alpha0.hat <- as.numeric( coef(out.outcome$fm0) );
  beta.hat <- as.numeric( coef(out.ps$fm) );
  n.alpha1 <- length( alpha1.hat );
  n.alpha0 <- length( alpha0.hat );
  n.beta <- length( beta.hat );

  # sandwich variance
  n <- nrow(dat);
  Amat <- Bmat <- 0;
  for (i in 1:n) {
    Xi <- as.numeric( c(1, dat[i, names(coef(out.ps$fm))[-1] ] ) );  # variable in the propensity score model
    Vi <- as.numeric( c(1, dat[i, names(coef(out.outcome$fm1))[-1] ] ) );  # variable in outcome model
    Zi <- dat[ i, trt.var ] ;
    Yi <- dat[ i, out.var ] ;
    Wi <- W[i];
    Qi <- Q[i];
    omegai <- omega[i];

    ei <- calc.ps.Xbeta(Xmat=Xi, beta=beta.hat);
    ei.deriv1 <- calc.ps.deriv1(Xmat=Xi, beta=beta.hat);
    ei.deriv2 <- calc.ps.deriv2( Xi=Xi, beta=beta.hat);

    Wi.deriv.beta <- calc.W.derive.beta( Zi=Zi, Xi=Xi, omega.ei=omegai, beta.hat=beta.hat, ei=ei, Qi=Qi, weight=weight, delta=delta, K=K );
    omegai.deriv.beta <- calc.omega.derive.beta( Xi=Xi, beta.hat=beta.hat, ei=ei, weight=weight, delta=delta, K=K );

    if ( family == "gaussian" ) {
      m1.hat <- sum(Vi*alpha1.hat);
      m0.hat <- sum(Vi*alpha0.hat);

      m1.deriv.alpha1 <- Vi;
      m0.deriv.alpha0 <- Vi;

      s1.deriv.alpha1 <- -outer(Vi, Vi);
      s0.deriv.alpha0 <- -outer(Vi, Vi);
    }
    if ( family == "binomial" ) {
      tmp1 <- sum(Vi*alpha1.hat);
      tmp0 <- sum(Vi*alpha0.hat);
      m1.hat <- exp(tmp1)/(1+exp(tmp1));
      m0.hat <- exp(tmp0)/(1+exp(tmp0));

      m1.deriv.alpha1 <- m1.hat*(1-m1.hat)*Vi;
      m0.deriv.alpha0 <- m0.hat*(1-m0.hat)*Vi;

      s1.deriv.alpha1 <- -m1.hat*(1-m1.hat)*outer(Vi, Vi);
      s0.deriv.alpha0 <- -m0.hat*(1-m0.hat)*outer(Vi, Vi);
    }

    this.phi.row1 <- Wi*Zi*( Yi - m1.hat - mu1.hat );
    this.phi.row2 <- omegai*( m1.hat - mu2.hat );
    this.phi.row3 <- Wi*(1-Zi)*( Yi - m0.hat - mu3.hat );
    this.phi.row4 <- omegai*( m0.hat - mu4.hat );
    this.phi.row5 <- Zi*( Yi - m1.hat )*Vi;
    this.phi.row6 <- (1-Zi)*( Yi - m0.hat )*Vi;
    this.phi.row7 <- (Zi-ei)/ei/(1-ei)*ei.deriv1;
    this.phi <- c( this.phi.row1, this.phi.row2, this.phi.row3,
                   this.phi.row4, this.phi.row5, this.phi.row6, this.phi.row7 );
    Bmat <- Bmat + outer( this.phi, this.phi );

    quad1 <- diag( c( -Wi*Zi, -omegai, -Wi*(1-Zi), -omegai ) );
    quad2 <- matrix(0, nrow=n.alpha1+n.alpha0+n.beta, ncol=4);

    tmp <- c( -Wi*Zi*m1.deriv.alpha1, rep(0, n.alpha0), Zi*( Yi - m1.hat - mu1.hat )*Wi.deriv.beta,
              omegai*m1.deriv.alpha1, rep(0, n.alpha0), ( m1.hat - mu2.hat )*omegai.deriv.beta,
              rep(0, n.alpha1), -Wi*(1-Zi)*m0.deriv.alpha0, (1-Zi)*( Yi - m0.hat - mu3.hat )*Wi.deriv.beta,
              rep(0, n.alpha1), omegai*m0.deriv.alpha0, ( m0.hat - mu4.hat )*omegai.deriv.beta );
    quad3 <- matrix( tmp, byrow = TRUE, nrow = 4,
                     ncol = n.alpha1 + n.alpha0 + n.beta );

    quad4.blk1 <- cbind( Zi*s1.deriv.alpha1, matrix(0, nrow=n.alpha1, ncol=n.alpha0 + n.beta) );
    quad4.blk2 <- cbind( matrix(0, nrow=n.alpha0, ncol=n.alpha1 ), (1-Zi)*s0.deriv.alpha0, matrix(0, nrow=n.alpha0, ncol=n.beta) );
    quad4.blk3 <- cbind( matrix(0, nrow=n.beta, ncol=n.alpha1+n.alpha0), -ei*(1-ei)*outer(Xi, Xi) );
    quad4 <- rbind( quad4.blk1, quad4.blk2, quad4.blk3 );

    this.phi.deriv <- rbind( cbind(quad1, quad3) ,
                             cbind(quad2, quad4) );
    Amat <- Amat + this.phi.deriv;
  }
  Amat <- Amat/n;
  Bmat <- Bmat/n;
  Amat.inv <- solve(Amat) ;
  var.mat <- ( Amat.inv %*% Bmat %*% t(Amat.inv) )/n;
  #tmp1 <- c(1, 1, -1, rep(0, length = n.alpha1+n.alpha0+n.beta) );
  var.mat <- var.mat[ c(1:4), c(1:4) ];

  if ( family == "gaussian" ) {
    tmp <- c( 1, 1, -1, -1 );
    var.est <- rowVec(tmp) %*% var.mat %*% colVec(tmp);  # get variance of mu1+mu2-mu3;
    std <- sqrt( as.numeric(var.est) );
    ans <- list( est = est, std = std );
  }
  if ( family == "binomial" ) {
    # risk difference
    tmp.risk <- c( 1, 1, -1, -1 );
    var.risk <- rowVec( tmp.risk ) %*% var.mat %*% colVec( tmp.risk );
    std.risk <- sqrt( as.numeric( var.risk ) );

    ans <- list( est.risk = est.risk, std.risk = std.risk );
  }

  return( ans );
}

psw.spec.test.core <- function( X.mat, V.mat, V.name, trt, beta.hat, omega, Q, trans.type, weight, delta=0.002, K=4 ) {
  # check and test covariate balance between treatment and control arms.
  #
  # Args:
  #   X.mat: covarite matrxi for PS model
  #   V.mat: covariate matrxi used for balance test
  #   V.name: variables for balance test
  #   trt: treatment indicator
  #   beta.hat: fitted PS model parameters
  #   omega: omega function
  #   Q: denominator of PS weight function
  #   trans.type: transformation types applied to covariates, including log, logit, Fisher's Z, square root, and identity transformation.
  #   weight: weighting method.
  #   delta: closeness to non-differential point in omega function, delta=0.002 by default
  #   K: coefficient of trapezoidal weight, K is the slope of left trapezoidal edge, K=4 by default
  #
  # Returns:
  #   A list of balance test parameters


  n <- nrow( X.mat );  # number of sujects
  W <- omega/Q;  # generic weight

  # obtain balance test statistic.
  mu.B1.hat <- as.numeric( t(V.mat) %*% colVec( W * trt ) ) / sum( W * trt ) ;
  mu.B0.hat <- as.numeric( t(V.mat) %*% colVec( W * ( 1 - trt ) ) ) / sum( W * ( 1 - trt ) ) ;
  eta.B1.hat <- g.fun( x=mu.B1.hat, fun.type = trans.type ) ;
  eta.B0.hat <- g.fun( x=mu.B0.hat, fun.type = trans.type ) ;

  theta.hat <- c( eta.B1.hat, eta.B0.hat, beta.hat ) ;
  B.hat <- eta.B1.hat - eta.B0.hat ;

  n.mu.B1 <- length( mu.B1.hat ) ;
  n.mu.B0 <- length( mu.B0.hat ) ; # NOTE: n.mu.B1 = n.mu.B1 = length(B.hat)
  n.beta <- length( beta.hat ) ;
  n.theta <- length(theta.hat) ;
  D.mat <- cbind( diag(n.mu.B1), -diag(n.mu.B0),
                  matrix(0, nrow=n.mu.B1, ncol=n.beta) ) ;

  Meat.mat <- Bread.mat <- 0 ; # the meat and bread of the sandwich variance
  for ( i in 1 : n ) {
    Vi <- as.numeric( V.mat[i, ] );
    Xi <- as.numeric( X.mat[i, ] );
    Zi <- trt[i];
    Wi <- W[i];
    omegai <- omega[i];
    Qi <- Q[i];
    ei <- calc.ps.Xbeta(Xmat=Xi, beta=beta.hat);
    ei.deriv <- calc.ps.deriv1( Xmat = Xi, beta = beta.hat );
    omegai.deriv <- omega.derive.ei( ps = ei, weight = weight, delta = delta, K = K );
    Qi.deriv <- 2 * Zi - 1;
    Wi.deriv.beta <- ei.deriv * ( Qi*omegai.deriv - omegai*Qi.deriv ) / Qi^2;

    tmp1 <- Wi * Zi * ( Vi - g.inv( x = eta.B1.hat, fun.type = trans.type ) );
    tmp2 <- Wi * ( 1 - Zi ) * ( Vi - g.inv( x = eta.B0.hat, fun.type = trans.type ) );
    tmp3 <- ( Zi - ei ) * Xi;
    this.phi <- c( tmp1, tmp2, tmp3 );

    Meat.mat <- Meat.mat + outer( this.phi, this.phi );

    Block.11 <- - Wi * Zi *( diag( g.inv.deriv( x = eta.B1.hat, fun.type = trans.type ) ) );
    Block.12 <- matrix( 0, nrow = n.mu.B1, ncol = n.mu.B0 );
    Block.13 <- Zi * ( colVec( Vi - g.inv( x = eta.B1.hat, fun.type = trans.type ) ) %*% rowVec( Wi.deriv.beta ) );
    Block.21 <- matrix( 0, nrow = n.mu.B0, ncol = n.mu.B1 ) ;
    Block.22 <- - Wi * ( 1 - Zi ) * ( diag( g.inv.deriv( x = eta.B0.hat, fun.type = trans.type ) ) );
    Block.23 <- (1-Zi)*( colVec( Vi - g.inv( x = eta.B0.hat, fun.type = trans.type ) ) %*% rowVec( Wi.deriv.beta ) );
    Block.31 <- matrix( 0, nrow = n.beta, ncol = n.mu.B1 );
    Block.32 <- matrix( 0, nrow = n.beta, ncol = n.mu.B0 );
    Block.33 <- -ei * ( 1 - ei ) * outer( Xi, Xi ) ;
    Bread.mat <- Bread.mat + rbind( cbind( Block.11, Block.12, Block.13 ) ,
                                    cbind( Block.21, Block.22, Block.23 ) ,
                                    cbind( Block.31, Block.32, Block.33 ) );
  }
  B.mat <- Meat.mat/n;

  A.mat <- Bread.mat/n;
  A.mat.inv <- solve( A.mat );
  Sigma.theta.hat <- ( A.mat.inv %*% B.mat %*% t( A.mat.inv ) )/n;
  Sigma.B.hat <- D.mat %*% Sigma.theta.hat %*% t( D.mat );
  df <- qr(Sigma.B.hat)$rank;  # degree of freedom
  test.stat <- as.numeric( t(B.hat) %*% solve(Sigma.B.hat) %*% B.hat );
  pvalue <- pchisq( test.stat, df = df, lower.tail = FALSE );
  # chi-squared test for overall comparison

  names(eta.B1.hat) <- paste("eta.B1.", V.name, sep="");
  names(eta.B0.hat) <- paste("eta.B0.", V.name, sep="");
  names(B.hat) <- c( paste("B.hat.", V.name, sep="") )
  rownames(Sigma.B.hat) <- colnames(Sigma.B.hat) <- names(B.hat);
  # Add names to point and variance estimators

  return( list( weight = weight,
                V.name = V.name,
                g.B1.hat = eta.B1.hat,
                g.B0.hat = eta.B0.hat,
                B.hat = B.hat,
                var.B.hat = Sigma.B.hat,
                test.stat = test.stat,
                df = df,
                pvalue = pvalue
                ) );
}

mirror.hist.core <- function( ps.above, ps.below, wt.above, wt.below, add.weight = FALSE,
                              label.above, label.below, nclass = 50 ) {
  # Mirror histogram is used to plot the mirror histogram showing the propensity score distribution overlap between groups.
  #   ps.above propensity score for the control.
  #   ps.below propensity score for the treated.
  #   wt.above weight for the control.
  #   wt.below weight for the treated.
  #   add.weight add to the mirror histogram propensity score weights, \code{add.weight=FALSE} by default.
  #   label.above label for the control group.
  #   label.below label for the treated group.
  #   label.title plot title.
  #   nclass number of breaks, \code{nclass=50} by default.

  x0 <- ps.above ; wt0 <- wt.above ;
  x1 <- ps.below ; wt1 <- wt.below ;
  if ( is.null(nclass) ) {
    breaks <- hist( c(x0, x1), nclass=nclass, plot=F )$breaks ;
  } else {
    breaks <- hist( c(x0, x1), nclass=nclass, plot=F )$breaks ;
  }
  fm.x0 <- hist( x0, breaks=breaks, plot = F ) ;
  fm.x1 <- hist( x1, breaks=breaks, plot = F ) ;
  x.range <- range( c(x0, x1) ) ;
  y.range <- c( -max( fm.x1$counts ), max( fm.x0$counts ) ) ;

  plot( x=x.range, y=y.range, xaxt="n", yaxt="n", type="n" ,
        xlim=x.range, ylim=y.range, ylab="",
        xlab="Propensity score", cex.lab=2 ) ;
  axis( side=1, at=pretty( x.range ), cex.axis=2 ) ;
  axis( side=2, at=pretty( y.range ) , labels=abs( pretty( y.range ) ), cex.axis=2 ) ;
  title( ylab="Frequency", cex.lab=2, line=4 );
  abline( h=0, lty=1 ) ;

  # plot the histogram above the horizontal line
  fm <- fm.x0 ;
  for (i in 1:(length(breaks)-1)) {
    x <- c( fm$breaks[i], fm$breaks[i], fm$breaks[i+1], fm$breaks[i+1] ) ;
    y <- c(0, fm$counts[i], fm$counts[i], 0) ;
    polygon(x, y) ;

    # add weight
    if ( add.weight ) {
      tmp <- ( breaks[i] <= x0 & x0 <= breaks[i+1] ) ;
      y2 <- c(0, sum(wt0[tmp]), sum(wt0[tmp]), 0) ;
      polygon(x, y2, col="light green", border="black") ;
    }
  }
  # plot the histogram below the horizontal line
  fm <- fm.x1 ;
  for (i in 1:(length(breaks)-1)) {
    x <- c( fm$breaks[i], fm$breaks[i], fm$breaks[i+1], fm$breaks[i+1] ) ;
    y <- c(0, -fm$counts[i], -fm$counts[i], 0) ;
    polygon(x, y) ;

    # add weight
    if ( add.weight ) {
      tmp <- ( breaks[i] <= x1 & x1 <= breaks[i+1] ) ;
      y2 <- c(0, -sum(wt1[tmp]), -sum(wt1[tmp]), 0) ;
      polygon(x, y2, col="dark green", border="black") ;
    }
  }
}

########################################################################################
##
## Supporting functions
##
########################################################################################
rowVec <- function(x) {
  return( t(x) );
}

colVec <- function(x) {
  return( t(t(x)) );
}

outcome.model <- function(dat, form, trt.var, family) {
  # outcome.model to fit outcome with specified formula
  #
  # Args:
  #   dat: data based on which outcome model is to be fitted.
  #   form: outcome come model formula.
  #   trt.var treatment indicator
  #   family: outcome family
  #
  # Returns:
  #   Estimated outcome by treatment indicator and lm fitting

  fm1 <- glm( as.formula( form ), data = dat[ dat[ , trt.var ] == 1, ] );
  fm0 <- glm( as.formula( form ), data = dat[ dat[ , trt.var ] == 0, ] );
  Y.hat.m1 <- as.numeric( predict( fm1, newdata = dat ) );
  Y.hat.m0 <- as.numeric( predict( fm0, newdata = dat ) );
  return( list( Y.hat.m1 = Y.hat.m1 ,
                Y.hat.m0 = Y.hat.m0 ,
                fm1 = fm1 , fm0 = fm0 ) );
}

calc.omega <- function(ps, weight, delta=0.002, K=4) {
  # Calculate omega for each weighting method
  #
  # Args:
  #   ps: estimated proopensity score,
  #   weight: weighting method.
  #   delta: closeness to non-differentiable region
  #   K: wighting coefficient for trapezoidal weighting
  #
  # Return:
  #   A vector of length(ps)

  ans <- 0;
  if (weight == "ATE") {
    ans <- 1;
  } else if (weight == "MW") {
    ans <- calc.omega.MW( ps = ps, delta = delta );
  } else if (weight == "ATT") {
    ans <- ps;
  } else if (weight == "ATC") {
    ans <- 1-ps;
  } else if (weight == "OVERLAP") {
    ans <- 4*ps*(1-ps);
  } else if (weight == "TRAPEZOIDAL") {
    ans <- calc.omega.trapzd(ps=ps, delta=delta, K=K);
  } else {
    stop("Error in calc.omega: weight method does not exist!");
  }

  return( ans );
}

calc.ps.Xbeta <- function( Xmat, beta ) {
  # Calculate the propensity score, e(X, beta)
  #
  # Args:
  #   Xmat: a vector or a matrix with each column being X_i
  #   beta: propensity score model coefficients, of the same length as the nrow(X)
  #
  # Return:
  #   A scalar of matching weight

  Xmat <- as.matrix(Xmat);
  tmp <- as.numeric(rowVec(beta) %*% Xmat);
  tmp <- exp(tmp);
  names(tmp) <- NULL;

  return( tmp/(1 + tmp) );
}

calc.ps.deriv1 <- function(Xmat, beta) {
  # Calculate the derivative of propensity score w.r.t. beta.
  #
  # Args:
  #   Xmat: a matrix with each column being X_i,
  #   beta: propensity score coefficients, of the same length as the nrow(X).
  #
  # Return:
  #   A vector of the same dimension as beta.

  Xmat <- as.matrix(Xmat);
  tmp.ps <- calc.ps.Xbeta( Xmat=Xmat, beta=beta );
  ans <- tmp.ps*(1 - tmp.ps)*t(Xmat);  # Xmat is row vector from a single line of data
  names(ans) <- rownames(ans) <- NULL;

  return( t(ans) );
}


calc.ps.deriv2 <- function( Xi, beta ) {
  # Calculate the second derivative of propensity score w.r.t beta.
  #
  # Args:
  #   Xi: a vector (not a matrix!)
  #   beta:  a vector of coefficients, of the same length as Xi
  #
  # Return:
  #   A matrix of secondd derivative

  Xi <- colVec(Xi);
  tmp.ps <- calc.ps.Xbeta( Xmat=Xi, beta=beta );
  tmp.deriv1 <- calc.ps.deriv1( Xmat=Xi, beta=beta );
  ans <- Xi %*% rowVec(tmp.deriv1);
  names(ans) <- rownames(ans) <- NULL;
  ans <- (1 - 2*tmp.ps)*ans;

  return( ans );
}

calc.omega.MW <- function (ps, delta) {
  # Calculate omega value for MW method
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differentiable region
  #
  # Return:
  # Value of omega, scalar

  ans <- 0;
  if (ps <= 0.5 - delta) {
    ans <- 2*ps;
  } else if (ps >= 0.5 + delta) {
    ans <- 2*(1 - ps);
  } else {
    ans <- approx.omega.MW(ps, delta);
  }

  return( ans );
}

approx.omega.MW <- function (ps, delta) {
  # Approximate omega function of MW at non-differentiable region, eta(ps)
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differentiable region
  #
  # Return:
  # Approximated omega at non-differentiable region, scalar

  A <- solve.A.MW(delta);
  ans <- rowVec(A) %*% c(1, ps, ps^2, ps^3);

  ans;
}

solve.A.MW <- function(delta) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of MW
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #
  # Return
  #   A vecotor of solved coefficients

  if ( delta < 0.00001 ) {
    stop("*** ERROR in solve.a: delta too small ***");
  }

  tmp1 <- 0.5 - delta;
  tmp2 <- 0.5 + delta;

  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE);
  C <- 2*c(tmp1, 1, tmp1, -1);
  A <- solve(D) %*% C;  # coefficients of cubic polynomial

  A;
}

calc.omega.trapzd <- function (ps, delta, K) {
  # Calculate omega value for trapezoidal weighting method
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differentiable region
  # K: trapezoidal weighting coffecient
  #
  # Return:
  # Value of omega, scalar

  ans <- 0;
  if ( (0 < ps) & (ps <= 1/K - delta) ) {
    ans <- K*ps;
  } else if ( (1/K + delta <= ps) & (ps <= 1 - 1/K - delta) ) {
    ans <- 1;
  } else if ( (1 - 1/K + delta <= ps) & (ps < 1) ) {
    ans <- K*(1 - ps);
  } else {
    ans <- approx.omega.trapzd(ps, delta, K);
  }

  ans;
}

approx.omega.trapzd <- function (ps, delta, K) {
  # Approximate omega function of trapezoidal weight at non-differentiable region, eta(ps)
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differentiable region
  # K: trapezoidal weight coefficient
  #
  # Return:
  # Approximated omega at non-differentiable region, scalar

  A <- 0;
  if ( (1/K - delta < ps) & (ps < 1/K + delta) ) {
    A <- solve.A.trapzd1st(delta=delta, K=K);
  } else {
    A <- solve.A.trapzd2nd(delta=delta, K=K);
  }
  ans <- rowVec(A) %*% c(1, ps, ps^2, ps^3);

  ans;
}

solve.A.trapzd1st <- function (delta, K) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of
  # trapezoidal weighting at the first non-differentiable pivot
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #   K: coefficient of trapezoidal weight
  #
  # Return
  #   A vecotor solved coefficients

  if ( delta < 0.00001 ) {
    stop("*** ERROR in solve.a: delta too small ***");
  }

  tmp1 <- 1/K - delta;
  tmp2 <- 1/K + delta;

  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE);
  C <- 2*c(K*tmp1, K, 1, 0);
  A <- solve(D) %*% C;  # coefficients of cubic polynomial

  A;
}

solve.A.trapzd2nd <- function (delta, K) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of
  # trapezoidal weighting at the second non-differentiable pivot
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #   K: coefficient of trapezoidal weight
  #
  # Return
  #   A vecotor solved coefficients

  if ( delta < 0.00001 ) {
    stop("*** ERROR in solve.a: delta too small ***");
  }

  tmp1 <- 1 - 1/K - delta;
  tmp2 <- 1 - 1/K + delta;

  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE);
  C <- 2*c(1, 0, K*(1/K - delta), -K);
  A <- solve(D) %*% C;  # coefficients of cubic polynomial

  A;
}

omega.derive.ei <- function(ps, weight, delta=0.002, K=4) {
  # Calculate the first derivative of omega function w.r.t propensity score (ei)
  #
  # Args:
  #    ps: estimated propensity score
  #    weight: selected weighting method
  #    delta: pre-defined closeness to non-differentiable regions
  #    K: weighting coefficient for trapezoidal weighting
  #
  # Returns:
  #    The first derivative of omega w.r.t to propensity score

  if (weight == "ATE") {
    ans <- 0;
  } else if (weight == "ATT") {
    ans <- 1;
  } else if (weight == "ATC") {
    ans <- -1;
  } else if (weight == "OVERLAP") {
    ans <- 4*(1 - 2*ps);
  } else if (weight == "MW") {
    ans <- omega.derive.ei.MW(ps, delta);
  } else if (weight == "TRAPEZOIDAL") {
    ans <- omega.derive.ei.trapzd(ps, delta, K);
  } else {
    stop( "User defined first-order derivative of omega function is not provided!" );
  }

  ans;
}

omega.derive.ei.MW <- function (ps, delta) {
  # calculate of the first derivative of omega function w.r.t propensity score with approximation
  # in MW method
  #
  # Args:
  #    ps: propensity score
  #    delta: pre-defined closeness to non-differentiable regions
  #
  # Returns:
  #    The first derivative of omega (MW method) w.r.t to propensity score,
  #    approximation at non-differentiable region is included.

  if ( (0 < ps) & (ps <= 0.5 - delta) ) {
    ans <- 2;
  } else if ( (0.5 + delta <= ps) & (ps <1) ) {
    ans <- -2;
  } else {
    A <- solve.A.MW(delta);
    ans <- A[2] + 2*A[3]*ps + 3*A[4]*ps^2;
  }

  ans;
}

omega.derive.ei.trapzd <- function (ps, delta, K) {
  # calculate of the first derivative of omega function w.r.t propensity score with approximation
  # in Trapezoidal weighting method
  #
  # Args:
  #    ps: propensity score
  #    delta: pre-defined closeness to non-differentiable regions
  #    K: coefficient of trapezoidal weighting
  #
  # Returns:
  #    The first derivative of omega (trapezoidal weighting method) w.r.t to propensity score,
  #    approximation at non-differentiable region is included.

  if ( (0 < ps) & (ps <= 1/K - delta) ) {
    ans <- K;
  } else if ( (1/K - delta < ps) & (ps < 1/K + delta) ) {
    A <- solve.A.trapzd1st(delta = delta, K = K);
    ans <- A[2] + 2*A[3]*ps + 3*A[4]*ps^2;
  } else if ( (1/K + delta <= ps) & (ps <= 1 - 1/K - delta) ) {
    ans <- 0;
  } else if ( (1 - 1/K + delta <= ps) & (ps < 1) ) {
    ans <- -K;
  } else {
    A <- solve.A.trapzd2nd(delta = delta, K = K);
    ans <- A[2] + 2*A[3]*ps + 3*A[4]*ps^2;
  }

  ans;
}

calc.W.derive.beta <- function(Zi, Xi, omega.ei, beta.hat, ei, Qi, weight, delta, K ) {
  # Calculate the first order derivative of weight (W) with respect to beta (coefficients in ps model)
  #
  # Args:
  #   Zi: treamtment indicator
  #   Xi: covariates of i-th subject
  #   omega.ei: value of omega for subject
  #   beta.hat: parameter estimation in ps model
  #   ei: propensity score
  #   Qi: denominator of weight formula
  #   weight: weight type
  #   delta: closeness to non-differential points
  #   K: slope of trapezoidal left edge
  #
  # Returns
  #   A vector, first order derivative of weight w.r.t. beta

  ei.deriv1 <- calc.ps.deriv1( Xmat=Xi, beta=beta.hat );
  omegaei.deriv <- omega.derive.ei( ps = ei, weight = weight, delta = delta, K = K );
  Qi.deriv.ei <- 2*Zi - 1;

  ans <- ei.deriv1*(Qi*omegaei.deriv - omega.ei*Qi.deriv.ei)/Qi^2;

  return( ans );
}

calc.omega.derive.beta <- function(Xi, beta.hat, ei, weight, delta=0.002, K=4 ) {
  # Calculate the first order derivative of omega with respect to beta (coefficients in ps model)
  #
  # Args:
  #   Xi: covariates of i-th subject
  #   beta.hat: parameter estimation in ps model
  #   ei: propensity score of i-th subject
  #   weight: weight type
  #   delta: closeness to non-differential points, which is 0.002 by default
  #   K: slope of trapezoidal left edge, K=4 by default
  #
  # Returns
  #   A vector, derivative of omega function w.r.t. beta

  ei.deriv1 <- calc.ps.deriv1( Xmat=Xi, beta=beta.hat );
  omegaei.deriv <- omega.derive.ei( ps=ei, weight=weight, delta=delta, K=K );
  ans <- omegaei.deriv * ei.deriv1;

  return( ans );
}

#' @import Hmisc
calc.std.diff <- function(var.value, wt, Z) {
  # calculate standardized difference between treatment arms
  #
  # Args:
  #   var.value: covariate value.
  #   wt: weights of each subject.
  #   Z: treatment indicator
  #
  #   Returns:
  #   Absolute standardized difference between treatment groups in 100% scale

  Z1.mean <- wtd.mean( x = var.value[Z==1], weights = wt[Z==1] );
  Z1.var <- wtd.var( x = var.value[Z==1], weights = wt[Z==1] );
  Z0.mean <- wtd.mean( x = var.value[Z==0], weights = wt[Z==0] );
  Z0.var <- wtd.var( x = var.value[Z==0], weights = wt[Z==0] );
  std.diff <- 100 * ( Z1.mean - Z0.mean ) / sqrt( ( Z1.var + Z0.var ) / 2 ) ;
  ans <- c( Z1.mean, sqrt(Z1.var), Z0.mean, sqrt(Z0.var), std.diff );

  return( ans );
}

diff.plot <- function( diff.before, diff.after, name, weight ) {
  # plot standardize difference for before and after weighting
  #
  # Args:
  #   diff.before: standardized mean differece before weighting
  #   diff.after: standardized mean difference after weighting
  #   name: covariate name
  #   weight: propensity score weighting method
  # Returns
  #   A plot is generated

  par( las=1, lwd = 2, mar=c(5, max( nchar(name) ), 4, 2), bty="n" );

  x.range <- range( c(diff.before, diff.after) );
  y.range <- c(1, length(name));
  ord <- order( diff.before, decreasing = T );
  plot( x=x.range, y=y.range, xaxt="n", yaxt="n", type="n",
        xlim=x.range, ylim=y.range, ylab="", xlab="Standardized difference" );
  axis( side=1, at=pretty( x.range ) );
  axis( side=2, at=length(name):1, labels=name[ord], tick=F );
  abline( v=0, lty=1, col="gray" );
  points( y = length(name):1, x = diff.before[ord], pch=4 );
  points( y = length(name):1, x = diff.after[ord], pch=21 );
  legend("topleft", legend=c("Unadjusted", weight), pch=c(4, 21) );
}

#' @import gtools
g.fun <- function( x, fun.type ) {
  # Perform monotone transformation for vector of covariates
  # Included transformation types are: log, logit, square root, Fisher's Z transformation
  # and identity transformation (default).
  #
  # Args:
  #   x: vector of covariates
  #   fun.type: transformation type
  #
  # Returns
  #   A vector of transformed values

  if ( length(x) != length(fun.type) ) {
    print("*** ERROR in g.fun() ***");
    return( rep(NA, length(x)) );
  } else {
    m <- length(x);
    ans <- rep(0, m);
    for (j in 1:m) {
      if ( fun.type[j] == "log" ) {
        ans[j] <- log( max( 1e-6, x[j] ) );
      } else if ( fun.type[j] == "logit" ) {
        ans[j] <- logit( min( max( 1e-6, x[j]), 1-1e-6 ) );
      } else if ( fun.type[j] == "Fisher" ) {
        ans[j] <- 0.5*( log( 1 + x[j] ) - log( 1 - x[j] ) );
      } else if (fun.type[j] == "sqrt") {
        ans[j] <- sqrt( x[j] );
      } else {  # identity transformation, by default
        ans[j] <- x[j] ;
      }
    }
    return( ans ) ;
  }
}

#' @import gtools
g.inv <- function( x, fun.type ) {
  # Calculate the inverse function of transformed covaraites
  # Inverse back from log, logit, square root, and Fisher's Z transformations.
  #
  # Args:
  #   x: vector of previously transformed covariates
  #   fun.type: transformation type
  #
  # Returns
  #   A vector of inverse function value.

  if ( length(x) != length(fun.type) ) {
    print("*** ERROR in g.inv() ***") ;
    return( rep(NA, length(x)) ) ;
  } else {
    m <- length(x) ;
    ans <- rep(0, m) ;
    for (j in 1:m) {
      if ( fun.type[j] == "log" ) {
        ans[j] <- exp( x[j] ) ;
      } else if ( fun.type[j] == "logit" ) {
        ans[j] <- inv.logit( x[j] ) ;
      } else if ( fun.type[j] == "Fisher" ) {
        ans[j] <- ( (exp(2*x[j])-1)/(exp(2*x)+1) );
      } else if ( fun.type[j] == "sqrt") {
        ans[j] <- x[j]^2;
      } else {  # identity transformation, by default
        ans[j] <- x[j] ;
      }
    }
    return( ans ) ;
  }
}

#' @import gtools
g.inv.deriv <- function( x, fun.type ) {
  # Calculate the first order derivative of inverse function for transformed covaraites
  #
  # Args:
  #   x: vector of previously transformed covariates
  #   fun.type: type of transformation
  #
  # Returns
  #   A vector of first order derivative.

  if ( length(x) != length(fun.type) ) {
    print("*** ERROR in g.inv.deriv() ***") ;
    return( rep(NA, length(x)) ) ;
  } else {
    m <- length(x) ;
    ans <- rep(0, m) ;
    for (j in 1:m) {
      if ( fun.type[j] == "log" ) {
        ans[j] <- exp( x[j] ) ;
      } else if ( fun.type[j] == "logit" ) {
        tmp <- inv.logit( x[j] ) ;
        ans[j] <- tmp*(1-tmp) ;
      } else if ( fun.type[j] == "Fisher" ) {
        tmp <- 2*exp( 2*x[j] );
        ans[j] <- tmp^2;
      } else if ( fun.type[j] == "sqrt" ) {
        ans[j] <- 2*x[j];
      } else {  # identity transformation, by default
        ans[j] <- 1 ;
      }
    }
    return( ans ) ;
  }
}



