#' @name EXTLindley
#' @aliases EXTLindley dextlindley pextlindley qextlindley rextlindley hextlindley
#'
#' @title Extended Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the extended Lindley distribution with parameters theta, alpha and beta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @note The \code{\link[stats]{uniroot}} function with default arguments is used to find out the quantiles.
#'
#' @references
#'
#' Bakouch, H. S., Al-Zahrani, B. M., Al-Shomrani, A. A., Marchi, V. A. A., Louzada, F., (2012). An extended Lindley distribution. \emph{Journal of the Korean Statistical Society}, \bold{41}, (1), 75-85.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta positive parameter.
#' @param beta greater than or equal to zero.
#' @param alpha \eqn{\rm I\!R^{-}\cup (0,1)}.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param L,U interval which \code{uniroot} searches for a root (quantile), L = 1e-4 and U = 50 are the default values.
#'
#' @return \code{dextlindley} gives the density, \code{pextlindley} gives the distribution function, \code{qextlindley} gives the quantile function, \code{rextlindley} generates random deviates and \code{hextlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[lamW]{lambertWm1}}, \code{\link[stats]{uniroot}}.
#'
#' @source [d-h-p-q-r]extlindley are calculated directly from the definitions. \code{rextlindley} uses the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta ,\alpha ,\beta )=\frac{\theta }{\left( 1+\theta \right) }\left( 1+\frac{\theta x}{1+\theta }\right) ^{\alpha -1}\left[ \beta \left(  1+\theta +\theta x\right) \left( \theta x\right) ^{\beta -1}-\alpha \right]e^{-\left( \theta x\right) ^{\beta }}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ,\alpha ,\beta )=1-\left( 1+\frac{\theta x}{1+\theta }\right)^{\alpha }e^{-\left( \theta x\right) ^{\beta }}}
#'
#' Quantile function
#' \deqn{\code{does not have a closed mathematical expression}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha ,\beta )=\frac{\beta \left( 1+\theta +\theta x\right)\theta ^{\beta }x^{\beta -1}-\alpha \theta }{\left( 1+\theta +\theta x\right) }}
#'
#' \bold{Particular cases:} \eqn{(\alpha = 1, \beta = 1)} the one-parameter Lindley distribution, \eqn{(\alpha = 0, \beta = 1)} the exponential distribution and for \eqn{\alpha = 0} the Weibull distribution. See Bakouch et al. (2012) for other particular cases.
#'
#' @examples
#' set.seed(1)
#' x <- rextlindley(n = 1000, theta = 5.0, alpha = -1.0, beta = 5.0)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' plot(S, dextlindley(S, theta = 5.0, alpha = -1.0, beta = 5.0), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pextlindley(q, theta = 5.0, alpha = -1.0, beta = 5.0, lower.tail = TRUE)
#' pextlindley(q, theta = 5.0, alpha = -1.0, beta = 5.0, lower.tail = FALSE)
#' qextlindley(p, theta = 5.0, alpha = -1.0, beta = 5.0, lower.tail = TRUE)
#' qextlindley(p, theta = 5.0, alpha = -1.0, beta = 5.0, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'extlindley', start = list(theta = 5.0, alpha = -1.0, beta = 5.0))
#' plot(fit)
#'
#
#' @rdname EXTLindley
#' @export
dextlindley <- function(x, theta, alpha, beta, log = FALSE)
{
  if(!((alpha == 1) || (alpha <= 0))) stop('alpha must be equal to 1 or less than or equal to 0')
  stopifnot(theta > 0, beta >= 0)
  if(log)
  {
	t3 <- log1p(x * theta + theta)
	t5 <- log1p(theta)
	t8 <- x ^ beta
	t9 <- theta ^ beta
	t10 <- t9 * t8
	t11 <- log(x)
	t12 <- beta + 1
	t13 <- x ^ t12
	t14 <- theta ^ t12
	t23 <- log(-alpha * theta * x + beta * t14 * t13 + beta * t8 * t14 + beta * t10)
	(t3 - t5) * alpha - t10 - t3 - t11 + t23
 }
 else
 {
    t1 <- x * theta
    t2 <- t1 + theta + 1
    t6 <- (1 / (1 + theta) * t2) ^ alpha
    t7 <- t1 ^ beta
    t8 <- exp(-t7)
    t10 <- beta * t7
    1 / x / t2 * (-alpha * theta * x + t1 * t10 + theta * t10 + t10) * t8 * t6
  }
}

#' @rdname EXTLindley
#' @export
pextlindley <- function(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  if(!((alpha == 1) || (alpha <= 0))) stop('alpha must be equal to 1 or less than or equal to 0')
  stopifnot(theta > 0, beta >= 0)
  if(lower.tail)
  {
    t1 <- theta * q
    t6 <- (1 + 1 / (1 + theta) * t1) ^ alpha
    t7 <- t1 ^ beta
    t8 <- exp(-t7)
    cdf<- -t8 * t6 + 1
  }
  else
  {
    t1 <- theta * q
    t6 <- (1 + 1 / (1 + theta) * t1) ^ alpha
    t7 <- t1 ^ beta
    t8 <- exp(-t7)
    cdf<- t8 * t6
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname EXTLindley
#' @export
qextlindley <- function(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE, L = 1e-4, U = 50)
{
  if(!((alpha == 1) || (alpha <= 0))) stop('alpha must be equal to 1 or less than or equal to 0')
  stopifnot(theta > 0, beta >= 0)
  if(lower.tail)
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pextlindley(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
  }
  else
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pextlindley(q, theta, alpha, beta, lower.tail = FALSE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
    if(log.p) return(log(qtf)) else return(qtf)
  }
}

#' @rdname EXTLindley
#' @export
rextlindley <- function(n, theta, alpha, beta, L = 1e-4, U = 50)
{
  if(!((alpha == 1) || (alpha <= 0))) stop('alpha must be equal to 1 or less than or equal to 0')
  stopifnot(theta > 0, beta >= 0)
  x  <- qextlindley(p = runif(n), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE, L, U)
  {is <- is.na(x);  na <- any(is)}
  if(na)
  {
    {y <- which(is); i <- 1; l <- length(y)}
    while(i <= l)
    {
      x[y[i]] <- qextlindley(p = runif(1), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE, L, U)
      if(!is.na(x[y[i]])) i <- i + 1
    }
  }
  x
}

#' @rdname EXTLindley
#' @export
hextlindley <- function(x, theta, alpha, beta, log = TRUE)
{
  if(!((alpha == 1) || (alpha <= 0))) stop('alpha must be equal to 1 or less than or equal to 0')
  stopifnot(theta > 0, beta >= 0)
  if(log)
  {
    t1 <- x * theta
    t3 <- log1p(t1 + theta)
    t4 <- log(x)
    t5 <- t1 ^ beta
    t6 <- beta * t5
    t12 <- log(-alpha * theta * x + t1 * t6 + theta * t6 + t6)
    -t3 - t4 + t12
  }
  else
  {
    t1 <- x * theta
    t2 <- t1 ^ beta
    t3 <- beta * t2
    1 / x / (t1 + theta + 1) * (-alpha * theta * x + t1 * t3 + theta * t3 + t3)
  }
}
