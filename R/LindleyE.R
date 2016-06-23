#' @importFrom lamW lambertWm1
#'
#' @name LindleyE
#' @aliases LindleyE dlindleye plindleye qlindleye rlindleye hlindleye
#'
#' @title Lindley Exponential Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the Lindley exponential distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Bhati, D., Malik, M. A., Vaman, H. J., (2015). Lindley-Exponential distribution: properties and applications. \emph{METRON}, \bold{73}, (3), 335–357.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{dlindleye} gives the density, \code{plindleye} gives the distribution function, \code{qlindleye} gives the quantile function, \code{rlindleye} generates random deviates and \code{hlindleye} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[lamW]{lambertWm1}}.
#'
#' @source [d-h-p-q-r]lindleye are calculated directly from the definitions. \code{rlindleye} uses the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta ,\alpha )={\frac{{\theta }^{2}\alpha {{e}^{-\alpha x}}\left(1-{{e}^{-\alpha x}}\right) ^{\theta -1}\left[ 1-\log \left( 1-{{e}^{-\alpha x}}\right) \right] }{1+\theta }}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ,\alpha )={\frac{\left( 1-{{e}^{-\alpha x}}\right) ^{\theta }\left[ 1+\theta -\theta \log \left( 1-{{e}^{-\alpha x}}\right) \right] }{1+\theta }}}
#'
#' Quantile function
#' \deqn{\code{see Bhati et al., 2015}}
#'
#' Hazard rate function
#' \deqn{\code{see Bhati et al., 2015}}
#'
#' @examples
#' set.seed(1)
#' x <- rlindleye(n = 1000, theta = 5.0, alpha = 0.2)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dlindleye(S, theta = 5.0, alpha = 0.2), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' plindleye(q, theta = 5.0, alpha = 0.2, lower.tail = TRUE)
#' plindleye(q, theta = 5.0, alpha = 0.2, lower.tail = FALSE)
#' qlindleye(p, theta = 5.0, alpha = 0.2, lower.tail = TRUE)
#' qlindleye(p, theta = 5.0, alpha = 0.2, lower.tail = FALSE)
#'
#' ## waiting times data (from Ghitany et al., 2008)
#' data(waitingtimes)
#' library(fitdistrplus)
#' fit <- fitdist(waitingtimes, 'lindleye', start = list(theta = 2.6, alpha = 0.15),
#'  lower = c(0.01, 0.01))
#' plot(fit)
#'
#' @rdname LindleyE
#' @export
dlindleye <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t4 <- log1p(theta)
	t5 <- log(alpha)
	t7 <- exp(-alpha * x)
	t8 <- 1 - t7
	t9 <- log(t8)
	t12 <- t8 ^ (theta - 1)
	t15 <- log((1 - t9) * t12 * t7)
	t16 <- 2 * t1 - t4 + t5 + t15
  }
  else
  {
	t1 <- theta ^ 2
	t6 <- exp(-alpha * x)
	t7 <- 1 - t6
	t8 <- log(t7)
	t12 <- t7 ^ (theta - 1)
	t15 <- t1 / (1 + theta) * (1 - t8) * t12 * alpha * t6
  }
}

#' @rdname LindleyE
#' @export
plindleye <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t2 <- exp(q * alpha)
	t3 <- t2 - 1
	t4 <- t3 ^ theta
	t6 <- theta * q * alpha
	t7 <- exp(-t6)
	t9 <- log(t3)
	cdf<- -t4 * t7 * (t9 * theta - t6 - theta - 1) / (1 + theta)
  }
  else
  {
	t1 <- q * alpha
	t2 <- exp(t1)
	t3 <- t2 - 1
	t4 <- t3 ^ theta
	t5 <- t4 * theta
	t8 <- exp(-theta * q * alpha)
	t11 <- log(t3)
	cdf <- (-t5 * t1 * t8 + t5 * t11 * t8 - t4 * t8 - t5 * t8 + theta + 1) / (1 + theta)
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname LindleyE
#' @export
qlindleye <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1  <- 1 + theta
	t3  <- exp(-t1)
	t5  <- lambertWm1(-p * t1 * t3)
	t7  <- 1 / theta
	t9  <- exp(-(t5 + theta + 1) * t7)
	t12 <- log(0.1e1 / (-1 + t9))
	qtf <- (t12 * theta - t5 - theta - 1) / alpha * t7
  }
  else
  {
	t1  <- 1 + theta
	t4  <- exp(-t1)
	t6  <- lambertWm1(t1 * (p - 1) * t4)
	t8  <- 1 / theta
	t10 <- exp(-(t6 + theta + 1) * t8)
	t13 <- log(0.1e1 / (t10 - 1))
	qtf <- -(-t13 * theta + t6 + theta + 1) / alpha * t8
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname LindleyE
#' @export
rlindleye <- function(n, theta, alpha)
{
  stopifnot(theta > 0, alpha > 0)
  qlindleye(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
}

#' @rdname LindleyE
#' @export
hlindleye <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	dlindleye(x, theta, alpha, log = TRUE) - plindleye(x, theta, alpha, lower.tail = FALSE, log.p = TRUE)
  }
  else
  {
	dlindleye(x, theta, alpha, log = FALSE) / plindleye(x, theta, alpha, lower.tail = FALSE, log.p = FALSE)
  }
}

