#' @importFrom lamW lambertWm1
#'
#' @name GAMLindley
#' @aliases GAMLindley dgamlindley pgamlindley qgamlindley rgamlindley hgamlindley
#'
#' @title Gamma Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the Gamma Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Nedjar, S. and Zeghdoudi (2016). On gamma Lindley distribution: Properties and simulations. \emph{Journal of Computational and Applied Mathematics}, \bold{298}, 167-174.
#'
#' Zeghdoudi, H, and Nedjar, S. (2015) Gamma Lindley distribution and its application. \emph{Journal of Applied Probability and Statistics}, \bold{11}, (1), 1-11.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a mixture of gamma and one-parameter Lindley distributions, otherwise from the quantile function.
#'
#' @return \code{dgamlindley} gives the density, \code{pgamlindley} gives the distribution function, \code{qgamlindley} gives the quantile function, \code{rgamlindley} generates random deviates and \code{hgamlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[lamW]{lambertWm1}}.
#'
#' @source [d-h-p-q-r]gamlindley are calculated directly from the definitions. \code{rgamlindley} uses either a mixture of gamma and one-parameter Lindley distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta ,\alpha )=\frac{{\theta }^{2}}{\alpha \left( 1+\theta \right) }\left[ \left( \alpha +\alpha \theta -\theta \right) x+1\right] e{^{-\theta x}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ,\alpha )=\frac{1}{\alpha \left( 1+\theta \right) }\left[\left( \alpha +\alpha \theta -\theta \right) \left( 1+\theta x\right) +\theta \right] e{^{-\theta  x}}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta ,\alpha )=-\frac{\alpha \left( 1+\theta \right) }{\theta \left[ \left( \alpha +\alpha \theta -\theta \right) \right] }-\frac{1}{\theta }W_{-1}\left( {\frac{\left( 1+\theta \right) \alpha \left(p-1\right) }{\alpha +\alpha \theta -\theta }}e{{^{-{\frac{\left( 1+\theta \right) \alpha }{\alpha \theta +\alpha -\theta }}}}}\right) }
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha )=\frac{{\theta }^{2}\left[ \left( \alpha +\alpha \theta -\theta \right) x+1\right] }{{\theta }\left( \alpha +\alpha \theta-\theta \right) x+\alpha \left( 1+\theta \right) }}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rgamlindley(n = 1000, theta = 1.5, alpha = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dgamlindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pgamlindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pgamlindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qgamlindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qgamlindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'gamlindley', start = list(theta = 1.5, alpha = 1.5))
#' plot(fit)
#'
#' @rdname GAMLindley
#' @export
dgamlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1  <- log(theta)
	t4  <- log(alpha)
	t6  <- log1p(theta)
	t11 <- log1p((alpha * theta + alpha - theta) * x)
	-theta * x + 2 * t1 + t11 - t4 - t6
  }
  else
  {
	t1 <- theta ^ 2
	t8 <- exp(-theta * x)
	t1 * ((alpha * theta + alpha - theta) * x + 1) * t8 / (alpha * (1 + theta))
  }
}

#' @rdname GAMLindley
#' @export
pgamlindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1  <- theta * q
	t2  <- exp(-t1)
	t3  <- t2 * alpha
	t4  <- theta ^ 2
	cdf <- -(-t2 * q * t4 + t3 * q * t4 - alpha * theta + t3 * t1 + t3 * theta - alpha + t3) / alpha / (1 + theta)
  }
  else
  {
	t2  <- exp(-theta * q)
	t3  <- alpha * q
	t4  <- theta ^ 2
	cdf <- t2 * (alpha * theta - q * t4 + t3 * t4 + t3 * theta + alpha) / (1 + theta) / alpha
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname GAMLindley
#' @export
qgamlindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t2  <- (1 + theta) * alpha
	t3  <- 1 / theta
	t6  <- 1 / (alpha * theta + alpha - theta)
	t12 <- exp(-t2 * t6)
	t15 <- lambertWm1(t2 * (p - 1) * t6 * t12)
	qtf <- -t2 * t3 * t6 - t15 * t3
  }
  else
  {
	t2  <- (1 + theta) * alpha
	t3  <- 1 / theta
	t6  <- 1 / (alpha * theta + alpha - theta)
	t12 <- exp(-t2 * t6)
	t15 <- lambertWm1(-t2 * p  * t6 * t12)
	qtf <- -t2 * t3 * t6 - t15 * t3
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname GAMLindley
#' @export
rgamlindley <- function(n, theta, alpha, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0)
  if(mixture)
  {
	x <- rbinom(n, size = 1, prob = (alpha - 1) / alpha)
	x * rgamma(n, shape = 2, rate = theta) + (1 - x) * rlindley(n, theta, mixture = T)
  }
  else
  {
    qgamlindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname GAMLindley
#' @export
hgamlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
  	dgamlindley(x, theta, alpha, log = TRUE) - pgamlindley(x, theta, alpha, lower.tail = FALSE, log.p = TRUE)
  }
  else
  {
  	dgamlindley(x, theta, alpha, log = FALSE) / pgamlindley(x, theta, alpha, lower.tail = FALSE, log.p = FALSE)
  }
}
