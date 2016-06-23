#' @importFrom lamW lambertWm1
#'
#' @name EXTPLindley
#' @aliases EXTPLindley dextplindley pextplindley qextplindley rextplindley hextplindley
#'
#' @title Extended Power Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the extended power Lindley distribution with parameters theta, alpha and beta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Alkarni, S. H., (2015). Extended power Lindley distribution: A new statistical model for non-monotone survival data. \emph{European Journal of Statistics and Probability}, \bold{3}, (3), 19-34.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha,beta positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dextplindley} gives the density, \code{pextplindley} gives the distribution function, \code{qextplindley} gives the quantile function, \code{rextplindley} generates random deviates and \code{hextplindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[lamW]{lambertWm1}}.
#'
#' @source [d-h-p-q-r]extplindley are calculated directly from the definitions. \code{rextplindley} uses either a two-component mixture of gamma distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha,\beta )={\frac{\alpha \theta ^{2}}{\theta +\beta }}(1+\beta x^{\alpha })\ x^{\alpha -1}\ e^{-\theta x^{\alpha }}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha,\beta )=1-\left( 1+{\frac{\beta \theta x^{\alpha }}{\theta +\beta }}\right) \ e^{-\theta x^{\alpha }}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta ,\alpha ,\beta )={\left[ -\frac{1}{\theta }-\frac{1}{\beta }-{\frac{1}{\theta }}W_{-1}{\left( \frac{1}{\beta }\left( p-1\right) \left(  \beta +\theta \right) e{{^{-\left( {\frac{\beta +\theta }{\beta }}\right) }}}\right) }\right] }^{\frac{1}{\alpha }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha ,\beta )={\frac{\alpha {\theta }^{2}\left( 1+\beta {x}^{\alpha }\right) {x}^{\alpha -1}}{\left( \beta +\theta \right) {\left(1+{\frac{\beta \theta {x}^{\alpha }}{\beta +\theta }}\right) }}}}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular cases:} \eqn{\beta = 1} the power Lindley distribution, \eqn{\alpha = 1} the two-parameter Lindley distribution and \eqn{(\alpha = 1, \beta = 1)} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rextplindley(n = 1000, theta = 1.5, alpha = 1.5, beta = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dextplindley(S, theta = 1.5, alpha = 1.5, beta = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pextplindley(q, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = TRUE)
#' pextplindley(q, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = FALSE)
#' qextplindley(p, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = TRUE)
#' qextplindley(p, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'extplindley', start = list(theta = 1.5, alpha = 1.5, beta = 1.5))
#' plot(fit)
#'
#' @rdname EXTPLindley
#' @export
dextplindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
	t1 <- log(alpha)
	t2 <- log(theta)
	t5 <- log(x)
	t7 <- x ^ alpha
	t14 <- log(0.1e1 / (beta + theta) * (beta * t7 + 1))
	t1 + 2 * t2 + (alpha - 1) * t5 - theta * t7 + t14
  }
  else
  {
	t1 <- theta ^ 2
	t6 <- x ^ alpha
	t10 <- x ^ (alpha - 1)
	t13 <- exp(-theta * t6)
	alpha * t1 / (beta + theta) * (beta * t6 + 1) * t10 * t13
  }
}

#' @rdname EXTPLindley
#' @export
pextplindley <- function(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
	t4  <- q ^ alpha
	t9  <- exp(-theta * t4)
	cdf <- 1 - (1 + beta * theta / (beta + theta) * t4) * t9
  }
  else
  {
	t4  <- q ^ alpha
	t9  <- exp(-theta * t4)
	cdf <- (1 + beta * theta / (beta + theta) * t4) * t9
 }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname EXTPLindley
#' @export
qextplindley <- function(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
	t1 <- 0.1e1 / alpha
	t2 <- theta ^ t1
	t5 <- beta + theta
	t7 <- 0.1e1 / beta
	t9 <- exp(-t5 * t7)
	t12 <- lambertWm1((p - 1) * t5 * t7 * t9)
	t16 <- ((-beta * t12 - beta - theta) * t7) ^ t1
	qtf <- 0.1e1 / t2 * t16
  }
  else
  {
	t1 <- 0.1e1 / alpha
	t2 <- theta ^ t1
	t4 <- beta + theta
	t6 <- 0.1e1 / beta
	t8 <- exp(-t4 * t6)
	t11 <- lambertWm1(-p * t4 * t6 * t8)
	t15 <- ((-beta * t11 - beta - theta) * t6) ^ t1
	qtf <- 0.1e1 / t2 * t15
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname EXTPLindley
#' @export
rextplindley <- function(n, theta, alpha, beta, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (beta + theta))
    (p * rgamma(n, shape = 1, rate = theta) + (1 - p) * rgamma(n, shape = 2, rate = theta)) ^ (1 / alpha)
  }
  else
  {
    qextplindley(p = runif(n), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname EXTPLindley
#' @export
hextplindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t4 <- beta + theta
    t5 <- log(t4)
    t6 <- x ^ alpha
    t9 <- log1p(beta * t6)
    t11 <- log(x)
    t18 <- log1p(beta * theta / t4 * t6)
    t1 + 2 * t2 - t5 + t9 + (alpha - 1) * t11 - t18
  }
  else
  {
    t1 <- theta ^ 2
    t4 <- 0.1e1 / (beta + theta)
    t6 <- x ^ alpha
    t10 <- x ^ (alpha - 1)
    alpha * t1 * t4 * (beta * t6 + 1) * t10 / (beta * theta * t4 * t6 + 1)
  }
}
