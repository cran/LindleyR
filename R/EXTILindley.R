#' @importFrom lamW lambertWm1
#' @importFrom stats dgamma
#' @importFrom stats rgamma
#' @importFrom stats qgamma
#'
#' @name EXTILindley
#' @aliases EXTILindley dextilindley pextilindley qextilindley rextilindley hextilindley
#'
#' @title Extended Inverse Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the extended inverse Lindley distribution with parameters theta, alpha and beta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Alkarni, S. H., (2015). Extended inverse Lindley distribution: properties and application. \emph{SpringerPlus}, \bold{4}, (1), 690-703.
#'
#' Mead, M. E., (2015). Generalized inverse gamma distribution and its application in reliability. \emph{Communication in Statistics - Theory and Methods}, \bold{44}, 1426-1435.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha,beta positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a two-component mixture of inverse-gamma distributions, otherwise from the quantile function.
#'#'
#' @return \code{dextilindley} gives the density, \code{pextilindley} gives the distribution function, \code{qextilindley} gives the quantile function, \code{rextilindley} generates random deviates and \code{hextilindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[lamW]{lambertWm1}}.
#'
#' @source [d-h-p-q-r]extilindley are calculated directly from the definitions. \code{rextilindley} uses either a two-component mixture of generalized inverse gamma distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta ,\alpha ,\beta )=\frac{\beta \theta ^{2}}{\theta +\alpha }\left( \frac{\alpha +x^{\beta }}{x^{2\beta +1}}\right) e^{-\frac{\theta }{  x^{\beta }}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ,\alpha ,\beta )=\left( 1+\frac{\theta \alpha }{\left( \theta +\alpha \right) }\frac{1}{{x}^{\beta }}\right) e{{^{-{\frac{\theta }{  x^{\beta }}}}}}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta ,\alpha ,\beta) =\left[ -\frac{1}{\theta }-\frac{1}{\alpha }-\frac{1}{\theta }W_{-1}{\left( -\frac{p}{\alpha }\left( \theta+\alpha \right) {e{^{-\left( {\frac{\theta +\alpha }{\alpha }}\right) }}}\right) }\right] ^{-\frac{1}{\beta }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha ,\beta )=\frac{\beta \theta ^{2}\left( \alpha+x^{\beta }\right) e^{-\frac{\theta }{x^{\beta }}}}{\left( \theta +\alpha\right) x^{2\beta +1}\left[ 1-\left( 1+\frac{\theta \alpha }{\left( \theta+\alpha \right) }\frac{1}{{x}^{\beta }}\right) e{{^{-{\frac{\theta }{x^{\beta }}}}}}\right] }}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular cases:} \eqn{\alpha = 1, \beta = 1} the inverse Lindley distribution, \eqn{\alpha = 1} the generalized inverse Lindley distribution and for \eqn{\alpha = 0} the inverse Weibull distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rextilindley(n = 10000, theta = 5, alpha = 20, beta = 10)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' plot(S, dextilindley(S, theta = 5, alpha = 20, beta = 20), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pextilindley(q, theta = 5, alpha = 20, beta = 10, lower.tail = TRUE)
#' pextilindley(q, theta = 5, alpha = 20, beta = 10, lower.tail = FALSE)
#' qextilindley(p, theta = 5, alpha = 20, beta = 10, lower.tail = TRUE)
#' qextilindley(p, theta = 5, alpha = 20, beta = 10, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'extilindley', start = list(theta = 5, alpha = 20, beta = 10))
#' plot(fit)
#
#' @rdname EXTILindley
#' @export
dextilindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
    t1 <- log(theta)
    t4 <- log(alpha + theta)
    t5 <- x ^ beta
    t10 <- x ^ (2 * beta + 1)
    t14 <- exp(-theta / t5)
    t17 <- log(beta * (alpha + t5) / t10 * t14)
    2 * t1 - t4 + t17
  }
  else
  {
    t1 <- theta ^ 2
    t6 <- x ^ beta
    t10 <- x ^ (2 * beta + 1)
    t15 <- exp(-theta / t6)
    beta * t1 / (alpha + theta) * (alpha + t6) / t10 * t15
  }
}

#' @rdname EXTILindley
#' @export
pextilindley <- function(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
    t4  <- q ^ beta
    t5  <- 0.1e1 / t4
    t10 <- exp(-theta * t5)
    cdf <- (0.1e1 + theta * alpha / (alpha + theta) * t5) * t10
  }
  else
  {
    t4  <- q ^ beta
    t5  <- 0.1e1 / t4
    t10 <- exp(-theta * t5)
    cdf <- 0.1e1 - (0.1e1 + theta * alpha / (alpha + theta) * t5) * t10
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname EXTILindley
#' @export
qextilindley <- function(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
    t1  <- 0.1e1 / beta
    t2  <- theta ^ t1
    t3  <- alpha ^ t1
    t5  <- alpha + theta
    t7  <- 0.1e1 / alpha
    t9  <- exp(-t5 * t7)
    t12 <- lambertWm1(-p * t5 * t7 * t9)
    t15 <- (-alpha * t12 - alpha - theta) ^ (-t1)
    qtf <- t2 * t3 * t15
  }
  else
  {
    t1 <- 0.1e1 / beta
    t2 <- theta ^ t1
    t3 <- alpha ^ t1
    t6 <- alpha + theta
    t8 <- 0.1e1 / alpha
    t10 <- exp(-t6 * t8)
    t13 <- lambertWm1((p - 1) * t6 * t8 * t10)
    t16 <- (-alpha * t13 - alpha - theta) ^ t1
    qtf <- t2 * t3 / t16
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname EXTILindley
#' @export
rextilindley <- function(n, theta, alpha, beta, mixture = TRUE)
{
    if(mixture)
    {
      x <- rbinom(n, size = 1, prob = theta / (theta + alpha))
      x * rigen.gamma	(n, scale = theta, shape1 = 1, shape2 = beta) + (1 - x) * rigen.gamma(n, scale = theta, shape1 = 2, shape2 = beta)
    }
    else qextilindley(p = runif(n), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
}

#' @rdname EXTILindley
#' @export
hextilindley <- function(x, theta, alpha, beta, log = TRUE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
    t1 <- log(beta)
    t2 <- log(theta)
    t4 <- alpha + theta
    t5 <- log(t4)
    t6 <- x ^ beta
    t8 <- log(alpha + t6)
    t9 <- log(x)
    t15 <- x ^ (-beta)
    t17 <- exp(-theta * t15)
    t25 <- log(0.1e1 / (1 - 0.1e1 / t4 * t17 * (t15 * alpha * theta + alpha + theta)))
    t1 + 2 * t2 - t5 + t8 - 2 * t9 * beta - t9 - theta / t6 + t25
  }
  else
  {
    t1 <- theta ^ 2
    t4 <- 0.1e1 / (alpha + theta)
    t6 <- x ^ beta
    t10 <- x ^ (2 * beta + 1)
    t15 <- exp(-theta / t6)
    t16 <- x ^ (-beta)
    t18 <- exp(-theta * t16)
    beta * t1 * t4 * (alpha + t6) / t10 * t15 / (1 - t4 * t18 * (t16 * alpha * theta + alpha + theta))
  }
}
#################################################
dgen.gamma	<- function(x, scale, shape1, shape2)
{
  dgamma(scale * x ^ shape2, shape = shape1, scale = 1) * scale * shape2 * x ^ (shape2 - 1)
}

pgen.gamma	<- function(q, scale, shape1, shape2)
{
  pgamma(scale * q ^ shape2, shape = shape1, scale = 1)
}

qgen.gamma	<- function(p, scale, shape1, shape2)
{
  (qgamma(p, shape = shape1, scale = 1) / scale) ^ (1 / shape2)
}

rgen.gamma	<- function(n, scale, shape1, shape2)
{
  (rgamma(n, shape = shape1, scale = 1) / scale) ^ (1 / shape2)
}
#################################################
digen.gamma	<- function(x, scale, shape1, shape2)
{
  dgamma(scale * x ^ (-shape2), shape = shape1, scale = 1) * scale * shape2 * x ^ (-shape2 - 1)
}

pigen.gamma	<- function(q, scale, shape1, shape2)
{
  pgamma(scale * q ^ (-shape2), shape = shape1, scale = 1, lower.tail = FALSE)
}

qigen.gamma	<- function(p, scale, shape1, shape2)
{
  (qgamma(p, shape = shape1, scale = 1) / scale) ^ (-1 / shape2)
}

rigen.gamma	<- function(n, scale, shape1, shape2)
{
  (rgamma(n, shape = shape1, scale = 1) / scale) ^ (-1 / shape2)
}
