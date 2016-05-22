#' @importFrom VGAM rgengamma.stacy
#'
#' @importFrom LambertW W
#'
#' @name EXTILindley
#' @aliases EXTILindley dextilindley pextilindley qextilindley rextilindley hextilindley
#'
#' @title Extended Inverse Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the extended inverse Lindley distribution with parameters theta, alpha and beta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#' Alkarni, S. H., (2015). Extended inverse Lindley distribution: properties and application. \emph{SpringerPlus}, \bold{4}, (1), 690-703.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha,beta positive parameters.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical. If TRUE (default), random values are generated from a two-component mixture of inverse generalized gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dextilindley} gives the density, \code{pextilindley} gives the distribution function, \code{qextilindley} gives the quantile function, \code{rextilindley} generates random deviates and \code{hextilindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[LambertW]{W}}, \code{\link[VGAM]{gengamma.stacy}}.
#'
#' @source [dpqh]extilindley are calculated directly from the definitions. \code{rextilindley} uses either a two-component mixture of the inverse generalized gamma distributions or the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha,\beta )=\frac{\alpha \theta ^{2}}{\theta +\beta }\left( \frac{\beta +x^{\alpha }}{x^{2\alpha +1}}\right) e^{-\frac{\theta }{x^{\alpha }}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha,\beta )=\left( 1+\frac{\theta \beta }{\left(\theta +\beta \right) {x}^{\alpha }}\right) e{{^{-{\frac{\theta }{x^{\alpha }}}}}}}
#'
#' Quantile function
#' \deqn{Q\left( p\mid \theta,\alpha,\beta \right)=\left[ -\frac{1}{\theta }-\frac{1}{\beta }-\frac{1}{\theta }W_{-1}{\left( -p{\frac{\left( \theta +\beta \right) }{\beta }e{^{-\left( {\frac{\theta +\beta }{\beta }}\right) }}}\right) }\right] ^{-\frac{1}{\alpha }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha,\beta )=\frac{\alpha \theta ^{2}\left( \beta+x^{\alpha }\right) e^{-\frac{\theta }{x^{\alpha }}}}{\left( \theta +\beta\right) x^{2\alpha +1}\left[ 1-\left( 1+\frac{\theta \beta }{\left( \theta+\beta \right) {x}^{\alpha }}\right) e{{^{-{\frac{\theta }{x^{\alpha }}}}}}\right] }}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' @examples
#' set.seed(1)
#' x <- rextilindley(n = 1000, theta = 10, alpha = 10, beta = 10, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dextilindley(S, theta = 10, alpha = 10, beta = 10), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', breaks = seq(0, R[2] + 1, 0.1), add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pextilindley(q, theta = 10, alpha = 10, beta = 10, lower.tail = TRUE)
#' pextilindley(q, theta = 10, alpha = 10, beta = 10, lower.tail = FALSE)
#' qextilindley(p, theta = 10, alpha = 10, beta = 10, lower.tail = TRUE)
#' qextilindley(p, theta = 10, alpha = 10, beta = 10, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'extilindley', start = list(theta = 10, alpha = 10, beta = 10))
#' plot(fit)
#'
#
#' @rdname EXTILindley
#' @export
dextilindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t4 <- log(x)
    t7 <- x ^ alpha
    t14 <- log(0.1e1 / (beta + theta) * (beta + t7))
    t1 + 2 * t2 - 2 * t4 * alpha - t4 - theta / t7 + t14
  }
  else
  {
    t1 <- theta ^ 2
    t6 <- x ^ alpha
    t10 <- x ^ (2 * alpha + 1)
    t15 <- exp(-theta / t6)
    alpha * t1 / (beta + theta) * (beta + t6) / t10 * t15
  }
}

#' @rdname EXTILindley
#' @export
pextilindley <- function(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
    t3 <- q ^ (-alpha)
    t5 <- exp(-theta * t3)
    O  <- 0.1e1 / (beta + theta) * t5 * (t3 * beta * theta + beta + theta)
  }
  else
  {
    t3 <- q ^ (-alpha)
    t5 <- exp(-theta * t3)
    O  <- 1 - 0.1e1 / (beta + theta) * t5 * (t3 * beta * theta + beta + theta)
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname EXTILindley
#' @export
qextilindley <- function(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
    t1 <- 0.1e1 / alpha
    t2 <- theta ^ t1
    t3 <- beta + theta
    t5 <- 0.1e1 / beta
    t7 <- exp(-t3 * t5)
    t10 <- W(-t5 * t7 * p * t3, branch = -1)
    t14 <- ((-beta * t10 - beta - theta) * t5) ^ t1
    O  <-  t2 / t14
  }
  else
  {
    t1 <- 0.1e1 / alpha
    t2 <- theta ^ t1
    t3 <- beta + theta
    t4 <- 0.1e1 / beta
    t6 <- exp(-t3 * t4)
    t11 <- W(t3 * t6 * (p - 1) * t4, branch = -1)
    t15 <- ((-beta * t11 - beta - theta) * t4) ^ t1
    O   <- t2 / t15
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname EXTILindley
#' @export
rextilindley <- function(n, theta, alpha, beta, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (theta + beta))
    1 / (p * rgengamma.stacy(n, scale = 1 / theta ^ (1 / alpha), d = alpha, k = 1) + (1 - p) * rgengamma.stacy(n, scale = 1 / theta ^ (1 / alpha), d = alpha, k = 2))
  }
  else
  {
    qextilindley(p = runif(n), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname EXTILindley
#' @export
hextilindley <- function(x, theta, alpha, beta, log = TRUE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t4 <- beta + theta
    t5 <- log(t4)
    t6 <- x ^ alpha
    t8 <- log(beta + t6)
    t9 <- log(x)
    t15 <- x ^ (-alpha)
    t17 <- exp(-theta * t15)
    t25 <- log(0.1e1 / (1 - 0.1e1 / t4 * t17 * (t15 * beta * theta + beta + theta)))
    t1 + 2 * t2 - t5 + t8 - 2 * t9 * alpha - t9 - theta / t6 + t25
  }
  else
  {
    t1 <- theta ^ 2
    t4 <- 0.1e1 / (beta + theta)
    t6 <- x ^ alpha
    t10 <- x ^ (2 * alpha + 1)
    t15 <- exp(-theta / t6)
    t16 <- x ^ (-alpha)
    t18 <- exp(-theta * t16)
    alpha * t1 * t4 * (beta + t6) / t10 * t15 / (1 - t4 * t18 * (t16 * beta * theta + beta + theta))
  }
}
