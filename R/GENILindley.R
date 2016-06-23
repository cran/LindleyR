#' @importFrom lamW lambertWm1
#'
#' @name GENILindley
#' @aliases GENILindley dgenilindley pgenilindley qgenilindley rgenilindley hgenilindley
#'
#' @title Generalized Inverse Lindley Distribution
#'
#' @note Barco et al. (2016) named the generalized inverse Lindley distribution as inverse power Lindley distribution.
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the generalized inverse Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Barco, K. V. P., Mazucheli, J. and Janeiro, V. (2016). The inverse power Lindley distribution. \emph{Communications in Statistics - Simulation and Computation}, (to appear).
#'
#' Sharma, V. K., Singh, S. K., Singh, U., Merovci, F., (2015). The generalized inverse Lindley distribution: A new inverse statistical model for the study of upside-down bathtub data. \emph{Communication in Statistics - Theory and Methods}, \bold{0}, 0, 0-0.
#'
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a two-component mixture of generalized inverse gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dgenilindley} gives the density, \code{pgenilindley} gives the distribution function, \code{qgenilindley} gives the quantile function, \code{rgenilindley} generates random deviates and \code{hgenilindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[lamW]{lambertWm1}}.
#'
#' @source [d-h-p-q-r]genilindley are calculated directly from the definitions. \code{rgenilindley} uses either a two-component mixture of generalized inverse gamma distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha )=\frac{\alpha \theta ^{2}}{1+\theta }\left( \frac{1+x^{\alpha }}{x^{2\alpha +1}}\right) e^{-\frac{\theta }{x^{\alpha }}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha )=\left( 1+\frac{\theta }{\left( 1+\theta \right) {x}^{\alpha }}\right) e{{^{-{\frac{\theta }{x^{\alpha }}}}}}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta ,\alpha) =\left( -1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left( -p\left( 1+\theta \right) e{^{-(1+\theta) }}\right) \right) ^{-  \frac{1}{\alpha }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha )={\frac{\alpha {\theta }^{2}\left( 1+{x}^{\alpha}\right) e{{^{-{\frac{\theta }{{x}^{\alpha }}}}}}}{\left( 1+\theta \right) {x}^{2\alpha +1}\left[ 1-\left( 1+\frac{\theta }{\left( 1+\theta \right) {x}^{\alpha }}\right) e{{^{-{\frac{\theta }{x^{\alpha }}}}}}\right] }}}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the inverse Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rgenilindley(n = 1000, theta = 10, alpha = 20, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' plot(S, dgenilindley(S, theta = 10, alpha = 20), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pgenilindley(q, theta = 10, alpha = 20, lower.tail = TRUE)
#' pgenilindley(q, theta = 10, alpha = 20, lower.tail = FALSE)
#' qgenilindley(p, theta = 10, alpha = 20, lower.tail = TRUE)
#' qgenilindley(p, theta = 10, alpha = 20, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'genilindley', start = list(theta = 10, alpha = 20))
#' plot(fit)
#'
#'
#' @rdname GENILindley
#' @export
dgenilindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t5 <- log1p(theta)
    t6 <- x ^ alpha
    t8 <- log1p(t6)
    t9 <- log(x)
    t1 + 2 * t2 - t5 + t8 - 2 * t9 * alpha - t9 - theta / t6
  }
  else
  {
    t1 <- theta ^ 2
    t6 <- x ^ alpha
    t10 <- x ^ (2 * alpha + 1)
    t15 <- exp(-theta / t6)
    alpha * t1 / (1 + theta) * (1 + t6) / t10 * t15
  }
}

#' @rdname GENILindley
#' @export
pgenilindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
    t3  <- q ^ (-alpha)
    t4  <- theta * t3
    t5  <- exp(-t4)
    cdf <- 1 / (1 + theta) * t5 * (theta + 1 + t4)
  }  else
  {
    t3  <- q ^ (-alpha)
    t4  <- theta * t3
    t5  <- exp(-t4)
    cdf <- 1 - 1 / (1 + theta) * t5 * (theta + 1 + t4)
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname GENILindley
#' @export
qgenilindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
    t1  <- 0.1e1 / alpha
    t2  <- theta ^ t1
    t3  <- 1 + theta
    t5  <- exp(-t3)
    t7  <- lambertWm1(-p * t3 * t5)
    t9  <- (-t7 - 1 - theta) ^ t1
    qtf <- t2 / t9
  }
  else
  {
    t1  <- 0.1e1 / alpha
    t2  <- theta ^ t1
    t3  <- 1 + theta
    t6  <- exp(-t3)
    t8  <- lambertWm1(t3 * (p - 1) * t6)
    t10 <- (-t8 - 1 - theta) ^ t1
    qtf <- t2 / t10
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname GENILindley
#' @export
rgenilindley <- function(n, theta, alpha, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0)
  if(mixture)
  {
    x <- rbinom(n, size = 1, prob = theta / (1 + theta))
    ((x * rgamma(n, shape = 1, scale = 1) + (1 - x) * rgamma(n, shape = 2, scale = 1)) / theta) ^ (-1 / alpha)
  }
  else
  {
    qgenilindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname GENILindley
#' @export
hgenilindley <- function(x, theta, alpha, log = TRUE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t4 <- 1 + theta
    t5 <- log(t4)
    t6 <- x ^ alpha
    t8 <- log1p(t6)
    t9 <- log(x)
    t15 <- x ^ (-alpha)
    t16 <- theta * t15
    t17 <- exp(-t16)
    t23 <- log(0.1e1 / (1 - 0.1e1 / t4 * t17 * (theta + 1 + t16)))
    t1 + 2 * t2 - t5 + t8 - 2 * t9 * alpha - t9 - theta / t6 + t23
  }
  else
  {
    t1 <- theta ^ 2
    t4 <- 1 / (1 + theta)
    t6 <- x ^ alpha
    t10 <- x ^ (2 * alpha + 1)
    t15 <- exp(-theta / t6)
    t16 <- x ^ (-alpha)
    t17 <- theta * t16
    t18 <- exp(-t17)
    alpha * t1 * t4 * (1 + t6) / t10 * t15 / (1 - t4 * t18 * (theta + 1 + t17))
  }
}
