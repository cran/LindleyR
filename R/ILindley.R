#' @importFrom actuar rinvgamma
#'
#' @importFrom LambertW W
#'
#' @name ILindley
#' @aliases ILindley dilindley hilindley pilindley qilindley rilindley
#'
#' @title Inverse Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the inverse Lindley distribution with parameter theta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Sharma, V. K., Singh, S. K., Singh, U., Agiwal, V., (2015). The inverse Lindley distribution: a stress-strength reliability model with application to head and neck cancer data. \emph{Journal of Industrial and Production Engineering}, \bold{32}, (3), 162-173.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta positive parameter.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical. If TRUE (default), random values are generated from a two-component mixture of inverse-gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dilindley} gives the density, \code{pilindley} gives the distribution function, \code{qilindley} gives the quantile function, \code{rilindley} generates random deviates and \code{hilindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[LambertW]{W}}, \code{\link[actuar]{rinvgamma}}.
#'
#' @source [dpqh]ilindley are calculated directly from the definitions. \code{rilindley} uses either a two-component mixture of inverse gamma distributions or the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta )=\frac{\theta ^{2}}{1+\theta }\left( \frac{1+x}{x^{3}}\right) e^{-\frac{\theta }{x}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta )=\left( 1+\frac{\theta }{x\left( 1+\theta \right) }\right) {e{^{-{\frac{\theta }{x}}}}}}
#'
#' Quantile function
#' \deqn{Q\left( p\mid \theta \right) =-{\frac{\theta }{W\left(-p\left( 1+\theta \right) e{^{-1-\theta }}\right) +1+\theta }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta )=\frac{\theta ^{2}\left( 1+x\right) {e{^{-{\frac{\theta }{x}}}}}}{x^{3}\left( 1+\theta \right) \left[ 1-\left( 1+\frac{\theta }{x\left(1+\theta \right) }\right) {e{^{-{\frac{\theta }{x}}}}}\right] }}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' @examples 
#' set.seed(1)
#' x <- rilindley(n = 1000, theta = 0.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dilindley(S, theta = 1.5), xlab = 'x', ylab = 'pdf', xlim = c(0, 5))
#' hist(x, prob = TRUE, breaks = seq(0, R[2] + 1, 0.5), main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pilindley(q, theta = 1.5, lower.tail = TRUE)
#' pilindley(q, theta = 1.5, lower.tail = FALSE)
#' qilindley(p, theta = 1.5, lower.tail = TRUE)
#' qilindley(p, theta = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'ilindley', start = list(theta = 1.5), lower = c(0))
#' plot(fit)
#'
#'
#' @rdname ILindley
#' @export
dilindley <- function(x, theta, log = FALSE)
{
  stopifnot(theta > 0)
  if(log)
  {
	t1 <- log(theta)
	t4 <- log1p(x)
	t5 <- log(x)
	t10 <- log1p(theta)
	2 * t1 + t4 - 3 * t5 - theta / x - t10
  }
  else
  {
	t1 <- theta ^ 2
	t4 <- x ^ 2
	t9 <- exp(-theta / x)
	t1 * (x + 1) / t4 / x * t9 / (1 + theta)
  }
}

#' @rdname ILindley
#' @export
pilindley <- function(q, theta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0)
  if(lower.tail)
  {
	t3 <- 0.1e1 / q
	t5 <- exp(-t3 * theta)
	O  <- (q * theta + q + theta) * t5 * t3 / (0.1e1 + theta)

  }
  else
  {
	t3 <- 0.1e1 / q
	t5 <- exp(-t3 * theta)
	O  <- 0.1e1 - (q * theta + q + theta) * t5 * t3 / (0.1e1 + theta)
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname ILindley
#' @export
qilindley <- function(p, theta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0)
  if(lower.tail)
  {
	t1 <- 0.1e1 + theta
	t3 <- exp(-t1)
	t5 <- W(-p * t1 * t3,  branch = -1)
	O  <- -theta / (t5 + 1 + theta)
  }
  else
  {
	t1 <- 1 + theta
	t4 <- exp(-t1)
	t6 <- W(t1 * (p - 1) * t4, branch = -1)
	O  <- -theta / (t6 + 1 + theta)
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname ILindley
#' @export
rilindley <- function(n, theta, mixture = TRUE)
{
  stopifnot(theta > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (1 + theta))
    p * rinvgamma(n, shape = 1, scale = theta) + (1 - p) * rinvgamma(n, shape = 2.0, scale = theta)
  }
  else qilindley(p = runif(n), theta, lower.tail = TRUE, log.p = FALSE)
}

#' @rdname ILindley
#' @export
hilindley <- function(x, theta, log = FALSE)
{
  stopifnot(theta > 0)
  if(log)
  {
    t1 <- log(theta)
    t4 <- log1p(x)
    t5 <- log(x)
    t7 <- 1 / x
    t8 <- theta * t7
    t9 <- 1 + theta
    t10 <- log(t9)
    t15 <- exp(-t8)
    t19 <- log(0.1e1 / (1 - (1 + theta / t9 * t7) * t15))
    2 * t1 + t4 - 3 * t5 - t8 - t10 + t19
  }
  else
  {
    t1 <- theta ^ 2
    t4 <- x ^ 2
    t8 <- 1 / x
    t10 <- exp(-theta * t8)
    t12 <- 0.1e1 / (1 + theta)
    t1 * (x + 1) / t4 / x * t10 * t12 / (1 - (theta * t12 * t8 + 1) * t10)
  }
}

