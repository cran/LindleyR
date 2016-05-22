#' @importFrom LambertW W
#'
#' @name QLindley
#' @aliases QLindley dqlindley pqlindley qqlindley rqlindley hqlindley
#'
#' @title Quasi Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the quasi Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Shanker, R. and Mishra, A. (2013). A quasi Lindley distribution. \emph{African Journal of Mathematics and Computer Science Research},  \bold{6}, (4), 64-71.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta positive parameter.
#' @param alpha greater than -1.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical. If TRUE (default), random values are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dqlindley} gives the density, \code{pqlindley} gives the distribution function, \code{qqlindley} gives the quantile function, \code{rqlindley} generates random deviates and \code{hqlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[LambertW]{W}}.
#'
#' @source [dpqh]qlindley are calculated directly from the definitions. \code{rqlindley} uses either a two-component mixture of gamma distributions or the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta ,\alpha )={\frac{\theta \left( \alpha +\theta x\right) {{e}^{-\theta x}}}{1+\alpha }}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ,\alpha )=1-{\frac{\left( 1+\alpha +\theta x\right) }{1+\alpha }{e}^{-\theta x}}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta ,\alpha )=-\frac{1}{{\theta }}-{\frac{\alpha }{\theta }}-\frac{1}{{\theta }}{W}_{-1}\left( \left( p-1\right) \left( 1+\alpha	\right) {{e}^{-1-\alpha }}\right) }
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha )=\frac{\theta \left( \alpha +\theta x\right) }{\left( 1+\alpha +\theta x\right) }}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular cases:} \eqn{\alpha = \theta} the one-parameter Lindley distribution and for \eqn{\alpha=0} the gamma distribution with shape = 2 and scale = \eqn{\theta}.
#'
#' @examples 
#' set.seed(1)
#' x <- rqlindley(n = 1000, theta = 1.5, alpha = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dqlindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pqlindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pqlindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qqlindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qqlindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'qlindley', start = list(theta = 1.5, alpha = 1.5))
#' plot(fit)
#'
#'
#' @rdname QLindley
#' @export
dqlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t3 <- log1p(alpha)
	t4 <- theta * x
	t6 <- log(t4 + alpha)
	t1 - t3 - t4 + t6
  }
  else
  {
	t1 <- theta * x
	t6 <- exp(-t1)
	theta * (t1 + alpha) / (1 + alpha) * t6
  }
}

#' @rdname QLindley
#' @export
pqlindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t2 <- exp(-theta * q)
	O  <- -(t2 * q * theta + t2 * alpha - alpha + t2 - 1) / (1 + alpha)
  }
  else
  {
	t2 <- exp(-theta * q)
	O  <- 1 + (t2 * q * theta + t2 * alpha - alpha + t2 - 1) / (1 + alpha)
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname QLindley
#' @export
qqlindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1 <- 0.1e1 / theta
	t3 <- 1 + alpha
	t5 <- exp(-t3)
	t7 <- W((p - 1) * t3 * t5, branch = -1)
	O  <- -t1 * alpha - t1 * t7 - t1
  }
  else
  {
	t1 <- 0.1e1 / theta
	t2 <- 1 + alpha
	t4 <- exp(-t2)
	t6 <- W(-p * t2 * t4, branch = -1)
	O  <- -t1 * alpha - t1 * t6 - t1
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname QLindley
#' @export
rqlindley <- function(n, theta, alpha, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = alpha / (1 + alpha))
    p * rgamma(n, shape = 1, rate = theta) + (1 - p) * rgamma(n, shape = 2, rate = theta)
  }
  else
  {
    qqlindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname QLindley
#' @export
hqlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t2 <- theta * x
	t4 <- log(t2 + alpha)
	t5 <- 1 + alpha
	t6 <- log(t5)
	t7 <- exp(-t2)
	t15 <- log1p((t7 * x * theta + t7 * alpha - alpha + t7 - 1) / t5)
	t1 + t4 - t6 - t2 - t15
  }
  else
  {
	t1 <- theta * x
	t5 <- 1 / (1 + alpha)
	t6 <- exp(-t1)
	theta * (t1 + alpha) * t5 * t6 / (1 + (t6 * x * theta + t6 * alpha - alpha + t6 - 1) * t5)
  }
}

