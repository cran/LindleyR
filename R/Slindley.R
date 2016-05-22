#' @importFrom LambertW W
#'
#' @name SLindley
#' @aliases SLindley dslindley pslindley qslindley rslindley hslindley
#'
#' @title Two-Parameter Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random generation and hazard rate function for the two-parameter Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Shanker, R., Sharma, S. and Shanker, R. (2013). A two-parameter Lindley distribution for modeling waiting and survival times data. \emph{Applied Mathematics}, \bold{4}, (2), 363-368.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameter.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical. If TRUE (default), random values are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dslindley} gives the density, \code{pslindley} gives the distribution function, \code{qslindley} gives the quantile function, \code{rslindley} generates random deviates and \code{hslindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[LambertW]{W}}.
#'
#' @source [dpqh]slindley are calculated directly from the definitions. \code{rslindley} uses either a two-component mixture of the gamma distributions or the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha )=\frac{{\theta }^{2}}{\theta +\alpha }\left(1+\alpha x\right) e^{-\theta x}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha )=1-\frac{\left( \alpha +\theta \alpha +\alpha \theta x\right) }{\theta +\alpha }e^{-\theta x}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta,\alpha )=-\frac{1}{\theta }-\frac{1}{\alpha }-\frac{1}{\theta }W_{-1}\left( \frac{1}{\alpha }(p-1)\left( \theta +\alpha \right)e^{-\left( {\frac{\alpha +\theta }{\alpha }}\right) }\right)}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta )=\frac{\theta ^{2}}{\left( \alpha +\theta \alpha +\alpha\theta x\right) }(1+\alpha x)}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the one-parameter Lindley distribution.
#'
#' @examples 
#' set.seed(1)
#' x <- rslindley(n = 1000, theta = 1.5, alpha = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dslindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pslindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pslindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qslindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qslindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'slindley', start = list(theta = 1.5, alpha = 1.5))
#' plot(fit)
#'
#'
#' @rdname SLindley
#' @export
dslindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t4 <- log(theta + alpha)
	t7 <- log1p(alpha * x)
	-theta * x + 2 * t1 - t4 + t7
  }
  else
  {
	t1 <- theta ^ 2
	t8 <- exp(-theta * x)
	t1 / (theta + alpha) * (alpha * x + 1) * t8
  }
}

#' @rdname SLindley
#' @export
pslindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t8 <- exp(-theta * q)
	O  <- 1 - (alpha * theta * q + alpha + theta) / (theta + alpha) * t8
  }
  else
  {
  	t8 <- exp(-theta * q)
	O  <- (alpha * theta * q + alpha + theta) / (theta + alpha) * t8
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname SLindley
#' @export
qslindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1 <- 0.1e1 / theta
	t3 <- theta + alpha
	t5 <- 0.1e1 / alpha
	t7 <- exp(-t3 * t5)
	t10 <- W((p - 1) * t3 * t5 * t7, branch = -1)
	O   <- -t1 * t10 - t1 - t5
  }
  else
  {
	t1 <- 0.1e1 / theta
	t2 <- theta + alpha
	t4 <- 0.1e1 / alpha
	t6 <- exp(-t2 * t4)
	t9 <- W(-p * t2 * t4 * t6, branch = -1)
	O  <- -t1 * t9 - t1 - t4
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname SLindley
#' @export
rslindley <- function(n, theta, alpha, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (theta + alpha))
    p * rgamma(n, shape = 1, rate = theta) + (1 - p) * rgamma(n, shape = 2, rate = theta)
  }
  else
  {
    qslindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname SLindley
#' @export
hslindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t1 <- log(theta)
    t5 <- log1p(alpha * x)
    t9 <- log(alpha * theta * x + alpha + theta)
    2 * t1 + t5 - t9
  }
  else
  {
    t1 <- theta ^ 2
    t1 * (alpha * x + 1) / (alpha * theta * x + alpha + theta)
  }
}

