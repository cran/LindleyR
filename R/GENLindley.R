#' @name GENLindley
#' @aliases GENLindley dgenlindley pgenlindley qgenlindley rgenlindley hgenlindley
#'
#' @title Generalized Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the generalized Lindley distribution with parameters theta, alpha and beta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @note The \code{\link[stats]{uniroot}} function with default arguments is used to find out the quantiles.
#'
#' @references
#'
#' Zakerzadeh, H., Dolati, A., (2009). Generalized Lindley distribution. \emph{Journal of Mathematical Extension}, \bold{3}, (2), 13–25.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha,beta positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param L,U interval which \code{uniroot} searches for a root (quantile), L = 1e-4 and U = 50 are the default values.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dgenlindley} gives the density, \code{pgenlindley} gives the distribution function, \code{qgenlindley} gives the quantile function, \code{rgenlindley} generates random deviates and \code{hgenlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[lamW]{lambertWm1}}, \code{\link[stats]{uniroot}}.
#'
#' @source [d-h-p-q-r]genlindley are calculated directly from the definitions. \code{rgenlindley} uses either a two-component mixture of the gamma distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f( x\mid \theta,\alpha,\beta) =\frac{\theta ^{\alpha +1}}{\left( \theta +\beta \right) \Gamma \left( \alpha +1\right) }x^{\alpha	-1}\left( \alpha +\beta x\right) e^{-\theta x}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha,\beta )=\sum_{j=0}^{1}\left\vert j-\frac{\theta }{\left( \theta +\beta \right) }\right\vert \frac{\Gamma \left( \alpha -j,\theta x\right) }{\Gamma \left( \alpha -j\right) }}%
#'
#' Quantile function
#' \deqn{\code{does not have a closed mathematical expression}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha,\beta )=\frac{\theta ^{\alpha +1}x^{\alpha -1}\left(\alpha +\beta x\right) e^{-\theta x}}{\left( \theta +\beta \right) \Gamma   \left( \alpha +1\right) \sum\limits_{j=0}^{1}\left\vert j-\frac{\theta }{    \left( \theta +\beta \right) }\right\vert \frac{\Gamma \left( \alpha -j,\theta x\right) }{\Gamma \left( \alpha -j\right) }}}
#'
#' where \eqn{\Gamma \left( a,b\right)} is the lower incomplete gamma function.
#'
#' \bold{Particular cases:} \eqn{(\alpha=1, \beta = 1)} the one-parameter Lindley distribution, \eqn{\alpha=1} the two-parameter Lindley distribution, \eqn{(\alpha=1,\beta=0)} the exponential distribution, \eqn{\beta = 0} the gamma distribution and for \eqn{\beta=\alpha} the weighted Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rgenlindley(n = 1000, theta = 1.5, alpha = 1.5, beta = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dgenlindley(S, theta = 1.5, alpha = 1.5, beta = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pgenlindley(q, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = TRUE)
#' pgenlindley(q, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = FALSE)
#' qgenlindley(p, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = TRUE)
#' qgenlindley(p, theta = 1.5, alpha = 1.5, beta = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'genlindley', start = list(theta = 1.5, alpha = 1.5, beta = 1.5))
#' plot(fit)
#
#'
#' @rdname GENLindley
#' @export
dgenlindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
	t2 <- log(theta)
	t5 <- log(x)
	t8 <- log(beta + theta)
	t11 <- log(beta * x + alpha)
	t13 <- gamma(alpha)
	t14 <- log(t13)
	t15 <- log(alpha)
	(alpha + 1) * t2 + (alpha - 1) * t5 - t8 + t11 - theta * x - t14 - t15
  }
  else
  {
	t1 <- alpha + 1
	t2 <- theta ^ t1
	t4 <- x ^ (alpha - 1)
	t12 <- exp(-theta * x)
	t14 <- gamma(t1)
	t2 * t4 / (beta + theta) * (beta * x + alpha) * t12 / t14
  }
}

#' @rdname GENLindley
#' @export
pgenlindley <- function(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
	cdf <- theta / (theta + beta) * pgamma(q, shape = alpha, rate = theta, lower.tail = TRUE) + beta / (beta + theta) * pgamma(q, shape = alpha + 1, rate = theta, lower.tail = TRUE)
  }
  else
  {
	cdf <- theta / (theta + beta) * pgamma(q, shape = alpha, rate = theta, lower.tail = FALSE) + beta / (beta + theta) * pgamma(q, shape = alpha + 1, rate = theta, lower.tail = FALSE)
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname GENLindley
#' @export
qgenlindley <- function(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pgenlindley(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
    if(log.p) return(log(qtf)) else return(qtf)
  }
  else
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pgenlindley(q, theta, alpha, beta, lower.tail = FALSE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
   if(log.p) return(log(qtf)) else return(qtf)
  }
}

#' @rdname GENLindley
#' @export
rgenlindley <- function(n, theta, alpha, beta, mixture = TRUE, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(mixture)
  {
    x <- rbinom(n, size = 1, prob = theta / (theta + beta))
    x <- x * rgamma(n, shape = alpha, rate = theta) + (1 - x) * rgamma(n, shape = alpha + 1, rate = theta)
  }
  else
  {
    x  <- qgenlindley(p = runif(n), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE, L, U)
    {is <- is.na(x); na <- any(is)}
    if(na)
    {
      {y <- which(is); i <- 1; l <- length(y)}
      while(i <= l)
      {
        x[y[i]] <- qgenlindley(p = runif(1), theta, alpha, beta, lower.tail = TRUE, log.p = FALSE, L, U)
        if(!is.na(x[y[i]])) i <- i + 1
      }
    }
  }
  x
}

#' @rdname GENLindley
#' @export
hgenlindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
	dgenlindley(x, theta, alpha, beta, log = TRUE) - pgenlindley(x, theta, alpha, beta, lower.tail = FALSE, log.p = TRUE)
  }
  else
  {
	dgenlindley(x, theta, alpha, beta, log = FALSE) / pgenlindley(x, theta, alpha, beta, lower.tail = FALSE, log.p = FALSE)
  }
}


