#' @name TLindley
#' @aliases TLindley dtlindley ptlindley qtlindley rtlindley htlindley
#'
#' @title Transmuted Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the transmuted Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @note The \code{\link[stats]{uniroot}} function with default arguments is used to find out the quantiles.
#'
#' @references
#'
#' Merovci, F., (2013). Transmuted Lindley distribution. \emph{International Journal of Open Problems in Computer Science and Mathematics}, \bold{63}, (3), 63-72.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta positive parameter.
#' @param alpha \eqn{-1 \leq \alpha \leq +1}.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param L,U interval which \code{uniroot} searches for a root (quantile), L = 1e-4 and U = 50 are the default values.
#'
#' @return \code{dtlindley} gives the density, \code{ptlindley} gives the distribution function, \code{qtlindley} gives the quantile function, \code{rtlindley} generates random deviates and \code{htlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[stats]{uniroot}}.
#'
#' @source [d-h-p-q-r]tlindley are calculated directly from the definitions. \code{rtlindley} uses the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta ,\alpha )={\frac{{\theta }^{2}\left( 1+x\right) e{^{-\theta x}}}{1+\theta }\left[ 1-\alpha +2\alpha \left( 1+{\frac{\theta x}{1+\theta }}\right) e{^{-\theta x}}\right] }}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ,\alpha )=\left( 1+\alpha \right) \left[ 1-\left( 1+{\frac{\theta x}{1+\theta }}\right) e{^{-\theta x}}\right] -\alpha \left[1-\left( 1+{\frac{\theta x}{1+\theta }}\right) e{^{-\theta x}}\right]^{2}}
#'
#' Quantile function
#' \deqn{\code{does not have a closed mathematical expression}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha )={\frac{{\theta }^{2}\left( 1+x\right) e{^{-\theta x}\left[ 1-\alpha +2\alpha \left( 1+{\frac{\theta x}{1+\theta }} \right) e{^{-\theta x}}\right] }}{\left( 1+\theta \right) \left\{ \left( 1+\alpha \right) \left[ 1-\left( 1+{\frac{\theta x}{1+\theta }}\right) e{^{-\theta x}}\right] -\alpha \left[ 1-\left( 1+{\frac{\theta x}{1+\theta }}\right) e{^{-\theta x}}\right] ^{2}\right\} }}}
#'
#' \bold{Particular case:} \eqn{\alpha = 0} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rtlindley(n = 1000, theta = 1.5, alpha = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dtlindley(S, theta = 1.5, alpha = 0.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' ptlindley(q, theta = 1.5, alpha = 0.5, lower.tail = TRUE)
#' ptlindley(q, theta = 1.5, alpha = 0.5, lower.tail = FALSE)
#' qtlindley(p, theta = 1.5, alpha = 0.5, lower.tail = TRUE)
#' qtlindley(p, theta = 1.5, alpha = 0.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'tlindley', start = list(theta = 1.5, alpha = 0.5))
#' plot(fit)
#'
#'
#' @rdname TLindley
#' @export
dtlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > -1, alpha < 1)
  if(log)
  {
    t1 <- log(theta)
    t3 <- 1 + theta
    t4 <- log(t3)
    t6 <- log1p(x)
    t7 <- theta * x
    t12 <- exp(-t7)
    t16 <- log1p(-alpha + 2 * alpha * (1 + t7 / t3) * t12)
    2 * t1 - t4 + t6 - t7 + t16
  }
  else
  {
    t1 <- theta ^ 2
    t3 <- 1 / (1 + theta)
    t6 <- theta * x
    t7 <- exp(-t6)
    t1 * t3 * (1 + x) * t7 * (1 - alpha + 2 * alpha * (t6 * t3 + 1) * t7)
  }
}

#' @rdname TLindley
#' @export
ptlindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha >= -1, alpha <= 1)
  if(lower.tail)
  {
    t2  <- theta * q
    t7  <- exp(-t2)
    t9  <- 1 - (1 + t2 / (1 + theta)) * t7
    t11 <- t9 ^ 2
    cdf <- (1 + alpha) * t9 - alpha * t11
  }
  else
  {
    t2  <- theta * q
    t7  <- exp(-t2)
    t9  <- 1 - (1 + t2 / (1 + theta)) * t7
    t11 <- t9 ^ 2
    cdf <- 1 - (1 + alpha) * t9 + alpha * t11
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname TLindley
#' @export
qtlindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > -1, alpha <= 1)
  if(lower.tail)
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) ptlindley(q, theta, alpha, lower.tail = TRUE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
    if(log.p) return(log(qtf)) else return(qtf)
  }
  else
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) ptlindley(q, theta, alpha, lower.tail = FALSE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
    if(log.p) return(log(qtf)) else return(qtf)
  }
}

#' @rdname TLindley
#' @export
rtlindley <- function(n, theta, alpha, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha >= -1, alpha < 1)
  x  <- qtlindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE, L, U)
  is <- is.na(x)
  na <- any(is)
  if(na)
  {
    y <- which(is)
    i <- 1
    l <- length(y)
    while(i <= l)
    {
      x[y[i]] <- qtlindley(p = runif(1), theta, alpha, lower.tail = TRUE, log.p = FALSE, L, U)
      if(!is.na(x[y[i]])) i <- i + 1
    }
  }
  x
}

#' @rdname TLindley
#' @export
htlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha >= -1, alpha <= 1)
  if(log)
  {
	t1 <- log(theta)
	t3 <- 1 + theta
	t4 <- log(t3)
	t6 <- log1p(x)
	t7 <- theta * x
	t10 <- 1 + t7 / t3
	t12 <- exp(-t7)
	t18 <- -t10 * t12 + 1
	t20 <- t18 ^ 2
	t25 <- log((2 * alpha * t10 * t12 - alpha + 1) / (1 - (1 + alpha) * t18 + alpha * t20))
	2 * t1 - t4 + t6 - t7 + t25
  }
  else
  {
	t1 <- theta ^ 2
	t3 <- 1 / (1 + theta)
	t7 <- theta * x
	t8 <- exp(-t7)
	t10 <- t7 * t3 + 1
	t18 <- -t10 * t8 + 1
	t20 <- t18 ^ 2
	t1 * t3 * (1 + x) * t8 * (2 * alpha * t10 * t8 - alpha + 1) / (1 - (1 + alpha) * t18 + alpha * t20)
  }
}
