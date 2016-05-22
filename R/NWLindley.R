#' @importFrom LambertW W
#'
#' @name NWLindley
#' @aliases NWLindley dnwlindley pnwlindley qnwlindley rnwlindley hnwlindley
#'
#' @title New Weighted Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the new weighted Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @note The \code{\link[stats]{uniroot}} function with default arguments is used to find out the quantiles.
#'
#' @references
#' Asgharzadeh, A., Bakouch, H. S., Nadarajah, S., Sharafi, F., (2016). A new weighted Lindley distribution with application. \emph{Brazilian Journal of Probability and Statistics}, \bold{30}, 1-27.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param L,U interval which \code{uniroot} searches for a root (quantile), L = 1e-4 and U = 50 are the default values.
#'
#' @return \code{dnwlindley} gives the density, \code{pnwlindley} gives the distribution function, \code{qnwlindley} gives the quantile function, \code{rnwlindley} generates random deviates and \code{hnwlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[LambertW]{W}}, \code{\link[stats]{uniroot}}.
#'
#' @source [dpqh]nwlindley are calculated directly from the definitions. \code{rnwlindley} uses the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha )={\frac{{\theta }^{2}\left( 1+\alpha \right) ^{2}}{\alpha \left( \alpha \theta +\alpha +\theta +2\right) }}\left( 1+x\right) \left( 1-e{^{-\theta \alpha x}}\right) e{^{-\theta x}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha )=1-{\frac{\left( 1+\alpha \right) ^{2}\left( \theta x+\theta +1\right) e{^{-\theta x}}}{\alpha \left( \alpha \theta +\alpha +\theta +2\right) }}+{\frac{\left( \theta \alpha x+\alpha \theta +\theta x+\theta +1\right) e{^{-\theta x}}e{^{-\theta \alpha x}}}{\alpha \left(\alpha \theta +\alpha +\theta +2\right) }}}%
#'
#' Quantile function
#' \deqn{\code{does not have a closed mathematical expression}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha )=\frac{{\theta }^{2}\left( 1+\alpha \right) ^{2}\left( 1+x\right) \left( 1-e{^{-\theta \alpha x}}\right) e{^{-\theta x}}}{\left( 1+\alpha \right) ^{2}\left( \theta x+\theta +1\right) e{^{-\theta x}-}\left( \theta \alpha x+\alpha \theta +\theta x+\theta +1\right) e{^{-\theta x}}e{^{-\theta \alpha x}}}}
#'
#' @examples 
#' set.seed(1)
#' x <- rnwlindley(n = 1000, theta = 1.5, alpha = 1.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dnwlindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pnwlindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pnwlindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qnwlindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qnwlindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'nwlindley', start = list(theta = 1.5, alpha = 1.5))
#' plot(fit)
#'
#
#' @rdname NWLindley
#' @export
dnwlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t3 <- 1 + alpha
	t4 <- log(t3)
	t6 <- alpha * theta
	t11 <- log(t6 * t3 + alpha * (2 + alpha))
	t13 <- log1p(x)
	t15 <- exp(-t6 * x)
	t17 <- log1p(-t15)
	-theta * x + 2 * t1 - t11 + t13 + t17 + 2 * t4
  }
  else
  {
	t1 <- theta ^ 2
	t2 <- 1 + alpha
	t3 <- t2 ^ 2
	t5 <- alpha * theta
	t14 <- exp(-t5 * x)
	t18 <- exp(-theta * x)
	t1 * t3 / (t5 * t2 + alpha * (2 + alpha)) * (1 + x) * (1 - t14) * t18
  }
}

#' @rdname NWLindley
#' @export
pnwlindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t2 <- alpha * q * theta
	t3 <- exp(-t2)
	t5 <- theta * q
	t7 <- t3 * theta
	t11 <- exp(t5)
	t12 <- alpha ^ 2
	t13 <- t11 * t12
	t17 <- t11 * alpha
	t22 <- alpha * theta
	t25 <- t3 * alpha * t5 - t12 * q * theta + t3 * q * theta + t7 * alpha - t12 * theta + t13 * theta + t17 * theta - 2 * alpha - t12 + t13 + 2 * t17 - 2 * t2 - 2 * t22 + t3 - t5 + t7 - theta - 1
	t26 <- exp(-t5)
	O   <- t25 * t26 / alpha / (t22 + alpha + theta + 2)
  }
  else
  {
	t2 <- alpha * q * theta
	t3 <- exp(-t2)
	t5 <- theta * q
	t7 <- t3 * theta
	t11 <- alpha ^ 2
	t16 <- alpha * theta
	t19 <- t3 * alpha * t5 - t11 * q * theta + t3 * q * theta + t7 * alpha - t11 * theta - 2 * alpha - t11 - 2 * t16 - 2 * t2 + t3 - t5 + t7 - theta - 1
	t20 <- exp(-t5)
	O 	<- -t19 * t20 / alpha / (t16 + alpha + theta + 2)
  }
  if(log.p) return(log(O)) else return(O)
}


#' @rdname NWLindley
#' @export
qnwlindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pnwlindley(q, theta, alpha, lower.tail = TRUE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    O <- sapply(p, fx)
    if(log.p) return(log(O)) else return(O)
  }
  else
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pnwlindley(q, theta, alpha, lower.tail = FALSE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    O <- sapply(p, fx)
    if(log.p) return(log(O)) else return(O)
  }
}

#' @rdname NWLindley
#' @export
rnwlindley <- function(n, theta, alpha, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > 0)
  x  <- qnwlindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE, L, U)
  is <- is.na(x)
  na <- any(is)
  if(na)
  {
    y <- which(is)
    i <- 1
    l <- length(y)
    while(i <= l)
    {
      x[y[i]] <- qnwlindley(p = runif(1), theta, alpha, lower.tail = TRUE, log.p = FALSE, L, U)
      if(!is.na(x[y[i]])) i <- i + 1
    }
  }
  x
}

#' @rdname NWLindley
#' @export
hnwlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t3 <- 1 + alpha
	t4 <- log(t3)
	t6 <- alpha * theta
	t11 <- log(t6 * t3 + alpha * (2 + alpha))
	t13 <- log1p(x)
	t14 <- t6 * x
	t15 <- exp(-t14)
	t17 <- log1p(-t15)
	t18 <- theta * x
	t20 <- exp(-t18)
	t22 <- 1 / alpha
	t24 <- 0.1e1 / (t6 + alpha + theta + 2)
	t28 <- t3 ^ 2
	t35 <- log(-(t14 + t6 + t18 + theta + 1) * t20 * t22 * t24 * t15 + t28 * (t18 + theta + 1) * t20 * t22 * t24)
	2 * t1 + 2 * t4 - t11 + t13 + t17 - t18 - t35
  }
  else
  {
	t1 <- theta ^ 2
	t2 <- 1 + alpha
	t3 <- t2 ^ 2
	t5 <- alpha * theta
	t13 <- t5 * x
	t14 <- exp(-t13)
	t17 <- theta * x
	t18 <- exp(-t17)
	t21 <- 1 / alpha
	t23 <- 1 / (t5 + alpha + theta + 2)
	t1 * t3 / (t5 * t2 + alpha * (2 + alpha)) * (1 + x) * (1 - t14) * t18 / (-(t13 + t5 + t17 + theta + 1) * t18 * t21 * t23 * t14 + t3 * (t17 + theta + 1) * t18 * t21 * t23)
  }
}

