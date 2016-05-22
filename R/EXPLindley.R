#' @importFrom LambertW W
#'
#' @name EXPLindley
#' @aliases EXPLindley dexplindley pexplindley qexplindley rexplindley hexplindley
#'
#' @title Exponentiated Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the exponentiated Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#' Nadarajah, S., Bakouch, H. S., Tahmasbi, R., (2011). A generalized Lindley distribution. \emph{Sankhya B}, \bold{73}, (2), 331-359.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{dexplindley} gives the density, \code{pexplindley} gives the distribution function, \code{qexplindley} gives the quantile function, \code{rexplindley} generates random deviates and \code{hexplindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[LambertW]{W}}.
#'
#' @source [dpqh]explindley are calculated directly from the definitions. \code{rexplindley} uses the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha )=\frac{\alpha \theta ^{2}}{(1+\theta )} (1+x)e^{-\theta x}\left[ 1-\left( 1+\frac{\theta x}{1+\theta }\right) e^{-\theta x}\right] ^{\alpha -1}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha )=\left[ 1-\left( 1+\frac{\theta x}{1+\theta }\right) e^{-\theta x}\right] ^{\alpha }}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta,\alpha )=-1-\frac{1}{\theta }-{\frac{1}{\theta }}W_{-1}{\left( \left( 1+\theta \right) \left( p^{-\alpha }-1\right) e{^{-1-\theta }}\right) }}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha )={\frac{\alpha {\theta }^{2}\left( 1+x\right) {{e}^{-\theta x}}\left[ 1-\left( 1+{\frac{\theta ,}{1+\theta }}\right) e{^{-\theta x}}\right] ^{\alpha -1}}{\left( 1+\theta \right) \left\{ 1-\left[\left( 1-\left( 1+{\frac{\theta ,}{1+\theta }}\right) {{e}^{-\theta x}}\right) ^{\alpha }\right] \right\} }}}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rexplindley(n = 1000, theta = 1.5, alpha = 1.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dexplindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pexplindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pexplindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qexplindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qexplindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'explindley', start = list(theta = 1.5, alpha = 1.5))
#' plot(fit)
#'
#' @rdname EXPLindley
#' @export
dexplindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0L, alpha > 0L)
  if(log)
  {
	t1 <- log(alpha)
	t2 <- log(theta)
	t4 <- 1 + theta
	t5 <- log(t4)
	t7 <- log1p(x)
	t8 <- theta * x
	t12 <- exp(-t8)
	t16 <- (1 - t12 * (1 + 0.1e1 / t4 * t8)) ^ (alpha - 1)
	t17 <- log(t16)
	t1 + 2 * t2 - t5 + t7 - t8 + t17
  }
  else
  {
	t1 <- theta ^ 2
	t4 <- 1 / (1 + theta)
	t7 <- theta * x
	t8 <- exp(-t7)
	t15 <- (1 - t8 * (t4 * t7 + 1)) ^ (alpha - 1)
	t15 * t8 * (1 + x) * t4 * alpha * t1
  }
}

#' @rdname EXPLindley
#' @export
pexplindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1 <- theta * q
	t6 <- exp(-t1)
	O  <- (1 - (1 + t1 / (1 + theta)) * t6) ^ alpha
  }
  else
  {
	t1 <- theta * q
	t6 <- exp(-t1)
	t9 <- (1 - t6 * (1 + 1 / (1 + theta) * t1)) ^ alpha
	O  <- 1 - t9
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname EXPLindley
#' @export
qexplindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1  <- 1 + theta
	t2  <- log(p)
	t5  <- exp(0.1e1 / alpha * t2)
	t8  <- exp(-t1)
	t10 <- W(t8 * (t5 - 1) * t1, branch = -1)
	O   <- -1 / theta * (t10 + 1 + theta)
  }
  else
  {
	t1  <- 1 + theta
	t3  <- log1p(-p)
	t6  <- exp(0.1e1 / alpha * t3)
	t9  <- exp(-t1)
	t11 <- W(t9 * (t6 - 1) * t1, branch = -1)
	t13 <- exp(t11 + 1 + theta)
	O   <- -1 / theta / t13 * (theta * t13 + theta * t6 + t13 + t6 - theta - 1)
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname EXPLindley
#' @export
rexplindley <- function(n, theta, alpha)
{
  stopifnot(theta > 0, alpha > 0)
  qexplindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
}

#' @rdname EXPLindley
#' @export
hexplindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(x > 0, theta > 0, alpha > 0)
  if(log)
  {
      t1 <- log(alpha)
      t2 <- log(theta)
      t4 <- 1 + theta
      t5 <- log(t4)
      t7 <- log1p(x)
      t8 <- theta * x
      t12 <- exp(-t8)
      t14 <- 1 - t12 * (1 + 0.1e1 / t4 * t8)
      t16 <- t14 ^ (alpha - 1)
      t17 <- t14 ^ alpha
      t21 <- log(0.1e1 / (1 - t17) * t16)
      t1 + 2 * t2 - t5 + t7 - t8 + t21
  }
  else
  {
      t1 <- theta ^ 2
      t4 <- 1 / (1 + theta)
      t7 <- theta * x
      t8 <- exp(-t7)
      t13 <- 1 - t8 * (t4 * t7 + 1)
      t15 <- t13 ^ (alpha - 1)
      t16 <- t13 ^ alpha
      0.1e1 / (1 - t16) * t15 * t8 * (1 + x) * t4 * t1 * alpha
  }
}
