#' @importFrom LambertW W
#'
#' @name EXPPLindley
#' @aliases EXPPLindley dexpplindley pexpplindley qexpplindley rexpplindley hexpplindley
#'
#' @title Exponentiated Power Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the exponentiated power Lindley distribution with parameters theta, alpha and beta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Ashour, S. K., Eltehiwy, M. A., (2015). Exponentiated power Lindley distribution. \emph{Journal of Advanced Research}, \bold{6}, (6), 895-905.
#'
#' Warahena-Liyanage, G., Pararai, M., (2014). A generalized power Lindley distribution with applications. \emph{Asian Journal of Mathematics and Applications}, \bold{2014}, 1-23.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha,beta positive parameters.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical. If TRUE (default), random values are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dexpplindley} gives the density, \code{pexpplindley} gives the distribution function, \code{qexpplindley} gives the quantile function, \code{rexpplindley} generates random deviates and \code{hexpplindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[LambertW]{W}}.
#'
#' @source [dpqh]expplindley are calculated directly from the definitions. \code{rexpplindley} uses either a two-component mixture of gamma distributions or the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha,\beta )={\frac{\beta \alpha \theta ^{2}}{1 + \theta}}(1+x^{\alpha })x^{\alpha -1}e^{-\theta x^{\alpha }}\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) e^{-\theta x^{\alpha }}\right]^{\beta -1}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha,\beta )=\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) \ e^{-\theta x^{\alpha }}\right] ^{\beta }}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta,\alpha,\beta )=\left( -1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left( \left( 1+\theta \right) \left( p^{-\beta }-1\right)e^{-1-\theta }\right) \right) ^{\frac{1}{\alpha }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid\theta,\alpha,\beta )={\frac{\beta \alpha \theta ^{2}(1+x^{\alpha })x^{\alpha -1}e^{-\theta x^{\alpha }}\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) e^{-\theta x^{\alpha }}\right] ^{\beta -1}}{\left( \theta +1\right) \left\{ 1-\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) \ e^{-\theta x^{\alpha }}\right] ^{\beta }\right\} }}}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular cases:} \eqn{\alpha = 1} the exponentiated Lindley distribution, \eqn{\beta = 1} the power Lindley distribution and \eqn{(\alpha = 1, \beta = 1)} the one-parameter Lindley distribution. See Warahena-Liyanage and Pararai (2014) for other particular cases.
#'
#' @examples
#' set.seed(1)
#' x <- rexpplindley(n = 1000, theta = 1.5, alpha = 2.5, beta = 1.1, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dexpplindley(S, theta = 1.5, alpha = 2.5, beta = 1.1), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pexpplindley(q, theta = 1.5, alpha = 2.5, beta = 1.1, lower.tail = TRUE)
#' pexpplindley(q, theta = 1.5, alpha = 2.5, beta = 1.1, lower.tail = FALSE)
#' qexpplindley(p, theta = 1.5, alpha = 2.5, beta = 1.1, lower.tail = TRUE)
#' qexpplindley(p, theta = 1.5, alpha = 2.5, beta = 1.1, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'expplindley', start = list(theta = 1.5, alpha = 2.5, beta = 1.1))
#' plot(fit)
#' 
#'
#' @rdname EXPPLindley
#' @export
dexpplindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0,  beta > 0)
  if(log)
  {
    t1 <- log(beta)
    t2 <- log(alpha)
    t3 <- log(theta)
    t5 <- 1 + theta
    t6 <- log(t5)
    t7 <- x ^ alpha
    t9 <- log1p(t7)
    t11 <- log(x)
    t13 <- t7 * theta
    t18 <- exp(-t13)
    t22 <- (1 - t18 * (1 + t7 / t5 * theta)) ^ (beta - 1)
    t23 <- log(t22)
    t1 + t2 + 2 * t3 - t6 + t9 + t11 * (alpha - 1) - t13 + t23
  }
  else
  {
    t2 <- theta ^ 2
    t4 <- 1 / (1 + theta)
    t7 <- x ^ alpha
    t10 <- x ^ (alpha - 1)
    t13 <- exp(-t7 * theta)
    t20 <- (1 - t13 * (t7 * t4 * theta + 1)) ^ (beta - 1)
    t20 * t13 * t10 * (1 + t7) * t4 * t2 * beta * alpha
  }
}

#' @rdname EXPPLindley
#' @export
pexpplindley <- function(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
	t4 <- q ^ alpha
	t8 <- exp(-t4 * theta)
	O  <- (1 - t8 * (1 + t4 / (1 + theta) * theta)) ^ beta
  }
  else
  {
	t4  <- q ^ alpha
	t8  <- exp(-t4 * theta)
	t11 <- (1 - t8 * (1 + t4 / (1 + theta) * theta)) ^ beta
	O   <- 1 - t11
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname EXPPLindley
#' @export
qexpplindley <- function(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(lower.tail)
  {
	t1  <- 0.1e1 / alpha
	t2  <- theta ^ t1
	t4  <- 1 + theta
	t5  <- log(p)
	t8  <- exp(0.1e1 / beta * t5)
	t11 <- exp(-t4)
	t13 <- W(t11 * (t8 - 1) * t4, branch = -1)
	t15 <- (-t13 - 1 - theta) ^ t1
	O   <- t15 / t2
  }
  else
  {
	t1  <- 0.1e1 / alpha
	t2  <- theta ^ t1
	t5  <- log1p(-p)
	t8  <- exp(0.1e1 / beta * t5)
	t10 <- 1 + theta
	t13 <- exp(-t10)
	t15 <- W(t13 * (t8 - 1) * t10, branch = -1)
	t17 <- exp(t15 + 1 + theta)
	t20 <- (-theta * t17 - theta * t8 - t17 - t8 + theta + 1) ^ t1
	t23 <- exp(t15 * t1)
	t25 <- exp(t1)
	t29 <- exp(theta * t1)
	O   <- 0.1e1 / t29 / t25 / t23 * t20 / t2
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname EXPPLindley
#' @export
rexpplindley <- function(n, theta, alpha, beta, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(mixture)
  {
    rplindley(n, theta, alpha, mixture = TRUE) ^ (1 / beta)
  }
  else
  {
    rplindley(n, theta, alpha, mixture = FALSE) ^ (1 / beta)
  }
}

#' @rdname EXPPLindley
#' @export
hexpplindley <- function(x, theta, alpha, beta, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, beta > 0)
  if(log)
  {
	t1 <- log(alpha)
	t2 <- log(theta)
	t4 <- 1 + theta
	t5 <- log(t4)
	t6 <- x ^ alpha
	t8 <- log1p(t6)
	t10 <- log(x)
	t12 <- t6 * theta
	t17 <- exp(-t12)
	t19 <- 1 - t17 * (1 + t6 / t4 * theta)
	t21 <- t19 ^ (beta - 1)
	t23 <- t19 ^ beta
	t27 <- log(0.1e1 / (1 - t23) * t21 * beta)
	t1 + 2 * t2 - t5 + t8 + t10 * (alpha - 1) - t12 + t27
  }
  else
  {
	t2 <- theta ^ 2
	t4 <- 1 / (1 + theta)
	t7 <- x ^ alpha
	t10 <- x ^ (alpha - 1)
	t13 <- exp(-t7 * theta)
	t18 <- 1 - t13 * (t7 * t4 * theta + 1)
	t20 <- t18 ^ (beta - 1)
	t22 <- t18 ^ beta
	0.1e1 / (1 - t22) * t20 * t13 * t10 * (1 + t7) * t4 * t2 * beta * alpha
  }
}
