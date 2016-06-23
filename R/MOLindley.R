#' @importFrom lamW lambertWm1
#'
#' @name MOLindley
#' @aliases MOLindley dmolindley hmolindley pmolindley qmolindley rmolindley bladdercancer
#'
#' @title Marshall-Olkin Extended Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the Marshall-Olkin extended Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' do Espirito Santo, A. P. J., Mazucheli, J., (2015). Comparison of estimation methods for the Marshall-Olkin extended Lindley distribution. \emph{Journal of Statistical Computation and Simulation}, \bold{85}, (17), 3437-3450.
#'
#' Ghitany, M. E., Al-Mutairi, D. K., Al-Awadhi, F. A. and Al-Burais, M. M., (2012). Marshall-Olkin extended Lindley distribution and its application. \emph{International Journal of Applied Mathematics}, \bold{25}, (5), 709-721.
#'
#' Marshall, A. W., Olkin, I. (1997). A new method for adding a parameter to a family of distributions with application to the exponential and Weibull families. \emph{Biometrika}, \bold{84}, (3), 641.652.

#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{dmolindley} gives the density, \code{pmolindley} gives the distribution function, \code{qmolindley} gives the quantile function, \code{rmolindley} generates random deviates and \code{hmolindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[lamW]{lambertWm1}}, \code{\link[LindleyR]{Lindley}}.
#'
#' @source [d-h-p-q-r]molindley are calculated directly from the definitions. \code{rmolindley} uses the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha) =\frac{\alpha \theta^{2}(1+x)e^{-\theta x}}{(1+\theta )\left[ 1-\overline{\alpha }\left( 1+\frac{\theta x}{1+\theta }\right) e^{-\theta x}\right] ^{2}}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha) =1-\frac{\alpha \left( 1+\frac{\theta x}{1+\theta }\right) e^{-\theta x}}{1-\overline{\alpha }\left( 1+\frac{\theta x}{1+\theta }\right) e^{-\theta x}}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta,\alpha )=-1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left( \frac{(\theta +1)}{e^{1 + \theta}}\frac{(p-1)}{\left( 1-\overline{\alpha }p\right) }\right)}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha) =\frac{\theta ^{2}\left( 1+x\right) }{\left( 1+\theta +\theta x\right) \left[ 1-\overline{\alpha }\left( 1+\frac{\theta x}{1+\theta }\right) e^{-\theta x}\right] }}
#'
#' where \eqn{\overline{\alpha}=(1 - \alpha)} and \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular case:} \eqn{\alpha=1} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rmolindley(n = 1000, theta = 5, alpha = 5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dmolindley(S, theta = 5, alpha = 5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pmolindley(q, theta = 5, alpha = 5, lower.tail = TRUE)
#' pmolindley(q, theta = 5, alpha = 5, lower.tail = FALSE)
#' qmolindley(p, theta = 5, alpha = 5, lower.tail = TRUE)
#' qmolindley(p, theta = 5, alpha = 5, lower.tail = FALSE)
#'
#' ## bladder cancer data (from Warahena-Liyanage and Pararai, 2014)
#' data(bladdercancer)
#' library(fitdistrplus)
#' fit <- fitdist(bladdercancer, 'molindley', start = list(theta = 0.1, alpha =  1.0))
#' plot(fit)
#'
#' @rdname MOLindley
#' @export
dmolindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(alpha)
	t2 <- log(theta)
	t4 <- 1 + theta
	t5 <- log(t4)
	t7 <- log1p(x)
	t8 <- theta * x
	t14 <- exp(-t8)
	t17 <- (1 - t14 * (1 + 0.1e1 / t4 * t8) * (1 - alpha)) ^ 2
	t19 <- log(0.1e1 / t17)
	t1 + 2 * t2 - t5 + t7 - t8 + t19
  }
  else
  {
	t1 <- theta ^ 2
	t4 <- 1 / (1 + theta)
	t7 <- theta * x
	t8 <- exp(-t7)
	t16 <- (1 - t8 * (t4 * t7 + 1) * (1 - alpha)) ^ 2
	0.1e1 / t16 * t8 * (1 + x) * t4 * alpha * t1
  }
}

#' @rdname MOLindley
#' @export
pmolindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1 	<- theta * q
	t5 	<- 1 + t1 / (1 + theta)
	t7 	<- exp(-t1)
	cdf <- 1 - alpha * t5 * t7 / (1 - (1 - alpha) * t5 * t7)
  }
  else
  {
	t1  <- theta * q
	t5  <- 1 + t1 / (1 + theta)
	t7  <- exp(-t1)
	cdf <- alpha * t5 * t7 / (1 - (1 - alpha) * t5 * t7)
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname MOLindley
#' @export
qmolindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
	t1  <- 1 + theta
	t4  <- exp(-t1)
	t10 <- lambertWm1(t1 * (p - 1) * t4 / (alpha * p - p + 1))
	qtf <- -(t10 + theta + 1) / theta
  }
  else
  {
	t1  <- 1 + theta
	t3  <- exp(-t1)
	t9  <- lambertWm1(t1 * p * t3 / (alpha * p - alpha - p))
	qtf <- -(theta + t9 + 1) / theta
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname MOLindley
#' @export
rmolindley <- function(n, theta, alpha)
{
  stopifnot(theta > 0, alpha > 0)
  qmolindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
}

#' @rdname MOLindley
#' @export
hmolindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
	t1 <- log(theta)
	t3 <- 1 + theta
	t4 <- log(t3)
	t6 <- log1p(x)
	t7 <- theta * x
	t10 <- 1 + 0.1e1 / t3 * t7
	t11 <- log(t10)
	t14 <- exp(-t7)
	t18 <- log(0.1e1 / (1 - t14 * t10 * (1 - alpha)))
	2 * t1 - t4 + t6 - t11 + t18
  }
  else
  {
	t1 <- theta ^ 2
	t3 <- 1 / (1 + theta)
	t7 <- theta * x
	t9 <- t3 * t7 + 1
	t11 <- exp(-t7)
	1 / t9 / (1 - t11 * t9 * (1 - alpha)) * (1 + x) * t3 * t1
  }
}
