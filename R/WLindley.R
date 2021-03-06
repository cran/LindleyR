#' @importFrom stats rgamma runif rbinom uniroot pgamma integrate pexp
#'
#' @name WLindley
#' @aliases WLindley dwlindley pwlindley qwlindley rwlindley hwlindley carbonfibers
#'
#' @title Weighted Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the weighted Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @note The \code{\link[stats]{uniroot}} function with default arguments is used to find out the quantiles.
#'
#' @references
#'
#' Al-Mutairi, D. K., Ghitany, M. E., Kundu, D., (2015). Inferences on stress-strength reliability from weighted Lindley distributions. \emph{Communications in Statistics - Theory and Methods}, \bold{44}, (19), 4096-4113.
#'
#' Bashir, S., Rasul, M., (2015). Some properties of the weighted Lindley distribution. \emph{EPRA Internation Journal of Economic and Business Review}, \bold{3}, (8), 11-17.
#'
#' Ghitany, M. E., Alqallaf, F., Al-Mutairi, D. K. and Husain, H. A., (2011). A two-parameter weighted Lindley distribution and its applications to survival data. \emph{Mathematics and Computers in Simulation}, \bold{81}, (6), 1190-1201.
#'
#' Mazucheli, J., Louzada, F., Ghitany, M. E., (2013). Comparison of estimation methods for the parameters of the weighted Lindley distribution. \emph{Applied Mathematics and Computation}, \bold{220}, 463-471.
#'
#' Mazucheli, J., Coelho-Barros, E. A. and Achcar, J. (2016). An alternative reparametrization on the weighted Lindley distribution. \emph{Pesquisa Operacional}, (to appear).
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#' @param L,U interval which \code{uniroot} searches for a root (quantile), L = 1e-4 and U = 50 are the default values.
#'
#' @return \code{dwlindley} gives the density, \code{pwlindley} gives the distribution function, \code{qwlindley} gives the quantile function, \code{rwlindley} generates random deviates and \code{hwlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[lamW]{lambertWm1}}, \code{\link[stats]{uniroot}}, \code{\link[LindleyR]{DWLindley}}.
#'
#' @source [d-h-p-q-r]wlindley are calculated directly from the definitions. \code{rwlindley} uses either a two-component mixture of the gamma distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f( x\mid \theta,\alpha) =\frac{\theta ^{\alpha +1}}{\left( \theta+\alpha \right) \Gamma \left( \alpha \right) }x^{\alpha -1}\left( 1+x\right)e^{-\theta x}  \label{density-weighted-lindley}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha) =1 - \frac{\left( \theta +\alpha \right)\Gamma \left( \alpha,\theta x\right) +\left( \theta x\right) ^{\alpha}e^{-\theta x}}{\left( \theta +\alpha \right) \Gamma \left( \alpha \right) }}
#'
#' Quantile function
#' \deqn{\code{does not have a closed mathematical expression}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta,\alpha) =\frac{\theta ^{\alpha +1}x^{\alpha-1}\left( 1+x\right) e^{-\theta x}}{\left( \theta +\alpha \right) \Gamma\left( \alpha,\theta x\right) +\left( \theta x\right) ^{\alpha }e^{-\theta x}}}
#'
#' where \eqn{\Gamma \left(\alpha,\theta x\right) = \int_{\theta x}^{\infty}x^{\alpha -1}e^{-x}dx} is the upper incomplete gamma function.
#'
#' \bold{Particular case:} \eqn{\alpha=1} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rwlindley(n = 1000, theta = 1.5, alpha = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dwlindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pwlindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pwlindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qwlindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qwlindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' ## carbon fibers data (from Ghitany et al., 2013)
#' data(carbonfibers)
#' library(fitdistrplus)
#' fit <- fitdist(carbonfibers, 'wlindley', start = list(theta = 0.1, alpha = 0.1))
#' plot(fit)
#'
#' @rdname WLindley
#' @export
dwlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t2  <- log(theta)
    t5  <- log(alpha + theta)
    t7  <- log1p(x)
    t9  <- log(x)
    t12 <- gamma(alpha)
    t13 <- log(t12)
    (alpha + 1) * t2 - t5 + t7 + (alpha - 1) * t9 - theta * x - t13
  }
  else
  {
    t2  <- theta ^ (alpha + 1)
    t9  <- x ^ (alpha - 1)
    t11 <- exp(-theta * x)
    t13 <- gamma(alpha)
    t2 / (alpha + theta) * (1 + x) * t9 * t11 / t13
  }
}

#' @rdname WLindley
#' @export
pwlindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  p <- theta / (theta + alpha)
  if(lower.tail)
  {
    cdf <- p * pgamma(q, shape = alpha, rate = theta, lower.tail = TRUE) + (1 - p) * pgamma(q, shape = alpha + 1, rate = theta, lower.tail = TRUE)
  }
  else
  {
    cdf <- p * pgamma(q, shape = alpha, rate = theta, lower.tail = FALSE) + (1 - p) * pgamma(q, shape = alpha + 1, rate = theta, lower.tail = FALSE)
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname WLindley
#' @export
qwlindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pwlindley(q, theta, alpha, lower.tail = TRUE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
    if(log.p) return(log(qtf)) else return(qtf)
  }
  else
  {
    fx <- function(p)
    {
      tryCatch(uniroot(function(q) pwlindley(q, theta, alpha, lower.tail = FALSE, log.p = FALSE) - p, lower = L, upper = U)$root, error = function(e) NaN)
    }
    qtf <- sapply(p, fx)
    if(log.p) return(log(qtf)) else return(qtf)
  }
}

#' @rdname WLindley
#' @export
rwlindley <- function(n, theta, alpha, mixture = TRUE, L = 1e-4, U = 50)
{
  stopifnot(theta > 0, alpha > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (theta + alpha))
    x <- p * rgamma(n, shape = alpha, rate = theta) + (1 - p) * rgamma(n, shape = alpha + 1, rate = theta)
  }
  else
  {
    x  <- qwlindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE, L, U)
	  is <- is.na(x)
	  na <- any(is)
	  if(na)
	  {
		  y <- which(is)
    	i <- 1
		  l <- length(y)
    	while(i <= l)
		  {
        x[y[i]] <- qwlindley(p = runif(1), theta, alpha, lower.tail = TRUE, log.p = FALSE, L, U)
        if(!is.na(x[y[i]])) i <- i + 1
    	}
    }
	}
	x
}

#' @rdname WLindley
#' @export
hwlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t2  <- log(theta)
    t5  <- log1p(x)
    t7  <- log(x)
    t9  <- theta * x
    # upper incomplete gamma function
    t11 <- gamma(alpha) * pgamma(q = t9, shape = alpha, rate = 1, lower.tail = FALSE)
    t13 <- t9 ^ alpha
    t14 <- exp(-t9)
    t18 <- log(0.1e1 / ((alpha + theta) * t11 + t13 * t14))
    (alpha + 1) * t2 + t5 + (alpha - 1) * t7 - t9 + t18
  }
  else
  {
    t2 <- theta ^ (alpha + 1)
    t6 <- x ^ (alpha - 1)
    t7 <- theta * x
    t8 <- exp(-t7)
    # upper incomplete gamma function
    t11 <- gamma(alpha) * pgamma(q = t7, shape = alpha, rate = 1, lower.tail = FALSE)
    t13 <- t7 ^ alpha
    t2 * (1 + x) * t6 * t8 / ((alpha + theta) * t11 + t13 * t8)
  }
}
