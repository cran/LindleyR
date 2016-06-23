#' @name DWLindley
#' @aliases DWLindley ddwlindley pdwlindley qdwlindley rdwlindley
#'
#' @title Discrete Weighted Lindley Distribution
#'
#' @description Probability mass function, distribution function, quantile function and random number generation for the discrete weighted Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Ricardo P. de Oliveira \email{rpuziol.oliveira@gmail.com}

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
#' @param x,q vector of integer positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameter.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{ddwlindley} gives the probability mass function, \code{pdwlindley} gives the distribution function, \code{qdwlindley} gives the quantile function and \code{rdwlindley} generates random deviates.
#' @return Invalid arguments will return an error message.
#'
#' @seealso \code{\link[LindleyR]{WLindley}}.
#'
#' @source [d-p-q-r]dwlindley are calculated directly from the definitions. \code{rdwlindley} uses the discretize values.
#'
#' @details
#' Probability mass function
#' \deqn{P(X=x\mid \theta ,\alpha )=\frac{1}{\left( \theta +\alpha \right) \Gamma \left( \alpha \right) }\sum\limits_{i=0}^{1}\left( -1\right) ^{i}\left\{ \left( \theta +\alpha \right) \Gamma \left[ \alpha ,\theta \left( x+i\right) \right] +\left[\theta \left( x+i\right) \right] ^{\alpha }e^{-\theta \left( x+i\right)}\right\} }
#'
#' where \eqn{\Gamma \left(\alpha,\theta x\right) = \int_{\theta x}^{\infty}x^{\alpha -1}e^{-x}dx} is the upper incomplete gamma function.
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the one-parameter discrete Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rdwlindley(n = 1000, theta = 1.5, alpha = 1.5)
#' plot(table(x) / sum(table(x)))
#' points(unique(x),ddwlindley(unique(x), theta = 1.5, alpha = 1.5))
#'
#' ## fires in Greece data (from Bakouch et al., 2014)
#' data(fires)
#' library(fitdistrplus)
#' fit <- fitdist(fires, 'dwlindley', start = list(theta = 0.30, alpha = 1.0), discrete = TRUE)
#' gofstat(fit, discrete = TRUE)
#' plot(fit)
#'
#' @rdname DWLindley
#' @export
ddwlindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, x >= 0)
  if(!isTRUE(all(x == floor(x)))) stop("'x' must only contain nonnegative integers")
  if(log)
  {
    log(pwlindley(q = x, theta, alpha, lower.tail = FALSE) - pwlindley(q = x + 1, theta, alpha, lower.tail = FALSE))
  }
  else
  {
	  pwlindley(q = x, theta, alpha, lower.tail = FALSE) - pwlindley(q = x + 1, theta, alpha, lower.tail = FALSE)
  }
}

#' @rdname DWLindley
#' @export
pdwlindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    cdf <- sapply(1 : length(q), function(i) sum(ddwlindley(x = 0 : q[i], theta, alpha, log = FALSE)))
  }
  else
  {
    cdf <- 0.1e1 - sapply(1 : length(q), function(i) sum(x = ddwlindley(0 : q[i], theta, alpha, log = FALSE)))
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname DWLindley
#' @export
qdwlindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    qtf <- floor(qwlindley(p, theta, alpha, lower.tail = TRUE, log.p = FALSE))
  }
  else
  {
    qtf <- floor(qwlindley(p, theta, alpha, lower.tail = FALSE, log.p = FALSE))
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname DWLindley
#' @export
rdwlindley <- function(n, theta, alpha)
{
  qdwlindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
}
