#' @name DPLindley
#' @aliases DPLindley ddplindley pdplindley qdplindley rdplindley
#'
#' @title Discrete Power Lindley Distribution
#'
#' @description Probability mass function, distribution function, quantile function and random number generation for the discrete power Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Ricardo P. de Oliveira \email{rpuziol.oliveira@gmail.com}

#' @references
#'
#' Ghitany, M. E., Al-Mutairi, D. K., Balakrishnan, N. and Al-Enezi, L. J., (2013). Power Lindley distribution and associated inference. \emph{Computational Statistics and Data Analysis}, \bold{64}, 20-33.
#'
#' Mazucheli, J., Ghitany, M. E. and Louzada, F., (2013). Power Lindley distribution: Diferent methods of estimation and their applications to survival times data. \emph{Journal of Applied Statistical Science}, \bold{21}, (2), 135-144.
#'
#' @param x,q vector of integer positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameter.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @return \code{ddplindley} gives the probability mass function, \code{pdplindley} gives the distribution function, \code{qdplindley} gives the quantile function and \code{rdplindley} generates random deviates.
#' @return Invalid arguments will return an error message.
#'
#' @seealso \code{\link[LindleyR]{PLindley}}.
#'
#' @source [d-p-q-r]dplindley are calculated directly from the definitions. \code{rdplindley} uses the discretize values.
#'
#' @details
#' Probability mass function
#' \deqn{P(X=x\mid \theta ,\alpha )=\sum\limits_{i=0}^{1}\left( -1\right) ^{i}\left( 1+{\frac{\theta }{\theta +1}}\left( x+i\right) ^{\alpha }\right) \ e^{-\theta \left( x+i\right) ^{\alpha}}}
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the one-parameter discrete Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rdplindley(n = 1000, theta = 1.5, alpha = 0.5)
#' plot(table(x) / sum(table(x)))
#' points(unique(x),ddplindley(unique(x), theta = 1.5, alpha = 0.5))
#'
#' ## fires in Greece data (from Bakouch et al., 2014)
#' data(fires)
#' library(fitdistrplus)
#' fit <- fitdist(fires, 'dplindley', start = list(theta = 0.30, alpha = 1.0), discrete = TRUE)
#' gofstat(fit, discrete = TRUE)
#' plot(fit)
#'
#' @rdname DPLindley
#' @export
ddplindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0, x >= 0)
  if(!isTRUE(all(x == floor(x)))) stop("'x' must only contain nonnegative integers")
  if(log)
  {
	  log(pplindley(q = x, theta, alpha, lower.tail = FALSE) - pplindley(q = x + 1, theta, alpha, lower.tail = FALSE))
  }
  else
  {
	  pplindley(q = x, theta, alpha, lower.tail = FALSE) - pplindley(q = x + 1, theta, alpha, lower.tail = FALSE)
  }
}

#' @rdname DPLindley
#' @export
pdplindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    cdf <- sapply(1 : length(q), function(i) sum(ddplindley(x = 0 : q[i], theta, alpha, log = FALSE)))
  }
  else
  {
    cdf <- 0.1e1 - sapply(1 : length(q), function(i) sum(x = ddplindley(x = 0 : q[i], theta, alpha, log = FALSE)))
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname DPLindley
#' @export
qdplindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    qtf <- floor(qplindley(p, theta, alpha, lower.tail = TRUE, log.p = FALSE))
  }
  else
  {
    qtf <- floor(qplindley(p, theta, alpha, lower.tail = FALSE, log.p = FALSE))
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname DPLindley
#' @export
rdplindley <- function(n, theta, alpha)
{
  qdplindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
}
