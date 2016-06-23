#' @name DLindley
#' @aliases DLindley ddlindley pdlindley qdlindley rdlindley fires
#'
#' @title One-Parameter Discrete Lindley Distribution
#'
#' @description Probability mass function, distribution function, quantile function and random number generation for the one-parameter discrete Lindley distribution with parameter theta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Bakouch, H. S., Jazi, M. A. and Nadarajah, S. (2014). A new discrete distribution. \emph{Statistics: A Journal of Theoretical and Applied Statistics}, \bold{48}, 1, 200-240.
#'
#' Gomez-Deniz, E. and Calderín-Ojeda, E. (2013). The discrete Lindley distribution: properties and applications. \emph{Journal of Statistical Computation and Simulation}, \bold{81}, 11, 1405-1416.
#'
#' @param x,q vector of integer positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta positive parameter.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#' @return \code{ddlindley} gives the probability mass function, \code{pdlindley} gives the distribution function, \code{qdlindley} gives the quantile function and \code{rdlindley} generates random deviates.
#' @return Invalid arguments will return an error message.
#'
#' @seealso \code{\link[LindleyR]{Lindley}}.
#'
#' @source [d-p-q-r]dlindley are calculated directly from the definitions. \code{rdlindley} uses the discretize values.
#'
#' @details
#' Probability mass function
#' \deqn{P\left( X=x\mid \theta \right) =\sum\limits_{i=0}^{1}\left( -1\right) ^{i}\left( 1+\frac{\theta }{1+\theta }\left( x+i\right) \right) e^{-\theta \left( x+i\right) }}
#'
#' @examples
#' set.seed(1)
#' x <- rdlindley(n = 1000, theta = 1.5)
#' plot(table(x) / sum(table(x)))
#' points(unique(x),ddlindley(unique(x), theta = 1.5))
#'
#' ## fires in Greece data (from Bakouch et al., 2014)
#' data(fires)
#' library(fitdistrplus)
#' fit <- fitdist(fires, 'dlindley', start = list(theta = 0.30), discrete = TRUE)
#' gofstat(fit, discrete = TRUE)
#' plot(fit)
#'
#' @rdname DLindley
#' @export
ddlindley <- function(x, theta, log = FALSE)
{
  stopifnot(theta > 0, x >= 0)
  if(!isTRUE(all(x == floor(x)))) stop("'x' must only contain nonnegative integers")
  if(log)
  {
	  log(plindley(q = x, theta, lower.tail = FALSE) - plindley(q = x + 1, theta, lower.tail = FALSE))
  }
  else
  {
	  plindley(q = x, theta, lower.tail = FALSE) - plindley(q = x + 1, theta, lower.tail = FALSE)
  }
}

#' @rdname DLindley
#' @export
pdlindley <- function(q, theta, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    cdf <- sapply(1:length(q), function(i) sum(ddlindley(x = 0:q[i], theta, log = FALSE)))
  }
  else
  {
    cdf <- 0.1e1 - sapply(1:length(q), function(i) sum(x = ddlindley(0:q[i], theta, log = FALSE)))
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname DLindley
#' @export
qdlindley <- function(p, theta, lower.tail = TRUE, log.p = FALSE)
{
  if(lower.tail)
  {
    qtf <- floor(qlindley(p, theta, lower.tail = TRUE, log.p = FALSE))
  }
  else
  {
    qtf <- floor(qlindley(p, theta, lower.tail = FALSE, log.p = FALSE))
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname DLindley
#' @export
rdlindley <- function(n, theta)
{
  qdlindley(p = runif(n), theta, lower.tail = TRUE, log.p = FALSE)
}
