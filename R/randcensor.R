#' @name randcensor
#' @aliases randcensor
#' @title Generate Independent Random Censored Lifetime
#'
#' @description Implements a function to draw censored random samples, with a desired censoring rate, when the event times are any continuous lifetime distribution supported by R. The one-parameter Lindley, uniform and exponential are the distributions that can be used as the censoring distributions.
#'
#' @note Finds the parameter of the censoring distribution using \code{\link[stats]{integrate}} and \code{\link[stats]{uniroot}}.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @references
#'
#' Klein, J. P., Moeschberger, M. L., (2003). \emph{Survival Analysis: Techniques for Censored and Truncated Data, 2nd Edition}. Springer-Verlag, New York.
#'
#' Lawless, J. F., (2003). \emph{Statistical models and methods for lifetime data, 2nd Edition}. Wiley Series in Probability and Statistics. John Wiley & Sons, Hoboken, NJ.
#'
#' Meeker, W. Q., Escobar, L. A., (1998). \emph{Statistical Methods for Reliability Data}. John Wiley and Sons, New York.
#'
#' @param n number o generated observations.
#' @param pcens desired censoring rate.
#' @param timedistr a character string \code{'name'} naming a lifetime distribution for which the corresponding density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding random function \code{qname}  are defined. The one-parameter \code{Lindley} distribution is taken as default.
#' @param censordistr a character string \code{'name'} naming the censoring distribution. 'lindley' (default) for the one-parameter Lindley distribution, 'exp' for the exponential distribution and 'unif' for the uniform distribution.
#' @param ... parameters that define the event time distribution (timedistr). Must be provided in the same way as it is in built in R functions.
#'
#' @return \code{randcensor} returns a list with the timedistr distribution, the censordistr distribution, the calculated parameter of the censordist distribution and n observations which is either the lifetime (delta = 1) or a censored lifetime (delta = 0).
#' @return Invalid arguments will return an error message.
#'
#' @seealso \code{\link[stats]{Distributions}}, \code{\link[fitdistrplus]{fitdistcens}}, \code{\link[stats]{integrate}}, \code{\link[LindleyR]{Lindley}}, \code{\link[stats]{uniroot}}.
#'
#' @examples 
#' x <- randcensor(n = 100, pcens = 0.2, timedistr = 'lindley', censordistr = 'lindley',
#'  theta = 1.5)
#' table(x$data['delta']) / 100
#'
#
#' x <- randcensor(n = 100, pcens = 0.2, timedistr = 'wlindley', censordistr = 'lindley',
#'  theta = 1.5, alpha = 0.5)
#' table(x$data['delta']) / 100
#'
#' x <- randcensor(n = 100, pcens = 0.2, timedistr = 'weibull', censordistr = 'lindley',
#'  shape = 0.5, scale = 1.5)
#' table(x$data['delta']) / 100
#'
#' x <- randcensor(n = 100, pcens = 0.2, timedistr = 'lnorm', censordistr = 'unif',
#'  meanlog = 1, sdlog = 1)
#' table(x$data['delta']) / 100
#'
#'
#' @rdname randcensor
#' @export
randcensor <- function(n, pcens = 0.1, timedistr = 'lindley', censordistr = 'lindley', ...)
{
  stopifnot(n >= 1, pcens < 1, pcens > 0)
  dtime <-   function(z, timedistr, ...) {do.call(paste("d", timedistr, sep = ""), list(z, ...))}
  ptime <-   function(z, timedistr, ...) {do.call(paste("p", timedistr, sep = ""), list(z, ...))}
  rtime <-   function(Z, timedistr, ...) {do.call(paste("r", timedistr, sep = ""), list(n, ...))}

  if(censordistr == 'exp')
  {
    fxcE <- function(parcens, timedistr, pcens, ...)
    {
      integrate(f  = function(z, parcens, timedistr){ pexp(q = z, rate = 1 / parcens) * dtime(z, timedistr, ...)}, lower = 0, upper = Inf, parcens, timedistr)$value - pcens
    }
    theta  <- uniroot(f = fxcE, interval = c(0, 100), timedistr, pcens, ...)$root
    x.cens <- rtime(n, censordistr, rate = 1 / theta)
  }

  else if(censordistr == 'unif')
  {
    fxcU   <- function(parcens, timedistr, pcens, ...)
    {
      integrate(f  = function(z, timedistr){ptime(z, timedistr, lower.tail = FALSE, ...)}, lower = 0, upper = parcens, timedistr)$value / parcens - pcens
    }
    theta  <- uniroot(f = fxcU, interval = c(1e-4, 100), timedistr, pcens, ...)$root
    x.cens <- rtime(n, censordistr, min = 0, max = theta)
  }

  else if(censordistr == 'lindley')
  {
    fxcL   <- function(parcens, timedistr, pcens, ...)
    {
      integrate(f  = function(z, parcens, timedistr){ plindley(q = z, theta = parcens) * dtime(z, timedistr, ...)}, lower = 0, upper = Inf, parcens, timedistr)$value - pcens
    }
    theta  <- uniroot(f = fxcL, interval = c(1e-4, 100), timedistr, pcens, ...)$root
    x.cens <- rtime(n, censordistr, theta)
  }
  else stop('distribution not available')

  x.time <-  rtime(n, timedistr, ...)
  data   <-  data.frame(time = pmin(x.cens, x.time), delta = ifelse(x.time <= x.cens, 1, 0))
  list(timedistr = timedistr, censordistr = censordistr, parcens = theta, data = data)
}
