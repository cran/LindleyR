#' @importFrom lamW lambertWm1
#'
#' @name PLindley
#' @aliases PLindley dplindley pplindley qplindley rplindley hplindley carbonfibres
#'
#' @title Power Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random number generation and hazard rate function for the power Lindley distribution with parameters theta and alpha.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Ghitany, M. E., Al-Mutairi, D. K., Balakrishnan, N. and Al-Enezi, L. J., (2013). Power Lindley distribution and associated inference. \emph{Computational Statistics and Data Analysis}, \bold{64}, 20-33.
#'
#' Mazucheli, J., Ghitany, M. E. and Louzada, F., (2013). Power Lindley distribution: Diferent methods of estimation and their applications to survival times data. \emph{Journal of Applied Statistical Science}, \bold{21}, (2), 135-144.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta,alpha positive parameters.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical; If TRUE, (default), random deviates are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dplindley} gives the density, \code{pplindley} gives the distribution function, \code{qplindley} gives the quantile function, \code{rplindley} generates random deviates and \code{hplindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso  \code{\link[lamW]{lambertWm1}}, \code{\link[LindleyR]{DPLindley}}.
#'
#' @source [d-h-p-q-r]plindley are calculated directly from the definitions. \code{rplindley} uses either a two-component mixture of gamma distributions or the quantile function.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta,\alpha )={\frac{\alpha \theta ^{2}}{1 + \theta}}(1+x^{\alpha})\ x^{\alpha -1}\ e^{-\theta x^{\alpha }}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta,\alpha )=1-\left( 1+{\frac{\theta }{1 + \theta}}x^{\alpha }\right) \ e^{-\theta x^{\alpha }}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta,\alpha )=\left( -1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left( \left( 1+\theta \right) \left(p-1\right) e^{-(1+\theta) }\right) \right) ^{\frac{1}{\alpha }}}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta ,\alpha )={\frac{\alpha \theta ^{2}(1+x^{\alpha })x^{\alpha-1}}{\left( \theta +1\right) \left( 1+{\frac{\theta }{\theta +1}}x^{\alpha }\right) }} }
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' \bold{Particular case:} \eqn{\alpha = 1} the one-parameter Lindley distribution.
#'
#' @examples
#' set.seed(1)
#' x <- rplindley(n = 1000, theta = 1.5, alpha = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dplindley(S, theta = 1.5, alpha = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' pplindley(q, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' pplindley(q, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#' qplindley(p, theta = 1.5, alpha = 1.5, lower.tail = TRUE)
#' qplindley(p, theta = 1.5, alpha = 1.5, lower.tail = FALSE)
#'
#' ## carbon fibers data (from Ghitany et al., 2013)
#' data(carbonfibers)
#' library(fitdistrplus)
#' fit <- fitdist(carbonfibers, 'plindley', start = list(theta = 0.1, alpha = 0.1))
#' plot(fit)
#'
#' @rdname PLindley
#' @export
dplindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t5 <- log1p(theta)
    t6 <- x ^ alpha
    t8 <- log1p(t6)
    t10 <- log(x)
    t1 + 2 * t2 - t5 + t8 + (alpha - 1) * t10 - theta * t6
  }
  else
  {
    t1 <- theta ^ 2
    t6 <- x ^ alpha
    t9 <- x ^ (alpha - 1)
    t12 <- exp(-theta * t6)
    alpha * t1 / (1 + theta) * (1 + t6) * t9 * t12
  }
}

#' @rdname PLindley
#' @export
pplindley <- function(q, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
    t4 <- q ^ alpha
    t8 <- exp(-theta * t4)
    cdf<- 1 - (1 + theta / (1 + theta) * t4) * t8
  }
  else
  {
    t4 <- q ^ alpha
    t8 <- exp(-theta * t4)
    cdf<- (1 + theta / (1 + theta) * t4) * t8
  }
  if(log.p) return(log(cdf)) else return(cdf)
}

#' @rdname PLindley
#' @export
qplindley <- function(p, theta, alpha, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(lower.tail)
  {
    t1  <- 0.1e1 / alpha
    t2  <- theta ^ t1
    t4  <- 0.1e1 + theta
    t7  <- exp(-t4)
    t9  <- lambertWm1(t4 * (p - 0.1e1) * t7)
    t11 <- (-t9 - 1 - theta) ^ t1
    qtf <- 0.1e1 / t2 * t11
  }
  else
  {
    t1  <- 0.1e1 / alpha
    t2  <- theta ^ t1
    t4  <- 1 + theta
    t6  <- exp(-t4)
    t8  <- lambertWm1(-p * t4 * t6)
    t10 <- (-t8 - 0.1e1 - theta) ^ t1
    qtf <- 0.1e1 / t2 * t10
  }
  if(log.p) return(log(qtf)) else return(qtf)
}

#' @rdname PLindley
#' @export
rplindley <- function(n, theta, alpha, mixture = TRUE)
{
  stopifnot(theta > 0, alpha > 0)
  if(mixture)
  {
    rlindley(n, theta, mixture = TRUE) ^ (1 / alpha)
  }
  else
  {
    qplindley(p = runif(n), theta, alpha, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname PLindley
#' @export
hplindley <- function(x, theta, alpha, log = FALSE)
{
  stopifnot(theta > 0, alpha > 0)
  if(log)
  {
    t1 <- log(alpha)
    t2 <- log(theta)
    t4 <- 1 + theta
    t5 <- log(t4)
    t6 <- x ^ alpha
    t8 <- log1p(t6)
    t10 <- log(x)
    t16 <- log1p(theta / t4 * t6)
    t1 + 2 * t2 - t5 + t8 + (alpha - 1) * t10 - t16
  }
  else
  {
    t1 <- theta ^ 2
    t4 <- 1 / (1 + theta)
    t6 <- x ^ alpha
    t9 <- x ^ (alpha - 1)
    alpha * t1 * t4 * (1 + t6) * t9 / (theta * t4 * t6 + 1)
  }
}

