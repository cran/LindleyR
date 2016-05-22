#' @importFrom LambertW W
#
#' @name Lindley
#' @aliases Lindley dlindley plindley qlindley rlindley hlindley
#'
#' @title One-Parameter Lindley Distribution
#'
#' @description Density function, distribution function, quantile function, random numbers generation and hazard rate function for the one-parameter Lindley distribution with parameter theta.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
#'
#' @references
#'
#' Ghitany, M. E., Atieh, B., Nadarajah, S., (2008). Lindley distribution and its application. \emph{Mathematics and Computers in Simulation}, \bold{78}, 4, 49-506.
#'
#' Jodra, P., (2010). Computer generation of random variables with Lindley or Poisson-Lindley distribution via the Lambert W function. \emph{Mathematics and Computers in Simulation}, \bold{81}, (4), 851-859.
#'
#' Lindley, D. V., (1958). Fiducial distributions and Bayes' theorem. \emph{Journal of the Royal Statistical Society. Series B. Methodological}, \bold{20}, 102-107.
#'
#' Lindley, D. V., (1965). \emph{Introduction to Probability and Statistics from a Bayesian View-point, Part II: Inference}. Cambridge University Press, New York.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param theta positive parameter.
#' @param log,log.p logical. If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical. If TRUE (default) \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#' @param mixture logical. If TRUE (default), random values are generated from a two-component mixture of gamma distributions, otherwise from the quantile function.
#'
#' @return \code{dlindley} gives the density, \code{plindley} gives the distribution function, \code{qlindley} gives the quantile function, \code{rlindley} generates random deviates and \code{hlindley} gives the hazard rate function.
#' @return Invalid arguments will return an error message.
#
#' @seealso \code{\link[LambertW]{W}}, \code{\link[VGAM]{Lindley}}.
#'
#' @source [dpqh]lindley are calculated directly from the definitions. \code{rlindley}  uses either a two-component mixture of the gamma distributions or the inverse transform method.
#'
#' @details
#' Probability density function
#' \deqn{f(x\mid \theta )=\frac{\theta ^{2}}{(1+\theta )}(1+x)e^{-\theta x}}
#'
#' Cumulative distribution function
#' \deqn{F(x\mid \theta ) =1 - \frac{1+\theta +\theta x}{1+\theta }e^{-\theta x}}
#'
#' Quantile function
#' \deqn{Q(p\mid \theta )=-1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left((1+\theta)( p-1)e^{-1-\theta }\right)}
#'
#' Hazard rate function
#' \deqn{h(x\mid \theta )=\frac{\theta ^{2}}{1+\theta +\theta x}(1+x)}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
#'
#' @examples 
#' set.seed(1)
#' x <- rlindley(n = 1000, theta = 1.5, mixture = TRUE)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.1)
#' plot(S, dlindley(S, theta = 1.5), xlab = 'x', ylab = 'pdf')
#' hist(x, prob = TRUE, main = '', add = TRUE)
#'
#' p <- seq(from = 0.1, to = 0.9, by = 0.1)
#' q <- quantile(x, prob = p)
#' plindley(q, theta = 1.5, lower.tail = TRUE)
#' plindley(q, theta = 1.5, lower.tail = FALSE)
#' qlindley(p, theta = 1.5, lower.tail = TRUE)
#' qlindley(p, theta = 1.5, lower.tail = FALSE)
#'
#' library(fitdistrplus)
#' fit <- fitdist(x, 'lindley', start = list(theta = 1.5), lower = c(0))
#' plot(fit)
#'
#'
#' @rdname Lindley
#' @export
dlindley <- function(x, theta, log = FALSE)
{
  stopifnot(theta > 0)
  if(log)
  {
    t1 <- log(theta)
    t4 <- log1p(theta)
    t6 <- log1p(x)
    -theta * x + 2 * t1 - t4 + t6
  }
  else
  {
    t1 <- theta ^ 2
    t7 <- exp(-theta * x)
    t1 / (1 + theta) * (1 + x) * t7
  }
}

#' @rdname Lindley
#' @export
plindley <- function(q, theta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0)
  if(lower.tail)
  {
    t1 <- theta * q
    t6 <- exp(-t1)
    O  <- 1 - (t1 + theta + 1) / (1 + theta) * t6
  }
  else
  {
    t1 <- theta * q
    t6 <- exp(-t1)
    O  <- (t1 + theta + 1) / (1 + theta) * t6
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname Lindley
#' @export
qlindley <- function(p, theta, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(theta > 0)
  if(lower.tail)
  {
    t1 <- 1 + theta
    t4 <- exp(-t1)
    t6 <- W(t1 * (p - 1) * t4, branch = -1)
    O <- -(t6 + 1 + theta) / theta
  }
  else
  {
    t1 <- 1 + theta
    t3 <- exp(-t1)
    t5 <- W(-p * t1 * t3, branch = -1)
    O <- -(t5 + 1 + theta) / theta
  }
  if(log.p) return(log(O)) else return(O)
}

#' @rdname Lindley
#' @export
rlindley <- function(n, theta, mixture = TRUE)
{
  stopifnot(theta > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (1 + theta))
    p * rgamma(n, shape = 1, rate = theta) + (1 - p) * rgamma(n, shape = 2, rate = theta)
  }
  else
  {
    qlindley(p = runif(n), theta, lower.tail = TRUE, log.p = FALSE)
  }
}

#' @rdname Lindley
#' @export
hlindley <- function(x, theta, log = FALSE)
{
  stopifnot(theta > 0)
  if(log)
  {
    t1 <- log(theta)
    t5 <- log1p(theta * x + theta)
    t7 <- log1p(x)
    2 * t1 - t5 + t7
  }
  else
  {
    t1 <- theta ^ 2
    t1 / (theta * x + theta + 1) * (1 + x)
  }
}



