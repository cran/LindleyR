% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ILindley.R
\name{ILindley}
\alias{ILindley}
\alias{dilindley}
\alias{hilindley}
\alias{pilindley}
\alias{qilindley}
\alias{rilindley}
\title{Inverse Lindley Distribution}
\source{
[d-h-p-q-r]ilindley are calculated directly from the definitions. \code{rilindley} uses either a two-component mixture of inverse gamma distributions or the quantile function.
}
\usage{
dilindley(x, theta, log = FALSE)

pilindley(q, theta, lower.tail = TRUE, log.p = FALSE)

qilindley(p, theta, lower.tail = TRUE, log.p = FALSE)

rilindley(n, theta, mixture = TRUE)

hilindley(x, theta, log = FALSE)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{theta}{positive parameter.}

\item{log, log.p}{logical; If TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}

\item{mixture}{logical; If TRUE, (default), random deviates are generated from a two-component mixture of inverse-gamma distributions, otherwise from the quantile function.}
}
\value{
\code{dilindley} gives the density, \code{pilindley} gives the distribution function, \code{qilindley} gives the quantile function, \code{rilindley} generates random deviates and \code{hilindley} gives the hazard rate function.

Invalid arguments will return an error message.
}
\description{
Density function, distribution function, quantile function, random number generation and hazard rate function for the inverse Lindley distribution with parameter theta.
}
\details{
Probability density function
\deqn{f(x\mid \theta )=\frac{\theta ^{2}}{1+\theta }\left( \frac{1+x}{x^{3}}\right) e^{-\frac{\theta }{x}}}

Cumulative distribution function
\deqn{F(x\mid \theta )=\left( 1+\frac{\theta }{x\left( 1+\theta \right) }\right) {e{^{-{\frac{\theta }{x}}}}}}

Quantile function
\deqn{Q(p\mid \theta) =-\left[ 1+\frac{1}{\theta }+\frac{1}{\theta }W_{-1}\left( -p\left( 1+\theta \right) e{^{-\left( 1+\theta \right) }} \right) \right] ^{-1}}

Hazard rate function
\deqn{h(x\mid \theta )=\frac{\theta ^{2}\left( 1+x\right) {e{^{-{\frac{\theta }{x}}}}}}{x^{3}\left( 1+\theta \right) \left[ 1-\left( 1+\frac{\theta }{x\left(1+\theta \right) }\right) {e{^{-{\frac{\theta }{x}}}}}\right] }}

where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.
}
\examples{
x <- seq(from = 0.1, to = 3, by = 0.05)
plot(x, dilindley(x, theta = 1.0), xlab = 'x', ylab = 'pdf')

p <- seq(from = 0.1, to = 0.9, by = 0.1)
q <- quantile(x, prob = p)
pilindley(q, theta = 1.5, lower.tail = TRUE)
pilindley(q, theta = 1.5, lower.tail = FALSE)
qilindley(p, theta = 1.5, lower.tail = TRUE)
qilindley(p, theta = 1.5, lower.tail = FALSE)

set.seed(1)
x <- rilindley(n = 100, theta = 1.0)
library(fitdistrplus)
fit <- fitdist(x, 'ilindley', start = list(theta = 1.0))
plot(fit)

}
\author{
Josmar Mazucheli \email{jmazucheli@gmail.com}

Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
}
\references{
Sharma, V. K., Singh, S. K., Singh, U., Agiwal, V., (2015). The inverse Lindley distribution: a stress-strength reliability model with application to head and neck cancer data. \emph{Journal of Industrial and Production Engineering}, \bold{32}, (3), 162-173.
}
\seealso{
\code{\link[lamW]{lambertWm1}}, \code{\link[actuar]{rinvgamma}}.
}

