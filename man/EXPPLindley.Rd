% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EXPPLindley.R
\name{EXPPLindley}
\alias{EXPPLindley}
\alias{dexpplindley}
\alias{hexpplindley}
\alias{pexpplindley}
\alias{qexpplindley}
\alias{rexpplindley}
\title{Exponentiated Power Lindley Distribution}
\source{
[d-h-p-q-r]expplindley are calculated directly from the definitions. \code{rexpplindley} uses the quantile function.
}
\usage{
dexpplindley(x, theta, alpha, beta, log = FALSE)

pexpplindley(q, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)

qexpplindley(p, theta, alpha, beta, lower.tail = TRUE, log.p = FALSE)

rexpplindley(n, theta, alpha, beta)

hexpplindley(x, theta, alpha, beta, log = FALSE)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{theta, alpha, beta}{positive parameters.}

\item{log, log.p}{logical; If TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\value{
\code{dexpplindley} gives the density, \code{pexpplindley} gives the distribution function, \code{qexpplindley} gives the quantile function, \code{rexpplindley} generates random deviates and \code{hexpplindley} gives the hazard rate function.

Invalid arguments will return an error message.
}
\description{
Density function, distribution function, quantile function, random number generation and hazard rate function for the exponentiated power Lindley distribution with parameters theta, alpha and beta.
}
\details{
Probability density function
\deqn{f(x\mid \theta,\alpha,\beta )={\frac{\beta \alpha \theta ^{2}}{1 + \theta}}(1+x^{\alpha })x^{\alpha -1}e^{-\theta x^{\alpha }}\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) e^{-\theta x^{\alpha }}\right]^{\beta -1}}

Cumulative distribution function
\deqn{F(x\mid \theta,\alpha,\beta )=\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) \ e^{-\theta x^{\alpha }}\right] ^{\beta }}

Quantile function
\deqn{Q(p\mid \theta ,\alpha ,\beta )=\left( -1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left( \left( 1+\theta \right) \left( p^{\frac{1}{\beta }}-1\right)e^{-\left( 1+\theta \right) }\right) \right) ^{\frac{1}{\alpha }}}

Hazard rate function
\deqn{h(x\mid\theta,\alpha,\beta )={\frac{\beta \alpha \theta ^{2}(1+x^{\alpha })x^{\alpha -1}e^{-\theta x^{\alpha }}\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{\theta+1}}\right) e^{-\theta x^{\alpha }}\right] ^{\beta -1}}{\left( \theta +1\right) \left\{ 1-\left[ 1-\left( 1+{\frac{\theta x^{\alpha }}{1 + \theta}}\right) \ e^{-\theta x^{\alpha }}\right] ^{\beta }\right\} }}}

where \eqn{W_{-1}} denotes the negative branch of the Lambert W function.

\bold{Particular cases:} \eqn{\alpha = 1} the exponentiated Lindley distribution, \eqn{\beta = 1} the power Lindley distribution and \eqn{(\alpha = 1, \beta = 1)} the one-parameter Lindley distribution. See Warahena-Liyanage and Pararai (2014) for other particular cases.
}
\note{
Warahena-Liyanage and Pararai (2014) named the exponentiated power Lindley distribution as generalized power Lindley distribution.
}
\examples{
set.seed(1)
x <- rexpplindley(n = 1000, theta = 11.0, alpha = 5.0, beta = 2.0)
R <- range(x)
S <- seq(from = R[1], to = R[2], by = 0.01)
plot(S, dexpplindley(S, theta = 11.0, alpha = 5.0, beta = 2.0), xlab = 'x', ylab = 'pdf')
hist(x, prob = TRUE, main = '', add = TRUE)

p <- seq(from = 0.1, to = 0.9, by = 0.1)
q <- quantile(x, prob = p)
pexpplindley(q, theta = 11.0, alpha = 5.0, beta = 2.0, lower.tail = TRUE)
pexpplindley(q, theta = 11.0, alpha = 5.0, beta = 2.0, lower.tail = FALSE)
qexpplindley(p, theta = 11.0, alpha = 5.0, beta = 2.0, lower.tail = TRUE)
qexpplindley(p, theta = 11.0, alpha = 5.0, beta = 2.0, lower.tail = FALSE)

## bladder cancer data (from Warahena-Liyanage and Pararai, 2014)
data(bladdercancer)
library(fitdistrplus)
fit <- fitdist(bladdercancer, 'expplindley', start = list(theta = 1, alpha =  1, beta = 1))
plot(fit)


}
\author{
Josmar Mazucheli \email{jmazucheli@gmail.com}

Larissa B. Fernandes \email{lbf.estatistica@gmail.com}
}
\references{
Ashour, S. K., Eltehiwy, M. A., (2015). Exponentiated power Lindley distribution. \emph{Journal of Advanced Research}, \bold{6}, (6), 895-905.

Warahena-Liyanage, G., Pararai, M., (2014). A generalized power Lindley distribution with applications. \emph{Asian Journal of Mathematics and Applications}, \bold{2014}, 1-23.
}
\seealso{
\code{\link[lamW]{lambertWm1}}.
}

