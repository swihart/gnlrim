#' @name stable_pdf_and_cdf
#' @aliases stable_cdf2 stable_pdf2
#' @title PDF and CDF of a skew stable distribution.
#' @description Evaluate the PDF or the CDF of the skew stable distribution with parameters
#' pars = c(alpha, beta, sigma, mu) at the points given in x.\cr\cr
#' _parametrization_ argument specifies the parametrization used for the distribution
#' as described by JP Nolan (1997). The default value is _parametrization_ = 0.\cr\cr
#' _tol_ sets the relative error tolerance (precision) to _tol_. The default value is tol = 1e-12.
#' @param x Vector of points where the pdf will be evaluated.
#' @param pars Vector with an initial estimation of the parameters. `pars_init = c(alpha, beta, sigma, mu)`, where
#' * alpha: shape / stability parameter, with 0 < alpha <= 2.
#' * beta: skewness parameter, with -1 <= beta <= 1.
#' * sigma: scale parameter, with 0 < sigma.
#' * mu: location parameter, with mu real.
#' For the `2` version, one can have out of bound parameters and
#' an automatic adjustment will be performed to set the parameter
#' to the nearest bound to avoid error message and returning NA.
#' @param parametrization Parametrization used for the skew stable distribution, as defined by JP Nolan (1997). By default, parametrization = 0.
#' @param tol Relative error tolerance (precision) of the calculated values. By default, tol = 1e-12.
#' @return A numeric vector.
#' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López\cr\cr
#'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
#' @references Nolan JP (1997). Numerical Calculation of Stable Densities and Distribution Functions. Stochastic Models, 13(4) 759-774.
#' @keywords distribution
#' @importFrom libstableR stable_pdf stable_cdf stable_q
#' @export
#' @examples
#' ##take out of bound pars and force them in bounds
#' ## ***internally***
#' pars <- c(2.5, 0.9, 1, 0)
#' x <- seq(-5, 10, 0.001)
#'
#' pdf <- stable_pdf2(x, pars)
#' cdf <- stable_cdf2(x, pars)
#'
#' ## following would error:
#' ## because pars[1] (alpha > 2)
#' \dontrun{
#' pdf <- stable_pdf2(x, pars)
#' cdf <- stable_cdf2(x, pars)
#' }
stable_pdf2 <- function(x, pars, parametrization=0L, tol=1e-12){
  libstableR::stable_pdf(x,
                         c(min(max(pars[1],0+1e-20),2),
                           min(max(pars[2],-1),1),
                           min(max(pars[3],0+1e-20),Inf),
                           pars[4]),
                         parametrization,
                         tol
  )
}
#' @export
stable_cdf2 <- function(x, pars, parametrization=0L, tol=1e-12){
  libstableR::stable_cdf(x,
                         c(min(max(pars[1],0+1e-20),2),
                           min(max(pars[2],-1),1),
                           min(max(pars[3],0+1e-20),Inf),
                           pars[4]),
                         parametrization,
                         tol
  )
}

