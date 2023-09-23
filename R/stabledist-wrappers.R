#' @name pstable2_and_dstable2
#' @aliases pstable2 dstable2
#' @title PDF and CDF of a skew stable distribution.
#' @description Wrappers for stabledist::pstable and stabledist::dstable
#' @param x Vector of points where the pdf will be evaluated.
#' @param alpha shape / stability parameter, with 0 < alpha <= 2.
#' @param beta skewness parameter, with -1 <= beta <= 1.
#' @param gamma scale parameter, with 0 < sigma.
#' @param delta location parameter, with delta real.
#' For the `2` version, one can have out of bound parameters and
#' an automatic adjustment will be performed to set the parameter
#' to the nearest bound to avoid error message and returning NA.
#' @param pm Parametrization used for the skew stable distribution, as defined by JP Nolan (1997). By default, parametrization = 0.
#' @param tol Relative error tolerance (precision) of the calculated values. By default, tol = 1e-12.
#' @param ... Other stabledist::[dp]stable() arguments.
#' @return A numeric vector.
#' @author Javier Royuela del Val, Federico Simmross Wattenberg and Carlos Alberola López\cr\cr
#'         Maintainer: Javier Royuela del Val <jroyval@@lpi.tel.uva.es>
#' @references Nolan JP (1997). Numerical Calculation of Stable Densities and Distribution Functions. Stochastic Models, 13(4) 759-774.
#' @keywords distribution
#' @importFrom stabledist dstable pstable qstable
#' @export
#' @examples
#' ##take out of bound pars and force them in bounds
#' ## ***internally***
#' pars <- c(2.5, 0.9, 1, 0)
#' x <- seq(-5, 10, 0.001)
#'
#' pdf <- dstable2(x, pars[1], pars[2], pars[3], pars[4])
#' cdf <- pstable2(x, pars[1], pars[2], pars[3], pars[4])
#'
#' ## following would error:
#' ## because pars[1] (alpha > 2)
#' \dontrun{
#' pdf <- dstable(x, pars[1], pars[2], pars[3], pars[4])
#' cdf <- pstable(x, pars[1], pars[2], pars[3], pars[4])
#' }
dstable2 <- function(x, alpha, beta, gamma, delta, pm=0L, tol=1e-12, ...){
  pars <- c(alpha,beta,gamma,delta)
  stabledist::dstable(x,
                      alpha=min(max(pars[1],0+1e-20),2),
                      beta =min(max(pars[2],-1),1),
                      gamma =   min(max(pars[3],0+1e-20),Inf),
                      delta =   pars[4],
                      pm=pm,
                      tol=tol,
                      ...
  )
}
#' @export
pstable2 <- function(x, alpha, beta, gamma, delta, pm=0L, tol=1e-12, ...){
  pars <- c(alpha,beta,gamma,delta)
  stabledist::pstable(x,
                      alpha=min(max(pars[1],0+1e-20),2),
                      beta =min(max(pars[2],-1),1),
                      gamma =   min(max(pars[3],0+1e-20),Inf),
                      delta =   pars[4],
                      pm=pm,
                      tol=tol,
                      ...
  )
}

