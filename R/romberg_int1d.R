##' Vectorized One-dimensional Romberg Numerical Integration
##'
##' \code{romberg_int1d} performs numerical romberg integration of 1-dimensional
##' functions. Vectorized to handle several functions at once. See examples.
##'
##'
##' @param f The function (of one variable) to integrate, returning either a
##' scalar or a vector.
##' @param a A scalar or vector (only Romberg) giving the lower bound(s). A
##' vector cannot contain both -Inf and finite values.
##' @param b A scalar or vector (only Romberg) giving the upper bound(s). A
##' vector cannot contain both Inf and finite values.
##' @param eps Precision.
##' @param max The maximum number of steps, by default set to 16.
##' @param d The number of extrapolation points so that 2d is the
##' order of integration, by default set to 5; d=2 is Simpson's rule.
##' @return The vector of values of the integrals of the function supplied.
##' @author Bruce Swihart (based on \code{rmutil::int} by J.K. Lindsey.)
##' @keywords math
##' @examples
##'
##' f <- function(x) sin(x)+cos(x)-x^2
##' romberg_int1d(f, a=0, b=2)
##'
##' #
##' f <- function(x) exp(-(x-2)^2/2)/sqrt(2*pi)
##' romberg_int1d(f, a=0:3)
##' romberg_int1d(f, a=0:3, d=2)
##' 1-pnorm(0:3, 2)
##' #
##' f <- function(x) dnorm(x)
##' romberg_int1d(f, a=-Inf, b=qnorm(0.975))
##'
##' # NOTE: see \code{rmutil::int} for TOMS614.
##' @export romberg_int1d
##' @useDynLib gnlrim, .registration = TRUE
###
### vectorized one-dimensional integration
###
romberg_int1d <- function(f, a=-Inf, b=Inf, eps=0.0001,
                max=NULL, d=NULL){
  type = "Romberg"
  #
  # function to call the C code
  #
  int1 <- function(ff, aa, bb){
    envir2 <- environment(fun=ff)
    z <- .Call("romberg_sexp",
               ff,
               as.double(aa),
               as.double(bb),
               len=as.integer(len),
               eps=as.double(eps),
               pts=as.integer(d),
               max=as.integer(max),
               err=integer(1),
               envir2,
               PACKAGE="gnlrim")
    z
  }
  #
  # check algorithm to be used and initialize parameters
  #
  type <- match.arg(type,c("Romberg","TOMS614"))
  if(is.null(max))max <- if(type=="Romberg") 16 else 100
  if(is.null(d))d <- if(type=="Romberg") 5 else 1
  #
  # check function and integration limits
  #
  if(length(formals(f))!=1)stop("f must have one argument")
  if(any(a==Inf))stop("lower bound cannot be Inf")
  if(any(b==-Inf))stop("upper bound cannot be -Inf")
  if((length(a)>1&&any(a==-Inf)&&all(a!=-Inf))||(length(b)>1&&any(b==Inf)&&all(b!=Inf)))stop("int cannot be vectorized with some limits infinite")
  #
  # determine length of vector to be integrated
  #
  if(all(a!=-Inf)){
    if(all(b!=Inf)){
      if(any(a>=b))stop("some a>=b")
      len <- length(f((a+b)/2))}
    else len <- length(f(a+1))}
  else if(all(b!=Inf))len <- length(f(b-1))
  else len <- length(f(0))
  if(len>1&&type!="Romberg")stop("vector functions only allowed with Romberg")
  #
  # if a vector and there are infinite limits, check that all limits are infinite
  #
  if(all(a!=-Inf)&&length(a)!=len){
    if(length(a)!=1)stop("a has incorrect length")
    else a <- rep(a,len)}
  if(all(b!=Inf)&&length(b)!=len){
    if(length(b)!=1)stop("b has incorrect length")
    else b <- rep(b,len)}
  if(type=="Romberg"){
    # invert function for infinite limits
    ff <- function(x) f(1/x)/(x*x)
    if(all(b==Inf)){
      if(all(a==-Inf))
        # both limits infinite
        z <- int1(ff,rep(-1,len),rep(0,len))+
          int1(f,rep(-1,len),rep(1,len))+
          int1(ff,rep(0,len),rep(1,len))
      else {
        # only upper limit infinite, cut in 2 pieces about 0
        if(any(a>0)){
          if(any(a<=0))a1 <- ifelse(a>0,a,1)
          else a1 <- a
          z1 <- int1(ff,rep(0,len), 1/a1)}
        else z1 <- rep(0,len)
        if(any(a<=0)){
          if(any(a>0))a1 <- ifelse(a<=0,a,0)
          else a1 <- a
          z2 <- int1(f,a1,rep(1,len))+
            int1(ff,rep(0,len),rep(1,len))}
        else z2 <- rep(0,len)
        z <- z1*(a>0)+z2*(a<=0)}}
    else if(all(a==-Inf)){
      # only lower limit infinite, cut in 2 pieces about 0
      if(any(b<0)){
        if(any(b>=0))b1 <- ifelse(b<0,b,1)
        else b1 <- b
        z1 <- int1(ff, 1/b1,rep(0,len))}
      else z1 <- rep(0,len)
      if(any(b>=0)){
        if(any(b<0))b1 <- ifelse(b>=0,b,0)
        else b1 <- b
        z2 <- int1(f,rep(-1,len), b1)+
          int1(ff,rep(-1,len),rep(0,len))}
      else z2 <- rep(0,len)
      z <- z1*(b<0)+z2*(b>=0)}
    else z <- int1(f, a, b)
    z}
  else {
    #
    # TOMS614
    #
  }
}
