##' Vectorized Two-dimensional Romberg Numerical Integration
##'
##' \code{romberg_int2d} performs vectorized numerical integration of a given
##' two-dimensional function.
##'
##'
##' @param f The function (of two variables) to integrate, returning either a
##' scalar or a vector.
##' @param a A two-element vector or a two-column matrix giving the lower
##' bounds. It cannot contain both -Inf and finite values.
##' @param b A two-element vector or a two-column matrix giving the upper
##' bounds. It cannot contain both Inf and finite values.
##' @param eps Precision.
##' @param max The maximum number of steps, by default set to 16.
##' @param d The number of extrapolation points so that 2k is the order of
##' integration, by default set to 5; d=2 is Simpson's rule.
##' @return The vector of values of the integrals of the function supplied.
##' @author Bruce Swihart (based on \code{rmutil::int2} by J.K. Lindsey, who adapted function from Gentleman and Ihaka (2000)
##'    Jr Comp Graph Stat 9, 491-508)
##' @references Gentleman and Ihaka (2000) Jr Comp Graph Stat 9, 491-508
##' @keywords math
##' @examples
##'
##' f <- function(x,y) sin(x)+cos(y)-x^2
##' romberg_int2d(f, a=c(0,1), b=c(2,4))
##' #
##' fn1 <- function(x, y) x^2+y^2
##' fn2 <- function(x, y) (1:4)*x^2+(2:5)*y^2
##' romberg_int2d(fn1, c(1,2), c(2,4))
##' romberg_int2d(fn2, c(1,2), c(2,4))
##' romberg_int2d(fn1, matrix(c(1:4,1:4),ncol=2), matrix(c(2:5,2:5),ncol=2))
##' romberg_int2d(fn2, matrix(c(1:4,1:4),ncol=2), matrix(c(2:5,2:5),ncol=2))
##'
##' @export romberg_int2d
##' @useDynLib gnlrim, .registration = TRUE
###
### vectorized two-dimensional integration
###
romberg_int2d <- function(f, a=c(-Inf,-Inf), b=c(Inf,Inf), eps=1.0e-6, max=16, d=5){
#
# function adapted from Gentleman and Ihaka (2000)
#    Jr Comp Graph Stat 9, 491-508
#
g <- function(y){
	fx <- function(x) f(x,y)
	romberg(fx,a[,2],b[,2])}
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
# function for Romberg integration
#
romberg <- function(f, a=-Inf, b=Inf){
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
#
# check dimensions of limits
#
if(is.vector(a,mode="numeric")&&length(a)==2)a <- matrix(a,ncol=2)
else if(is.matrix(a)){
	if(dim(a)[2]!=2)stop("a must be 2-column matrix")}
else stop("a must be a 2-element vector or a 2-column matrix")
if(is.vector(b,mode="numeric")&&length(b)==2)b <- matrix(b,ncol=2)
else if(is.matrix(b)){
	if(dim(b)[2]!=2)stop("b must be 2-column matrix")}
else stop("b must be a 2-element vector or a 2-column matrix")
if((dim(a)[1]>1&&((any(a[,1]==-Inf)&&all(a[,1]!=-Inf))||(any(a[,2]==-Inf)&&all(a[,1]!=-Inf))))||(dim(b)[1]>1&&((any(b[,1]==Inf)&&all(b[,1]!=Inf))||(any(b[,2]==Inf)&&all(b[,1]!=Inf)))))stop("romberg_int2d cannot have only some limits infinite")
if(length(formals(f))!=2)stop("f must have two arguments")
#
# determine length of vectors to be integrated
#
if(all(a!=-Inf)){
	if(all(b!=Inf)){
		if(any(a[,1]>=b[,1])||any(a[,2]>=b[,2]))stop("some a>=b")
		len <- length(f((a[,1]+b[,1])/2,(a[,2]+b[,2])/2))}
	else len <- length(f(a[,1]+1,a[,2]+1))}
else if(all(b!=Inf))len <- length(f(b[,1]-1,b[,2]-1))
else len <- length(f(0,0))
#
# if a matrix and there are infinite limits, check that all limits are infinite
#
if(any(a!=-Inf)&&dim(a)[1]!=len){
	if(dim(a)[1]!=1)stop("a has incorrect size")
	else a <- matrix(rep(a,len),ncol=2,byrow=TRUE)}
if(any(b!=Inf)&&dim(b)[1]!=len){
	if(dim(b)[1]!=1)stop("b has incorrect size")
	else b <- matrix(rep(b,len),ncol=2,byrow=TRUE)}
#
# integrate
#
romberg(g,a[,1],b[,1])}
