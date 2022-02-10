########################################$$$$#########
## 2022-02-10                                      ##
## START two random parameters                     ##
#####################################################
# Derived expression for gamma
g <- function(a) {
  iu <- complex(real=0, imaginary=1)
  return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}

bridgecloglog_rstable <- function(n, alpha, delta){
  mult <- (delta/alpha)^(1/alpha)
  X <- stabledist::rstable(n , alpha, beta=1, gamma=g(alpha), delta=0, pm=1)
  Z <- log(mult * X)
  Z
}

## add in beta random effect
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Logistic-BIVARIATE_NORM", alpha=1.89, gamma=1.2, delta=1){

  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }

  n <- n1 + n2
  u <- round(apply(G(n,v1,v2,rho),2, function(W) rep(W,each=J)),2)

  ## x <- c(rep(1, n1*J), rep(0, n2*J))
  ##  x <-c(runif(n1*J, 0.5,1.5), runif(n2*J, 0,0.60))
  x <-c(sample(c(1,2,3,4)/10,n1*J,TRUE), sample(c(1.2,2.2,3.2,4.2)/10,n2*J,TRUE))
  eta <- round(a0 + a1*x,2)

  eta_i <- round( (a0 + u[,1]) + (a1+u[,2])*x, 2)
  py1 <- round(H(eta_i),2)
  y <- rbinom(length(eta_i), 1, prob=py1 )

  data.frame(id=rep(1:n, each=J),
             j = rep(1:J),
             x1 = x,
             eta = eta,
             u_i = u,
             eta_i = eta_i,
             py1 = py1,
             y=y
  )

}


detach(summed_binom_dat)
set.seed(5)
binom_dat <-
#  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(2,2, J=100, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat
## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial, nAGQ = 0)



attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

(rand.int.rand.slopes.nonzero.corr <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + rand2*b_p),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=3, var2=3, corr12= 0.20),
           p_uppb = c(  0,   2, 4.00, 4.00, 0.90),
           p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
           distribution="binomial",
           nest=id,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-corr",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb'
    )
)

########################################$$$$#########
## 2022-02-10                                      ##
## END two random parameters                       ##
#####################################################

########################################$$$$#########
## 2022-02-09                                      ##
## START two random parameters                     ##
#####################################################
library(lme4)
# from ?glmer
## generalized linear mixed model
library(lattice)
xyplot(incidence/size ~ period|herd, cbpp, type=c('g','p','l'),
       layout=c(3,5), index.cond = function(x,y)max(y))
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
## using nAGQ=0 only gets close to the optimum
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 0))
## using  nAGQ = 9  provides a better evaluation of the deviance
## Currently the internal calculations use the sum of deviance residuals,
## which is not directly comparable with the nAGQ=0 or nAGQ=1 result.
## 'verbose = 1' monitors iteratin a bit; (verbose = 2 does more):
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, verbose = 1, nAGQ = 9))

## GLMM with individual-level variability (accounting for overdispersion)
## For this data set the model is the same as one allowing for a period:herd
## interaction, which the plot indicates could be needed.
cbpp$obs <- 1:nrow(cbpp)
(gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
                (1 | herd) +  (1|obs),
              family = binomial, data = cbpp))
anova(gm1,gm2)

# let's focus on gm2 since it has 2 random parameters:
summary(gm2)


## But first can I reproduce gm1?
cbpp$period2 <- (cbpp$period == 2) + 0L
cbpp$period3 <- (cbpp$period == 3) + 0L
cbpp$period4 <- (cbpp$period == 4) + 0L
cbpp$period_numeric <- as.numeric(cbpp$period)

attach(cbpp)
gm1.redux <-
  gnlrem(y=cbind(incidence, size - incidence),
         mu = ~ plogis(Intercept + period2*b2 + period3*b3 + period4*b4 + rand),
         pmu = c(Intercept=-1,b2=-1, b3=-1, b4=-1),
         pmix=c(var=0.78),
         p_uppb = c(  2,   2, 2, 2,  1.0),
         p_lowb = c( -2,  -2,-2,-2,  0.1),
         distribution="binomial",
         nest=herd,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         abs.tol.nlminb = 1e-6,
         xf.tol.nlminb =  1e-6,
         x.tol.nlminb =   1e-6,
         rel.tol.nlminb = 1e-6
  )

gm1.redux$value
summary(gm1)$logLik

summary(gm1a)$loglik

gm1a.redux <-
  gnlrem(y=cbind(incidence, size - incidence),
         mu = ~ plogis(Intercept + period2*b2 + period3*b3 + period4*b4 + rand),
         pmu = c(Intercept=-1,b2=-1, b3=-1, b4=-1),
         pmix=c(var=0.78),
         p_uppb = c(  2,   2, 2, 2,  1.0),
         p_lowb = c( -2,  -2,-2,-2,  0.1),
         distribution="binomial",
         nest=herd,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         abs.tol.nlminb = 1e-7,
         xf.tol.nlminb =  1e-7,
         x.tol.nlminb =   1e-7,
         rel.tol.nlminb = 1e-7,
         points = 10, #5 default
         steps = 20   # 10 is default
  )

 gm1.redux$value
gm1a.redux$value

## tips on covariance...
## https://stats.stackexchange.com/questions/414864/lme4glmer-get-the-covariance-matrix-of-the-fixed-and-random-effect-estimates



(gm2.redux <-
  gnlrem(y=cbind(incidence, size - incidence),
         mu = ~ plogis(Intercept + period2*b2 + period3*b3 + period4*b4 + rand1 + rand2),
         pmu = c(Intercept=-1,b2=-1, b3=-1, b4=-1),
         pmix=c(var1=0.78, var2=0.40),
         p_uppb = c(  2,   2, 2, 2,  1.0, 1.0),
         p_lowb = c( -2,  -2,-2,-2,  0.1, 0.1),
         distribution="binomial",
         nest=herd,
         random=c("rand1", "rand2"),
         mixture="bivariate-normal-indep",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
  )
)


## glmer with correlation between random intercept and random slope
glmer(cbind(incidence, size - incidence) ~ period_numeric + (period_numeric | herd), cbpp, binomial, nAGQ = 0)
# Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite
#                                                           Quadrature, nAGQ = 0) [glmerMod]
# Family: binomial  ( logit )
# Formula: cbind(incidence, size - incidence) ~ period_numeric + (period_numeric |      herd)
# Data: cbpp
# AIC      BIC   logLik deviance df.resid
# 194.3058 204.4325 -92.1529 184.3058       51
# Random effects:
#   Groups Name           Std.Dev. Corr
# herd   (Intercept)      1.0603
# period_numeric          0.3311   -0.86
# Number of obs: 56, groups:  herd, 15
# Fixed Effects:
#   (Intercept)  period_numeric
# -0.8765         -0.5843
##
##
#
#
## glmer wout correlation between random intercept and random slope
glmer(cbind(incidence, size - incidence) ~ period_numeric + (1 | herd) + (0+period_numeric|herd), cbpp, binomial, nAGQ = 0)
# Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite
#                                                           Quadrature, nAGQ = 0) [glmerMod]
# Family: binomial  ( logit )
# Formula: cbind(incidence, size - incidence) ~ period_numeric + (1 | herd) +
#   (0 + period_numeric | herd)
# Data: cbpp
# AIC      BIC   logLik deviance df.resid
# 194.6247 202.7261 -93.3124 186.6247       52
# Random effects:
#   Groups Name           Std.Dev.
# herd   (Intercept)    0.6634
# herd.1 period_numeric 0.0000
# Number of obs: 56, groups:  herd, 15
# Fixed Effects:
#   (Intercept)  period_numeric
# -0.8987         -0.5517



(rand.int.rand.slopes <-
    gnlrem(y=cbind(incidence, size - incidence),
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + rand2*b_p),
           pmu = c(Intercept=-1, b_p=-1),
           pmix=c(var1=0.78, var2=0.40),
           p_uppb = c(  2,   2, 4.00, 4.00),
           p_lowb = c( -2,  -2, 0.1, 0.1),
           distribution="binomial",
           nest=herd,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-indep",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb'
    )
)
# [1] 4
# Intercept       b_p      var1      var2
# -1.00     -1.00      0.78      0.40
# [1] 103.7967
# fn is  fn
# Looking for method =  nlminb
# Function has  4  arguments
# par[ 1 ]:  -2   <? -1   <? 2     In Bounds
# par[ 2 ]:  -2   <? -1   <? 2     In Bounds
# par[ 3 ]:  0.1   <? 0.78   <? 4     In Bounds
# par[ 4 ]:  0.1   <? 0.4   <? 4     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.39794   log bounds ratio= 0.01099538
# Method:  nlminb
# 0:     103.79673: -1.00000 -1.00000 0.780000 0.400000
# 1:     94.395935: -0.921281 -0.625287 0.773239 0.393239
# 2:     93.454422: -0.996482 -0.490285 0.443614 0.274224
# 3:     93.305688: -1.00875 -0.550801 0.437306 0.272707
# 4:     93.248260: -0.959084 -0.551177 0.401156 0.263686
# 5:     93.242240: -0.948738 -0.547427 0.374704 0.255650
# 6:     93.239082: -0.932206 -0.564129 0.366603 0.252365
# 7:     93.237854: -0.925698 -0.560545 0.368247 0.252491
# 8:     93.237583: -0.932887 -0.558137 0.367660 0.252801
# 9:     93.237511: -0.933923 -0.558984 0.369611 0.253577
# 10:     93.237500: -0.933290 -0.559070 0.368478 0.253123
# 11:     93.237499: -0.933371 -0.558983 0.368667 0.253196
# 12:     93.237499: -0.933370 -0.558997 0.368656 0.253192
# 13:     93.237499: -0.933369 -0.558997 0.368655 0.253192
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2
# -0.9333689 -0.5589967  0.3686549  0.2531917
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 93.2375
#
# $fevals
# function
# 19
#
# $gevals
# gradient
# 74
#
# $nitns
# [1] 13
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 1391.858
#
# Assemble the answers
# Intercept        b_p      var1      var2   value fevals gevals niter convcode kkt1
# nlminb -0.9333689 -0.5589967 0.3686549 0.2531917 93.2375     19     74    13        0   NA
# kkt2    xtime
# nlminb   NA 1391.858

(rand.int.rand.slopes.nonzero.corr <-
    gnlrem(y=cbind(incidence, size - incidence),
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + rand2*b_p),
           pmu = c(Intercept=-0.9333689, b_p=-0.558996),
           pmix=c(var1=0.3686549, var2=0.2531917, corr12= -0.10),
           p_uppb = c(  2,   2, 4.00, 4.00, 0.90),
           p_lowb = c( -2,  -2, 0.05, 0.04,-0.90),
           distribution="binomial",
           nest=herd,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-corr",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb'
    )
)
#
# fn is  fn
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -2   <? -0.9333689   <? 2     In Bounds
# par[ 2 ]:  -2   <? -0.558996   <? 2     In Bounds
# par[ 3 ]:  0.05   <? 0.3686549   <? 4     In Bounds
# par[ 4 ]:  0.04   <? 0.2531917   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? -0.1   <? 0.9     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.9700533   log bounds ratio= 0.3467875
# Method:  nlminb
# 0:     93.246926: -0.933369 -0.558996 0.368655 0.253192 -0.100000
# 1:     93.240001: -0.935611 -0.553697 0.352205 0.247219 -0.0946299
# 2:     93.238462: -0.937009 -0.558041 0.350772 0.246712 -0.0941815
# 3:     93.237903: -0.936588 -0.556481 0.346662 0.245240 -0.0928892
# 4:     93.237721: -0.936962 -0.558312 0.345884 0.244963 -0.0926473
# 5:     93.237611: -0.935713 -0.557602 0.344551 0.244488 -0.0922390
# 6:     93.237524: -0.932448 -0.559458 0.343021 0.243941 -0.0917914
# 7:     93.237500: -0.933385 -0.559064 0.342267 0.243669 -0.0915783
# 8:     93.237500: -0.933358 -0.558945 0.342223 0.243654 -0.0915649
# 9:     93.237499: -0.933394 -0.558982 0.342113 0.243615 -0.0915364
# 10:     93.237499: -0.933360 -0.559000 0.342104 0.243612 -0.0915342
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept         b_p        var1        var2      corr12
# -0.93336043 -0.55900009  0.34210449  0.24361159 -0.09153423
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 93.2375
#
# $fevals
# function
# 17
#
# $gevals
# gradient
# 73
#
# $nitns
# [1] 10
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 1167.854
#
# Assemble the answers
# Intercept        b_p      var1      var2      corr12   value fevals gevals niter
# nlminb -0.9333604 -0.5590001 0.3421045 0.2436116 -0.09153423 93.2375     17     73    10
# convcode kkt1 kkt2    xtime
# nlminb        0   NA   NA 1167.854

## BELOW HERE: tried to dichotomize lmer example for my glmer purposes
# (fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
# summary(fm1)
#
# (fm2 <- lmer(Reaction ~ Days + (1 | Subject) + (0+Days|Subject), sleepstudy))
# summary(fm2)
#
#
# (fm1 <- glmer(I(Reaction>400) ~ Days + (Days | Subject), data=sleepstudy, binomial, nAGQ=0))
# summary(fm1)
#
# (fm2 <- glmer(I(Reaction>200) ~ Days + (1 | Subject) + (0+Days|Subject), data=sleepstudy, binomial, nAGQ=0))
# summary(fm2)



########################################$$$$#########
## 2022-02-09                                      ##
## END two random parameters                       ##
#####################################################

########################################$$$$#########
## 2021-12-14                                      ##
## START: 0 mean intercept - cloglog bridge        ##
#####################################################
# Derived expression for gamma
g <- function(a) {
  iu <- complex(real=0, imaginary=1)
  return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}

bridgecloglog_rstable <- function(n, alpha, delta){
  mult <- (delta/alpha)^(1/alpha)
  X <- stabledist::rstable(n , alpha, beta=1, gamma=g(alpha), delta=0, pm=1)
  Z <- log(mult * X)
  Z
}

## add in beta random effect
sim_mrim_data <- function(n1, n2, J, a0, a1, v=1, mrim="CLOGLOG-CLOGLOG-LPS", alpha=1.89, gamma=1.2, delta=1){

  if(mrim=="PPN"){
    G <- function(n,v){rnorm(n, s=sqrt(v))}
    H <- function(x) pnorm(x)
  }

  if(mrim=="LLB"){
    G <- function(n,v){bridgedist::rbridge(n, scale=1/sqrt(1+3/pi^2*v))}
    H <- function(x) plogis(x)
  }

  if(mrim=="SSS"){
    G <- function(n,v){rstable(n, alpha, 0, v, 0, 0)}
    H <- function(x) pstable(x, alpha, 0, gamma, 0, 0)
  }

  if(mrim=="CLOGLOG-CLOGLOG-LPS" & alpha!=delta){
    #library(evd)
    #source("functions.R")
    G <- function(n, v){
      ## this combo should give variance 1 mean 0
      (aa <- 1/sqrt(1+6*pi^-2*v))
      (dd <- aa*exp(-digamma(1)*(aa-1)))
      print("IN HERE")
      bridgecloglog_rstable(n, aa, dd)
    }
    H <- function(x){
      #1-pgumbel(x, loc=0, scale=1)
      1-exp(-exp(x))
    }
  }
  if(mrim=="CLOGLOG-CLOGLOG-LPS" & alpha==delta){
    #library(evd)
    #source("functions.R")
    G <- function(n, v){
      ## this combo should give variance 1 mean 0
      (aa <- 1/sqrt(1+6*pi^-2*v))
      (dd <- aa)
      bridgecloglog_rstable(n, aa, dd)
    }
    H <- function(x){
      #1-pgumbel(x, loc=0, scale=1)
      1-exp(-exp(x))
    }
  }
  if(mrim=="binomial-beta-HGLM"){
    G <- function(n, v){rbeta(n, v, v)}
    H <- function(x) plogis(x)
  }

  n <- n1 + n2
  u <- round(rep(G(n,v), each=J),2)

  x <- c(rep(1, n1*J), rep(0, n2*J))

  eta <- round(a0 + a1*x,2)

  eta_i <- round(eta + u,2)
  py1 <- round(H(eta_i),2)
  y <- rbinom(length(eta_i), 1, prob=py1 )

  data.frame(id=rep(1:n, each=J),
             j = rep(1:J),
             group = x,
             eta = eta,
             u_i = u,
             eta_i = eta_i,
             py1 = py1,
             y=y
  )

}


detach(summed_binom_dat)
set.seed(5)
binom_dat <-
  sim_mrim_data(500,120, J=47, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","group")]

attach(summed_binom_dat)

ybind <- cbind(r,n_r)
x1 <- group

gnlrim.a.eq.d <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(a_eq_d=0.78),
         p_uppb = c(  10,   2,   1.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="cloglog-bridge-delta-eq-alpha",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         abs.tol.nlminb = 1e-6,
         xf.tol.nlminb =  1e-6,
         x.tol.nlminb =   1e-6,
         rel.tol.nlminb = 1e-6
  )

## I want to demonstrate that when
## a != d we need to add an intercept adjustment
## for marginalization of coeffs.
## we used the fit above to 1) validate simulated values
## and 2) know where to aim for marg coeff bc when a==d
## we have straightforward multiply b0m = a*b0c
gnlrim.mean.0 <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=-2,treatment=1),
         pmix=c(var=2),
         p_uppb = c(  10,   2,   100),
         p_lowb = c( -10,  -2,   0.01),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         abs.tol.nlminb = 1e-6,
         xf.tol.nlminb =  1e-6,
         x.tol.nlminb =   1e-6,
         rel.tol.nlminb = 1e-6
  )


rbind(gnlrim.a.eq.d[1:2],
      gnlrim.mean.0[1:2])


calc.a <- 1/sqrt(1+6*pi^-2*gnlrim.mean.0[3])
calc.d <- calc.a*exp(-digamma(1)*(calc.a-1))


## looks like these two are identical:
log(calc.d/calc.a)/calc.a
gnlrim.a.eq.d[1] - gnlrim.mean.0[1]

## so the implication is that we need to add
## an adjustment when we use 0-mean param:
gnlrim.mean.0[1]+log(calc.d/calc.a)/calc.a
gnlrim.a.eq.d[1]

## looks like these two are identical:
log(calc.d/calc.a)
(gnlrim.a.eq.d[1]*gnlrim.a.eq.d[3]) - (gnlrim.mean.0[1]*gnlrim.a.eq.d[3])



########################################$$$$#########
## 2021-12-14                                      ##
## END:   0 mean intercept - cloglog bridge        ##
#####################################################

########################################$$$$#########
## 2021-12-04                                      ##
## START: ZINBRI TRA/TBA  STOP THE COUNT example   ##
#####################################################
library(data.table)
library(glmmTMB)
## rewrite vectorize -- originally from macbook, readme_toc.rmd 12/23/2015:
## beta= -2.303 should give 90% TRA 100*(1-exp(-2.303))=90.00
createData <- function(nc=NUM,
                       nt=NUM,
                       mu=exp(2.57),
                       beta=0,#-2.303,##tra2beta(TRA),
                       sigf=1.042,
                       sigc=SIGC,
                       pi= 0.056,
                       theta=1.93,
                       testconc=1,
                       sampleName="s1",
                       feedNumber=1,
                       num_pints=1){## num_pints PER treatment and control

  ## feed effect
  alphai<-rep(rnorm(length(feedNumber),mean=0,sd=sigf), each=2)

  muci <- mu*exp(alphai)

  ## pint effect
  deltaij<-rnorm(2*length(feedNumber),mean=0,sd=sigc)

  mucij <-muci * exp(deltaij)

  ## make every second one the treatment; multiply beta.
  mucij <- mucij * rep(c(1,exp(beta)), length(feedNumber))


  ## negative binomial draws, nc in each com, treated following control
  y <- rnbinom(nc*length(feedNumber)*2, mu = rep(mucij, each=nc), size=theta)

  ## zero process
  z <- rbinom(length(y),1,1-pi)

  ## zero-inflated neg bin:
  y[z==0] <- 0


  ## trmt
  trmt<-c(rep(0,nc),rep(       1,nt))

  ## feed
  feed <- rep(feedNumber, each=2*nc)

  ## combine into data.table
  data.table(feed, trmt, y, key=c("feed","trmt","y"))

}

## pi=0, theta=1 should give us logistic
## beta= -2.303 should give 90% TRA 100*(1-exp(-2.303))=90.00
set.seed(88)
zinbri_dat <- copy(createData(nc=20, nt=20, mu=exp(2.57), beta= -2.303, sigf=1, sigc=0, pi=0.05, theta=1.93, feedNumber = 1:5))
zinbri_dat[, mean(y), by=c("feed", "trmt")]
zinbri_dat[, var(y), by=c("feed", "trmt")]
zinbri_dat[,{j=list(mean=mean(y),var=var(y))}, by=c("feed","trmt")]

sum(zinbri_dat$y)
## assume the lab techs counted all 21,134
## then you could fit this model:
zinbri <- glmmTMB(y~trmt +(1|feed) ,
                  ziformula=~0,
                  family=nbinom2(),
                  data=zinbri_dat
)
summary(zinbri)

## or, assume they just identified 0 or at least 1.
## so here they checked all 2000 mosquitoes but only counted
sum(zinbri_dat$y>0)
## 1476 eggs.
## you know the RCM theta and pi.  So...
detach(summed_binom_zinbri_dat)
summed_binom_zinbri_dat <-
  zinbri_dat[, {j=list(r=sum(y>0), n_r=sum(y==0))}, by=c("feed","trmt")]


attach(summed_binom_zinbri_dat)

ybind <- cbind(r,n_r)
x1 <- trmt

sim_binom_locked <-
  gnlrim(y=ybind,
         mu = ~ (1-piz)*(1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta),
         pmu = c(piz=0.05, Intercept=0,beta_trmt=0,theta=1.93),
         pmix=c(var=20),
         p_uppb = c(0.05,  10,   10, 1.93 ,200.0),
         p_lowb = c(0.05, -10,  -10, 1.93,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'#,
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5


  )


sim_binom_free <-
  gnlrim(y=ybind,
         mu = ~ (1-piz)*(1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta),
         pmu = c(piz=0.05, Intercept=0,beta_trmt=0,theta=1.93),
         pmix=c(var=20),
         p_uppb = c(0.99,  10,   10, 100,  200.0),
         p_lowb = c(0.01, -10,  -10, 0.01,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'#,
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5


  )

summary(zinbri)
rbind(
sim_binom_locked,
sim_binom_free)


########################################$$$$#########
## 2021-12-04                                      ##
## END:   ZINBRI TRA/TBA  STOP THE COUNT example   ##
#####################################################

########################################
## 2021-12-02                         ##
## START: ZINBRI TRA/TBA   relation   ##
########################################
library(data.table)
library(glmmTMB)
## rewrite vectorize -- originally from macbook, readme_toc.rmd 12/23/2015:
createData <- function(nc=NUM,
                       nt=NUM,
                       mu=MUC,
                       beta=BETA,##tra2beta(TRA),
                       sigf=SIGF,
                       sigc=SIGC,
                       pi= PI,
                       theta=THETA,
                       testconc=1,
                       sampleName="s1",
                       feedNumber=1,
                       num_pints=1){## num_pints PER treatment and control

  ## feed effect
  alphai<-rep(rnorm(length(feedNumber),mean=0,sd=sigf), each=2)

  muci <- mu*exp(alphai)

  ## pint effect
  deltaij<-rnorm(2*length(feedNumber),mean=0,sd=sigc)

  mucij <-muci * exp(deltaij)

  ## make every second one the treatment; multiply beta.
  mucij <- mucij * rep(c(1,exp(beta)), length(feedNumber))


  ## negative binomial draws, nc in each com, treated following control
  y <- rnbinom(nc*length(feedNumber)*2, mu = rep(mucij, each=nc), size=theta)

  ## zero process
  z <- rbinom(length(y),1,1-pi)

  ## zero-inflated neg bin:
  y[z==0] <- 0


  ## trmt
  trmt<-c(rep(0,nc),rep(       1,nt))

  ## feed
  feed <- rep(feedNumber, each=2*nc)

  ## combine into data.table
  data.table(feed, trmt, y, key=c("feed","trmt","y"))

}

## pi=0, theta=1 should give us logistic
## beta= -2.303 should give 90% TRA 100*(1-exp(-2.303))=90.00
set.seed(88)
zinbri_dat <- createData(nc=2000, nt=2000, mu=10, beta= -2.303, sigf=1, sigc=0, pi=0, theta=1, feedNumber = 1:50)
zinbri_dat[, mean(y), by=c("feed", "trmt")]
zinbri_dat[, var(y), by=c("feed", "trmt")]
zinbri_dat[,{j=list(mean=mean(y),var=var(y))}, by=c("feed","trmt")]

zinbri <- glmmTMB(y~trmt +(1|feed) ,
        ziformula=~0,
        family=nbinom2(),
        data=zinbri_dat
)
summary(zinbri)

detach(summed_binom_zinbri_dat)
summed_binom_zinbri_dat <-
  zinbri_dat[, {j=list(r=sum(y>0), n_r=sum(y==0))}, by=c("feed","trmt")]

summed_binom_zinbri_dat[,{j=list(mean=mean(y),var=var(y))}, by=c("feed","trmt")]
attach(summed_binom_zinbri_dat)

ybind <- cbind(r,n_r)
x1 <- trmt

sim_binom_zinbri_gnlrim <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*beta_trmt + rand),
         pmu = c(Intercept=0,beta_trmt=0),
         pmix=c(var=2),
         p_uppb = c(  10,   10, 200.0),
         p_lowb = c( -10,  -10,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )


## NOTE:  this was erroring because theta was coming first...
# sim_binom_lock1theta_error <-
#   gnlrim(y=ybind,
#          mu = ~ 1-(theta/( exp(Intercept + x1*beta_trmt + rand) + theta ) )^theta,
#          pmu = c(Intercept=0,beta_trmt=0,theta=1),
#          pmix=c(var=2),
#          p_uppb = c(  10,   10, 1 ,200.0),
#          p_lowb = c( -10,  -10, 1,   0.1),
#          distribution="binomial",
#          nest=feed,
#          random="rand",
#          mixture="normal-var",
#          ooo=TRUE,
#          compute_hessian = FALSE,
#          compute_kkt = FALSE,
#          trace=1,
#          method='nlminb'
#   )


##
sim_binom_lock1theta <-
  gnlrim(y=ybind,
         mu = ~ 1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta,
         pmu = c(Intercept=0,beta_trmt=0,theta=1),
         pmix=c(var=2),
         p_uppb = c(  10,   10, 1 ,200.0),
         p_lowb = c( -10,  -10, 1,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## samesies!
sim_binom_zinbri_gnlrim
sim_binom_lock1theta



## test beta-binom-TBA by locking theta=1 pi =0
sim_betabinomTBA <-
  gnlrim(y=ybind,

         mu = ~ (1-piz)*(1-(theta/( exp(Intercept + x1*beta_trmt + rand) + theta ) )^theta),
         pmu = c(piz=0, theta=1, Intercept=0, beta_trmt=0),
         pmix=c(var=1),
         p_uppb = c(0, 1, 10, 10, 20),
         p_lowb = c(0, 1,-10,-10,  0.1),

         distribution="beta-binomial-TBA",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## samesies!
sim_binom_zinbri_gnlrim
sim_binom_lock1theta
sim_betabinomTBA

## relax theta
sim_binom_freerangetheta <-
  gnlrim(y=ybind,
         mu = ~ 1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta,
         pmu = c(Intercept=0,beta_trmt=0,theta=1),
         pmix=c(var=2),
         p_uppb = c(  10,   10, 10 ,200.0),
         p_lowb = c( -10,  -10, 0.1,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )




## test beta-binom-TBA by locking theta=1 pi =0
sim_betabinomTBA_free_range <-
  gnlrim(y=ybind,

         mu = ~ (1-piz)*(1-(theta/( exp(Intercept + x1*beta_trmt + rand) + theta ) )^theta),
         pmu = c(piz=0, theta=1, Intercept=0, beta_trmt=0),
         pmix=c(var=2),
         p_uppb = c(0, 1  , 10, 10, 50),
         p_lowb = c(0, 1,-10,-10, 0.1),

         distribution="beta-binomial-TBA",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_binom_freerangetheta
sim_betabinomTBA_free_range






##
###
#### Let's generate data that is theta=0.43
##
##
#

## pi=0, theta=1 should give us logistic
set.seed(1110)
zinbri_dat <- createData(nc=2000, nt=2000, mu=10, beta= -2.303, sigf=1, sigc=0, pi=0, theta=0.43, feedNumber = 1:50)
zinbri_dat[, mean(y), by=c("feed", "trmt")]



detach(summed_binom_zinbri_dat)
summed_binom_zinbri_dat <-
  zinbri_dat[, {j=list(r=sum(y>0), n_r=sum(y==0))}, by=c("feed","trmt")]

attach(summed_binom_zinbri_dat)

ybind <- cbind(r,n_r)
x1 <- trmt

## relax theta
sim_binom_freerangetheta <-
  gnlrim(y=ybind,
         mu = ~ 1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta,
         pmu = c(Intercept=0,beta_trmt=0,theta=0.432),
         pmix=c(var=20),
         p_uppb = c(  10,   10, 0.432 ,200.0),
         p_lowb = c( -10,  -10, 0.432,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'#,
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5


  )
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept beta_trmt     theta       var
# 2.532540 -2.307376  0.432000  1.030043
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 577.3475



sim_betabinomTBA_free_range <-
  gnlrim(y=ybind,

         mu = ~ (1-piz)*(1-(theta/( exp(Intercept + x1*beta_trmt + rand) + theta ) )^theta),
         pmu = c(piz=0, theta=1, Intercept=0, beta_trmt=0),
         pmix=c(var=2),
         p_uppb = c(1, 9  , 10, 10, 50),
         p_lowb = c(0, 1/9,-10,-10, 0.1),

         distribution="beta-binomial-TBA",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_betabinomTBA_free_range
#
# piz     theta Intercept beta_trmt      var    value fevals gevals niter convcode kkt1 kkt2
# nlminb 0.01040509 0.4280494  2.643594 -2.377035 1.097934 575.4275     54    197    34        0   NA   NA
# xtime
# nlminb 126.956
# >

sim_binom_freerangetheta







zinbri <- glmmTMB(y~trmt +(1|feed) ,
                  ziformula=~0,
                  family=nbinom2(),
                  data=zinbri_dat
)
summary(zinbri)
# Family: nbinom2  ( log )
# Formula:          y ~ trmt + (1 | feed)
# Data: zinbri_dat
#
# AIC       BIC    logLik  deviance  df.resid
# 996911.9  996952.7 -498452.0  996903.9    199996
#
# Random effects:
#
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# feed   (Intercept) 1.046    1.023
# Number of obs: 200000, groups:  feed, 50
#
# Overdispersion parameter for nbinom2 family (): 0.432
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  2.528181   0.144717   17.47   <2e-16 ***
#   trmt        -2.301011   0.007611 -302.33   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## lock down estimates.
gnlrim(y=ybind,
       mu = ~ 1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta,
       pmu = c(Intercept=2.528181 ,beta_trmt=-2.301011,theta=0.432),
       pmix=c(var=1.046),
       p_uppb = c(2.528181,   -2.301011,0.432,  1.046),
       p_lowb = c(2.528181,   -2.301011,0.432,  1.046),
       distribution="binomial",
       nest=feed,
       random="rand",
       mixture="normal-var",
       ooo=TRUE,
       compute_hessian = FALSE,
       compute_kkt = FALSE,
       trace=1,
       method='nlminb'
)



# ## lock down estimates.
# gnlrim(y=ybind,
#        mu = ~ 1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta,
#        pmu = c(Intercept=0 ,beta_trmt=0,theta=1),
#        pmix=c(var=1.00),
#        p_uppb = c( 10,    10, 6.60,   100,  3),
#        p_lowb = c(-10,   -10,     0, -10,  0),
#        distribution="beta binomial",
#        nest=feed,
#        random="rand",
#        mixture="normal-var",
#        ooo=TRUE,
#        compute_hessian = FALSE,
#        compute_kkt = FALSE,
#        trace=1,
#        method='nlminb',
#        pshape=c(1)
# )


###
#### Try some zero-inflation?
###
##
#
## pi=0, theta=1 should give us logistic
set.seed(1112)
zinbri_dat <- createData(nc=2000, nt=2000, mu=10, beta= -2.303, sigf=1, sigc=0, pi=0.15, theta=0.43, feedNumber = 1:50)
zinbri_dat[, mean(y), by=c("feed", "trmt")]



detach(summed_binom_zinbri_dat)
summed_binom_zinbri_dat <-
  zinbri_dat[, {j=list(r=sum(y>0), n_r=sum(y==0))}, by=c("feed","trmt")]

attach(summed_binom_zinbri_dat)

ybind <- cbind(r,n_r)
x1 <- trmt

## relax theta
sim_binom_freerangetheta_zi <-
  gnlrim(y=ybind,
         mu = ~ (1-piz)*(1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta),
         pmu = c(piz=0.12, Intercept=0,beta_trmt=0,theta=1),
         pmix=c(var=20),
         p_uppb = c(1,  10,   10, 4.432 ,200.0),
         p_lowb = c(0, -10,  -10, 0.0432,   0.1),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'#,
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5


  )
# Successful convergence!
#   Save results from method  nlminb
# $par
# piz  Intercept  beta_trmt      theta        var
# 0.1639896  2.3434760 -2.2814760  0.4643283  1.2013067
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 583.2621

sim_betabinomTBA_free_range <-
  gnlrim(y=ybind,

         mu = ~ (1-piz)*(1-(theta/( exp(Intercept + x1*beta_trmt + rand) + theta ) )^theta),
         pmu = c(piz=0.12, theta=1, Intercept=0, beta_trmt=0),
         pmix=c(var=2),
         p_uppb = c(1, 9  , 10, 10, 50),
         p_lowb = c(0, 1/9,-10,-10, 0.1),

         distribution="beta-binomial-TBA",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )
# $par
# piz      theta  Intercept  beta_trmt        var
# 0.1639897  0.4643284  2.3434758 -2.2814760  1.2013063
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 583.2621


zinbri.zi <- glmmTMB(y~trmt +(1|feed) ,
                     ziformula=~1,
                     family=nbinom2(),
                     data=zinbri_dat
)
summary(zinbri.zi)
# amily: nbinom2  ( log )
# Formula:          y ~ trmt + (1 | feed)
# Zero inflation:     ~1
# Data: zinbri_dat
#
# AIC       BIC    logLik  deviance  df.resid
# 859992.4  860043.4 -429991.2  859982.4    199995
#
# Random effects:
#
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# feed   (Intercept) 1.189    1.09
# Number of obs: 200000, groups:  feed, 50
#
# Overdispersion parameter for nbinom2 family (): 0.433
#
# Conditional model:
#               Estimate Std. Error z value Pr(>|z|)
# (Intercept)  2.363158   0.154337   15.31   <2e-16 ***
#   trmt        -2.291876   0.008564 -267.63   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -1.73914    0.02839  -61.27   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > plogis(-1.73)
# [1] 0.1505876
# > exp(2.36)
# [1] 10.59095



###
#### Try some zero-inflation?
###
##
#
## pi=0, theta=1 should give us logistic
set.seed(3353)
zinbri_dat <- createData(nc=20, nt=20, mu=10, beta= -2.303, sigf=1, sigc=0, pi=0.15, theta=0.43, feedNumber = 1:3)
zinbri_dat[, mean(y), by=c("feed", "trmt")]



detach(summed_binom_zinbri_dat)
summed_binom_zinbri_dat <-
  zinbri_dat[, {j=list(r=sum(y>0), n_r=sum(y==0))}, by=c("feed","trmt")]

attach(summed_binom_zinbri_dat)

ybind <- cbind(r,n_r)
x1 <- trmt

## relax theta
sim_binom_freerangetheta_zi.small.n <-
  gnlrim(y=ybind,
         mu = ~ (1-piz)*(1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta),
         pmu = c(piz=0.12, Intercept=0,beta_trmt=0,theta=0.324),
         pmix=c(var=1.06),
         p_uppb = c(1,  10,   10,1000, 1.06),
         p_lowb = c(0, -10,  -10, 0, 1.06),
         distribution="binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'#,
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5


  )
# $par
# piz  Intercept  beta_trmt      theta        var
# 0.2633867  1.3524598 -1.2996150  4.4320000  1.4818467
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 14.77954

## relax theta
sim_binom_freerangetheta_zi.small.n.bb <-
  gnlrim(y=ybind,
         mu = ~ (1-piz)*(1-(1/( exp(Intercept + x1*beta_trmt + rand) + theta )*theta )^theta),
         pmu = c(piz=0.12, Intercept=0,beta_trmt=0,theta=0.324),
         pmix=c(var=1.06),
         p_uppb = c(1,   3,   0, 0.324,      1.06),
         p_lowb = c(0,   0,  -3, 0.324,   1.06),
         distribution="beta binomial",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         shape=~(exp(Intercept + x1*beta_trmt) + theta)^theta,
         common=TRUE
  )

zinbri.zi.small.n <- glmmTMB(y~trmt +(1|feed) ,
                     ziformula=~1,
                     family=nbinom2(),
                     data=zinbri_dat
)
summary(zinbri.zi.small.n)
# Family: nbinom2  ( log )
# Formula:          y ~ trmt + (1 | feed)
# Zero inflation:     ~1
# Data: zinbri_dat
#
# AIC      BIC   logLik deviance df.resid
# 553.6    567.6   -271.8    543.6      115
#
# Random effects:
#
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# feed   (Intercept) 1.062    1.03
# Number of obs: 120, groups:  feed, 3
#
# Overdispersion parameter for nbinom2 family (): 0.324
#
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   2.2149     0.6846   3.236  0.00121 **
#   trmt         -1.8963     0.3779  -5.018 5.23e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -3.472      7.242  -0.479    0.632


sim_betabinomTBA_free_range <-
  gnlrim(y=ybind,

         mu = ~ (1-piz)*(1-(theta/( exp(Intercept + x1*beta_trmt + rand) + theta ) )^theta),
         pmu = c(piz=0, theta=1, Intercept=0, beta_trmt=0),
         pmix=c(var=2),
         p_uppb = c(1, 900  , 10, 10, 50),
         p_lowb = c(0, 1/9,-10,-10, 0.1),

         distribution="beta-binomial-TBA",
         nest=feed,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_betabinomTBA_free_range

########################################
## 2021-12-02                         ##
## END: ZINBRI TRA/TBA   relation     ##
########################################

########################################
## 2021-11-27                         ##
## START: HGLM intercept relation     ##
########################################

## add in beta random effect
sim_mrim_data <- function(n1, n2, J, a0, a1, v, mrim="binomial-beta-HGLM", alpha=1.89, gamma=1.2, delta=1){

  if(mrim=="PPN"){
    G <- function(n,v){rnorm(n, s=sqrt(v))}
    H <- function(x) pnorm(x)
  }

  if(mrim=="LLB"){
    G <- function(n,v){bridgedist::rbridge(n, scale=1/sqrt(1+3/pi^2*v))}
    H <- function(x) plogis(x)
  }

  if(mrim=="SSS"){
    G <- function(n,v){rstable(n, alpha, 0, v, 0, 0)}
    H <- function(x) pstable(x, alpha, 0, gamma, 0, 0)
  }

  if(mrim=="CLOGLOG-CLOGLOG-LPS" & alpha!=delta){
    #library(evd)
    #source("functions.R")
    G <- function(n, v){
      ## this combo should give variance 1 mean 0
      (aa <- 1/sqrt(1+6*pi^-2*v))
      (dd <- aa*exp(-digamma(1)*(aa-1)))
      bridgecloglog_rstable(n, aa, dd)
    }
    H <- function(x){
      #1-pgumbel(x, loc=0, scale=1)
      1-exp(-exp(x))
    }
  }
  if(mrim=="CLOGLOG-CLOGLOG-LPS" & alpha==delta){
    #library(evd)
    #source("functions.R")
    G <- function(n, v){
      ## this combo should give variance 1 mean 0
      (aa <- 1/sqrt(1+6*pi^-2*v))
      (dd <- aa)
      bridgecloglog_rstable(n, aa, dd)
    }
    H <- function(x){
      #1-pgumbel(x, loc=0, scale=1)
      1-exp(-exp(x))
    }
  }
  if(mrim=="binomial-beta-HGLM"){
    G <- function(n, v){rbeta(n, v, v)}
    H <- function(x) plogis(x)
  }

  n <- n1 + n2
  u <- round(rep(G(n,v), each=J),2)

  x <- c(rep(1, n1*J), rep(0, n2*J))

  eta <- round(a0 + a1*x,2)

  eta_i <- round(eta + u,2)
  py1 <- round(H(eta_i),2)
  y <- rbinom(length(eta_i), 1, prob=py1 )

  data.frame(id=rep(1:n, each=J),
             j = rep(1:J),
             group = x,
             eta = eta,
             u_i = u,
             eta_i = eta_i,
             py1 = py1,
             y=y
  )

}
library(hglm)

detach(summed_binom_beta_dat)
set.seed(789) ## warnings for logit-bridge-phi
binom_beta_dat <-
  sim_mrim_data(400,100, J=300, a0 = -2, a1 = 1, v=1.2)
data.table::setDT(binom_beta_dat)

summed_binom_beta_dat <-
  binom_beta_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","group")]

attach(summed_binom_beta_dat)




ybind <- cbind(r,n_r)
x1 <- group

sim_binom_beta_gnlrim <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="beta-HGLM",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## but what about log(betaprime) ? it gives same intercept as hglm
log_betaprime<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + log(rand)),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="betaprime",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_binom_beta_hglm <-
  hglm(fixed = r/(r+n_r) ~ group,
       weights = (r+n_r),
       data = summed_binom_beta_dat,
       random = ~1 | id,
       family = binomial(),
       rand.family = Beta(),
       fix.disp = 1)
mean(sim_binom_beta_hglm$ranef)
(disp.hglm <- unique(sim_binom_beta_hglm$phi))
(alpha.hglm <- 0.5*(1/disp.hglm-1))
c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm)

print(
  rbind(
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## but what about Wang-Louis (2003)?
logit_bridge_var<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=0.75),
         p_uppb = c(  10,   2,   1.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

print(
  rbind(
    wang_louis_2003 = logit_bridge_var[1:4],
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)







## I want smaller alpha so I get bigger variance or random intecrpt...
## I mean, are HGLM coefs marginal?

library(hglm)

detach(summed_binom_beta_dat)
set.seed(10004)
binom_beta_dat <-
  sim_mrim_data(400,100, J=30, a0 = -2, a1 = 1, v=30, mrim = "LLB")
data.table::setDT(binom_beta_dat)

summed_binom_beta_dat <-
  binom_beta_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","group")]

attach(summed_binom_beta_dat)




ybind <- cbind(r,n_r)
x1 <- group

sim_binom_beta_gnlrim <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   1e-6),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="beta-HGLM",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## but what about log(betaprime) ? it gives same intercept as hglm
log_betaprime<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + log(rand)),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="betaprime",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_binom_beta_hglm <-
  hglm(fixed = r/(r+n_r) ~ group,
       weights = (r+n_r),
       data = summed_binom_beta_dat,
       random = ~1 | id,
       family = binomial(),
       rand.family = Beta(),
       fix.disp = 1)
mean(sim_binom_beta_hglm$ranef)
(disp.hglm <- unique(sim_binom_beta_hglm$phi))
(alpha.hglm <- 0.5*(1/disp.hglm-1))
c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm)

print(
  rbind(
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)

## but what about Wang-Louis (2003)?
logit_bridge_var<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=0.75),
         p_uppb = c(  10,   2,   1.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

print(
  rbind(
    wang_louis_2003 = logit_bridge_var[1:4],
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)





### Let's flip it.  I'm trying to get a sense
### if hglm show marginal betas.  Let's generate
### LLB sim data with a strong attenuation effect
library(hglm)

detach(summed_binom_beta_dat)
set.seed(789)
binom_beta_dat <-
  sim_mrim_data(400,100, J=300, a0 = -2, a1 = 1, v=30, mrim="LLB")## 1/sqrt(1+3/pi^2*30)= 0.314
data.table::setDT(binom_beta_dat)

summed_binom_beta_dat <-
  binom_beta_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","group")]

attach(summed_binom_beta_dat)

ybind <- cbind(r,n_r)
x1 <- group




sim_binom_beta_gnlrim <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="beta-HGLM",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## but what about log(betaprime) ? it gives same intercept as hglm
log_betaprime<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + log(rand)),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="betaprime",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_binom_beta_hglm <-
  hglm(fixed = r/(r+n_r) ~ group,
       weights = (r+n_r),
       data = summed_binom_beta_dat,
       random = ~1 | id,
       family = binomial(),
       rand.family = Beta(),
       fix.disp = 1)
mean(sim_binom_beta_hglm$ranef)
(disp.hglm <- unique(sim_binom_beta_hglm$phi))
(alpha.hglm <- 0.5*(1/disp.hglm-1))
c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm)

print(
  rbind(
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## but what about Wang-Louis (2003)?
logit_bridge_var<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=0.25),
         p_uppb = c(  10,   2,   1.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )


print(
  rbind(
    wang_louis_2003 = logit_bridge_var[1:4],
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## cool!  The marginal effect should be 0.30 * 1 = 0.30.  Look at row 2.
## and...what happened for rows 3 and 4?  They aren't similar to each other
## and they don't have a beta near 0.31.
# Intercept treatment alp1_eq_alp2 value
# wang_louis_2003              -2.01     1.070        0.298  2287
# gnlrim.rand.beta             -1.03     0.311        0.202 44485
# gnlrim.rand.logbetaprime     -2.32     0.943        0.296  2367
# hglm                         -2.07     1.215       -0.311    NA



## Again!  That was fun.  Choose a phi=0.75.
library(hglm)

detach(summed_binom_beta_dat)
set.seed(1664)
binom_beta_dat <-
  sim_mrim_data(800,200, J=300, a0 = -2, a1 = 1, v=2.5579, mrim="LLB")## 1/sqrt(1+3/pi^2*2.5579)= 0.750
data.table::setDT(binom_beta_dat)

summed_binom_beta_dat <-
  binom_beta_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","group")]

attach(summed_binom_beta_dat)

ybind <- cbind(r,n_r)
x1 <- group




sim_binom_beta_gnlrim <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="beta-HGLM",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## but what about log(betaprime) ? it gives same intercept as hglm
log_betaprime<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + log(rand)),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="betaprime",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_binom_beta_hglm <-
  hglm(fixed = r/(r+n_r) ~ group,
       weights = (r+n_r),
       data = summed_binom_beta_dat,
       random = ~1 | id,
       family = binomial(),
       rand.family = Beta(),
       fix.disp = 1)
mean(sim_binom_beta_hglm$ranef)
(disp.hglm <- unique(sim_binom_beta_hglm$phi))
(alpha.hglm <- 0.5*(1/disp.hglm-1))
c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm)

print(
  rbind(
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## but what about Wang-Louis (2003)?
logit_bridge_var<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=0,treatment=0),
         pmix=c(alp1_eq_alp2=0.25),
         p_uppb = c(  10,   2,   1.0),
         p_lowb = c( -10,  -2,   0.1),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )


print(
  rbind(
    wang_louis_2003 = logit_bridge_var[1:4],
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## 0.746 * 1.05 = 0.7833 which is close to 0.821
#                          Intercept treatment alp1_eq_alp2 value
# wang_louis_2003              -2.01     1.054        0.746  5394
# gnlrim.rand.beta             -1.92     0.821        0.269 24169
# gnlrim.rand.logbetaprime     -2.04     1.118        1.316  5417
# hglm                         -2.03     1.112        0.640    NA



## one more with strong attenuation effect, please!?
library(hglm)

detach(summed_binom_beta_dat)
set.seed(15551)
binom_beta_dat <-
  sim_mrim_data(800,200, J=300, a0 = -2, a1 = 1, v=325, mrim="LLB")## 1/sqrt(1+3/pi^2*325) = 0.100
data.table::setDT(binom_beta_dat)

summed_binom_beta_dat <-
  binom_beta_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","group")]

attach(summed_binom_beta_dat)

ybind <- cbind(r,n_r)
x1 <- group




sim_binom_beta_gnlrim <-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=-0.1,treatment=0.1),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.00),
         p_lowb = c( -10,  -2,   0.01),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="beta-HGLM",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

## but what about log(betaprime) ? it gives same intercept as hglm
log_betaprime<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + log(rand)),
         pmu = c(Intercept=-1.0,treatment=0.10),
         pmix=c(alp1_eq_alp2=2),
         p_uppb = c(  10,   2, 200.0),
         p_lowb = c( -10,  -2,   0.001),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="betaprime",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )

sim_binom_beta_hglm <-
  hglm(fixed = r/(r+n_r) ~ group,
       weights = (r+n_r),
       data = summed_binom_beta_dat,
       random = ~1 | id,
       family = binomial(),
       rand.family = Beta(),
       fix.disp = 1)
mean(sim_binom_beta_hglm$ranef)
(disp.hglm <- unique(sim_binom_beta_hglm$phi))
(alpha.hglm <- 0.5*(1/disp.hglm-1))
c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm)

print(
  rbind(
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## but what about Wang-Louis (2003)?
logit_bridge_var<-
  gnlrim(y=ybind,
         mu = ~ plogis(Intercept + x1*treatment + rand),
         pmu = c(Intercept=-1.8,treatment=1.1), ##had to smarten up
         pmix=c(alp1_eq_alp2=0.25),
         p_uppb = c(  10,   2,   1.000),
         p_lowb = c( -10,  -2,   0.001),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'
  )


print(
  rbind(
    wang_louis_2003 = logit_bridge_var[1:4],
    gnlrim.rand.beta=sim_binom_beta_gnlrim[1:4],
    gnlrim.rand.logbetaprime=log_betaprime[1:4],
    hglm=c(as.numeric(sim_binom_beta_hglm$fixef), alpha.hglm, NA)
  ),
  digits=3)
## for low=phi (strong attenuation) the hglm goes off.
# Intercept treatment alp1_eq_alp2  value
# wang_louis_2003             -1.695    1.0403       0.0889   2539
# gnlrim.rand.beta            -0.586    0.0612       0.1918 124287
# gnlrim.rand.logbetaprime    -0.730    0.1369       0.1981   3670
# hglm                        -0.914    0.6622      -0.4052     NA

########################################
## 2021-11-27                         ##
## END:   HGLM intercept relation     ##
########################################




########################################
## 2021-11-15                         ##
## START: Senn Fertility data example ##
########################################

## found the data (not where linked in paper)
## added to package.  In Makubate & Senn (DOI: 10.1002/sim.3981)
## they fit logit-normal models for 2 datasets.
## we convert them to cloglog-normal so that
## we can then make them cloglog-LogPS for
## bridging.  We also fit a logit-bridge
## model to show how the Table VII point could have
## been made.  These examples include data from the
## supplemental materials of Makubate & Senn,
## and slightly edited code.

#Analysis of Gregoriou et al data
#Gregoriou, O, Vitoratos, N, Papadias, C, Konidaris, S, Gargaropoulos, A, Rizos, D.
#Pregnancy rates in gonadotrophin stimulated cycles with timed intercourse or intrauterine insemination
#for the treatment of male subfertility, Eur J Obstet Gynecol Reprod Biol 1996; 64: 213-216.

#Load lme4 library
library(lme4)
#Remember to define filepath as required
#Read data
#NB 'period' is the period in which treatment was given
## greg <- read.table("../data_files/Greg.txt", header=T)
greg #Print data
#Proceed to fit various models
fit1 <- glmer(response~(1|patient),family=binomial,data=greg)#null model
fit2 <- glmer(response~treat+(1|patient),family=binomial,data=greg)#treatment only
fit3 <- glmer(response~period+(1|patient),family=binomial,data=greg)#period only
fit4 <- glmer(response~treat+period+(1|patient),family=binomial,data=greg)#full model
summary(fit2)#Summary of model with treatment only
summary(fit4)#Summary of model with treatment and period
#carry out analysis of deviance to check effect of adding treatment to model with period
anova(fit3,fit4)


## now change link from logit to cloglog; retain normal random intercept:
fit5 <- glmer(response~treat+(1|patient),family=binomial(link=cloglog),data=greg)#full model
summary(fit5)

## now change link from logit to cloglog; retain normal random intercept; boost nAGQ
fit5_100 <- glmer(response~treat+(1|patient),family=binomial(link=cloglog),
                  nAGQ=100,data=greg)#full model
summary(fit5_100)

## From ?glmer:
## nAGQ0initStep
## Run an initial optimization phase with nAGQ = 0.
## While the initial optimization usually provides a good starting point for
## subsequent fitting (thus increasing overall computational speed), setting
## this option to FALSE can be useful in cases where the initial phase results
## in bad fixed-effect estimates (seen most often in binomial models
## with link="cloglog" and offsets).
fit5_0 <- glmer(response~treat+(1|patient),family=binomial(link=cloglog),
                  nAGQ=100,
                  control=glmerControl(nAGQ0initStep=FALSE, tol=1e-8),
                  data=greg)#full model
summary(fit5_0)

## reproduce the cloglog-normal model below in gnlrim:
attach(greg)

y2 <- cbind(greg$response,1-greg$response)

## now do the mean-0 parameterization.
## there is no alpha or delta; just the variance.
start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
cloglog_norm <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta_treat_c*treat + rand)),
         pmu=c(beta0_c=-3, beta_treat_c=1.36),
         pmix=c(var=2.7),
         p_uppb = c(  10,   10,  110.298   ),
         p_lowb = c( -10,  -10,    0.290),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="normal-var",#"cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5
  )
cloglog_norm
summary(fit5_0)$coeff

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

## now that we have similar numbers, switch r.e. distribution to LogPS


## now do the mean-0 parameterization.
## there is no alpha or delta; just the variance.
start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_0_mean <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta_treat_c*treat + rand)),
         pmu=c(beta0_c=-3, beta_treat_c=1.36),
         pmix=c(var=2.7),
         p_uppb = c(  10,   10,  110.298   ),
         p_lowb = c( -10,  -10,    0.290),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="cloglog-bridge-0-mean",#"normal-var",#"cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
          xf.tol.nlminb=2.2e-5, ##2.2e-14,
          x.tol.nlminb=1.5e-5, ##1.5e-8,
          rel.tol.nlminb=1e-5
  )
LogPS_alpha_0_mean

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

## calculate alpha based off variance
VARIANCE <- 1.776157
(aa <- 1/sqrt(1+6*pi^-2* VARIANCE  ))
## back-calculate delta that gave 0 mean:
del <- function(bmean,aa){ aa*exp(aa*bmean-digamma(1)*(aa-1))}
dd <- del(0,aa)

## we can see what a purely marginal fit
## of these clustered data would yield;
## as a general rule of thumb
## they should be similar to a gnlrim marginal coefficients;
## but not necessarily exactly the same -- highlights the importance
## of accounting of intra-class correlation and clustering
marg_fit <-
  glm(response~#-1 + I(1-group) + group,
        1+treat,
      data=greg,
      family=binomial("cloglog"))

(bm0 <- marg_fit$coefficients[1])
(bm1 <- marg_fit$coefficients[2])


data.frame(marg_fit_intecept=bm0,
           glrim_marg_intercept=as.numeric(aa*LogPS_alpha_0_mean[1] + log(dd/aa)))
data.frame(marg_fit_slope=bm1,
           "glrim_marg_slope"=as.numeric(aa*LogPS_alpha_0_mean[2]))


## now do the alpha equal delta -- this allows to
## also do the marginal parameterization because it simplifies
## the intercept transformation
## phi here is alpha

start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_eq_delta_cond <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta_treat_c*treat + rand)),
         pmu=c(beta0_c=-3, beta_treat_c=1.36),
         pmix=c(alp_eq_del=0.5),
         p_uppb = c(  10,   10,    1-1e-5),
         p_lowb = c( -10,  -10,    0+1e-5),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="cloglog-bridge-delta-eq-alpha",#"normal-var",#"cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5
  )
finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))
LogPS_alpha_0_mean[1:6]
LogPS_alpha_eq_delta_cond[1:6]


data.frame(marg_fit_intecept=bm0,
           glrim_marg_intercept=as.numeric(aa*LogPS_alpha_eq_delta_cond[1] + log(aa/aa))) ## log(dd/aa) = 0 in this case
data.frame(marg_fit_slope=bm1,
           "glrim_marg_slope"=as.numeric(aa*LogPS_alpha_eq_delta_cond[2]))


## for marginal param:
start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_eq_delta_marg <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp( (beta0_m + beta_treat_m*treat)/alp_eq_del + rand)),
         pmu=c(beta0_m=-3, beta_treat_m=1.36,alp_eq_del=0.5),
         pmix=c(alp_eq_del=0.5),
         p_uppb = c(  10,   10,    1-1e-5),
         p_lowb = c( -10,  -10,    0+1e-5),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="cloglog-bridge-delta-eq-alpha",#"normal-var",#"cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5
  )
finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))
LogPS_alpha_0_mean[1:6]
LogPS_alpha_eq_delta_cond[1:6]
LogPS_alpha_eq_delta_marg[1:6]



data.frame(marg_fit_intecept=bm0,
           glrim_marg_intercept=as.numeric(LogPS_alpha_eq_delta_marg[1]))
data.frame(marg_fit_slope=bm1,
           "glrim_marg_slope"=as.numeric(LogPS_alpha_eq_delta_marg[2]))





## In Makubate & Senn (DOI: 10.1002/sim.3981), they carry out
## simulations to show marginalized/parallel track
## Below, we  fit a logit-bridge
## model to show how the Table VII point could have
## been made.

## Table VII could have been fit with Wang & Louis (2003) logit-bridge:
cond_parm <-
  gnlrim(y=y2,
         mu = ~plogis( (beta0_c + beta_treat_c*treat) + rand),
         pmu=c(beta0_c=-3, beta_treat_c=1.36),
         pmix=c(phi=0.5),
         p_uppb = c(  10,   10,  1-1e-5   ),
         p_lowb = c( -10,  -10,  0+1e-5),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5
  )



marg_parm <-
  gnlrim(y=y2,
         mu = ~plogis( (beta0_m + beta_treat_m*treat)/phi + rand),
         pmu=c(beta0_m=-3, beta_treat_m=1.36, phi=0.646),
         pmix=c(phi=0.646),
         p_uppb = c(  10,   10,  0.9   ),
         p_lowb = c( -10,  -10,  0.2),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="logit-bridge-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5
  )

cond_parm[1:8]
marg_parm[1:8]


detach(greg)
## above was the "greg" example from Senn Fertility
## below  is the "cohlen" example

#Analysis of Cohlen et al data
#Cohlen, BJ, Velde, ERT, Looman, CWN, Eijckemans, R, Habbema, JDF.
#Crossover or parallel design in infertility trials? The discussion continues
#Fertility and Sterility 1998; 70: 40-45.

#Load lme4 library
library(lme4)
#Remember to define filepath as required
#Read data
#NB 'term' is the period in which treatment was given
#cohlen <- read.table("../data_files/Cohlen.txt", header=T)
head(cohlen) #Print data
#Proceed to fit various models
fit1 <- glmer(response~(1|patient),family=binomial,data=cohlen)#null model
fit2 <- glmer(response~treatment+(1|patient),family=binomial,data=cohlen)#treatment only
fit3 <- glmer(response~period2 + period3 + period4 + period5 + period6 +(1|patient),family=binomial,data=cohlen)#period only, as a factor
fit4 <- glmer(response~treatment + period2 + period3 + period4 + period5 + period6+(1|patient),family=binomial,data=cohlen)#full model,period as a factor
fit5 <- glmer(response~period+(1|patient),family=binomial,data=cohlen)#period only, Period having a linear effect
fit6 <- glmer(response~treatment+period+(1|patient),family=binomial,data=cohlen)#full model,period having a linear effect
summary(fit2)#Summary of model with treatment only
summary(fit4)#Summary of model with treatment and term
#carry out analysis of deviance to check effect of adding treatment to model with term
anova(fit3,fit4)

## do cloglog

fit4_cloglog <- glmer(response~treatment + period2 + period3 + period4 + period5 + period6+(1|patient),
                      family=binomial(link=cloglog),data=cohlen)#full model
summary(fit4_cloglog)

fit4_020 <- glmer(response~treatment +
                    # period2 +
                    # period3 +
                    # period4 +
                    # period5 +
                    # period6+
                    (1|patient),
                  family=binomial(link=cloglog),data=cohlen, nAGQ=20)
summary(fit4_020)

## nAGQ0initStep
## Run an initial optimization phase with nAGQ = 0.
## While the initial optimization usually provides a good starting point for
## subsequent fitting (thus increasing overall computational speed), setting
## this option to FALSE can be useful in cases where the initial phase results
## in bad fixed-effect estimates (seen most often in binomial models
## with link="cloglog" and offsets).

# fit4_020 <- glmer(response~treatment + period2 + period3 + period4 + period5 + period6+(1|patient),
#                   family=binomial(link=cloglog),
#                   nAGQ=20,
#                   control=glmerControl(nAGQ0initStep=FALSE, tol=1e-5),
#                   data=cohlen)#full model
# summary(fit4_020)

## now do the mean-0 parameterization.
## there is no alpha or delta; just the variance.
library(gnlrim)
attach(cohlen)

y2 <- cbind(cohlen$response,1-cohlen$response)

start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
cloglog_normal_cohlen <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta_treat_c*treatment +
                            # beta_period2*period2 +
                            # beta_period3*period3 +
                            # beta_period4*period4 +
                            # beta_period5*period5 +
                            # beta_period6*period6 +
                            rand)),
         pmu=c(beta0_c=0, beta_treat_c=0#,
               # beta_period2=0,
               # beta_period3=0,
               # beta_period4=0,
               # beta_period5=0,
               # beta_period6=0
               ),
         pmix=c(scl=0.9568),
         p_uppb = c(  10,   10,  20   ),
         p_lowb = c( -10,  -10,  0+1e-5),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb'#,
         # abs.tol.nlminb=1e-4,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-4, ##2.2e-14,
         # x.tol.nlminb=1.5e-4, ##1.5e-8,
         # rel.tol.nlminb=1e-4
  )
finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

cloglog_normal_cohlen
summary(fit4_020)

start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_0_mean <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta_treat_c*treatment + rand)),
         pmu=c(beta0_c=0, beta_treat_c=0),
         pmix=c(scl=1.292),
         p_uppb = c(  10,   10,  110.298   ),
         p_lowb = c( -10,  -10,    0+1e-5),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="cloglog-bridge-0-mean",#"normal-var",#"cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-8, ##2.2e-14,
         # x.tol.nlminb=1.5e-8, ##1.5e-8,
         # rel.tol.nlminb=1e-8
  )

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

## calculate alpha based off variance
VARIANCE <- 0.07141041
(aa <- 1/sqrt(1+6*pi^-2* VARIANCE  ))
## back-calculate delta that gave 0 mean:
del <- function(bmean,aa){ aa*exp(aa*bmean-digamma(1)*(aa-1))}
dd <- del(0,aa)



marg_fit <-
  glm(response~#-1 + I(1-group) + group,
        1+treat,#+period,
      data=greg,
      family=binomial("cloglog"))

(bm0 <- marg_fit$coefficients[1])
(bm1 <- marg_fit$coefficients[2])

data.frame(marg_fit_intecept=bm0,
           glrim_marg_intercept=as.numeric(aa*LogPS_alpha_0_mean[1] + log(dd/aa)))
data.frame(marg_fit_slope=bm1,
           "glrim_marg_slope"=as.numeric(aa*LogPS_alpha_0_mean[2]))

## now do the alpha equal delta -- this allows to
## also do the marginal parameterization because it simplifies
## the intercept transformation
## phi here is alpha

start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_eq_delta_cond <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta_treat_c*treatment + rand)),
         pmu=c(beta0_c=0, beta_treat_c=0),
         pmix=c(alp_eq_del=0.9),
         p_uppb = c(  10,   10,    1-1e-5),
         p_lowb = c( -10,  -10,    0+1e-5),
         # pmu=c(beta0_c=-2.425, beta_treat_c=1.1605, beta_period_c=-0.1645),
         # pmix=c(scl=0.58),
         # p_uppb = c(-2.42, 1.161, -0.164,  0.59   ),
         # p_lowb = c(-2.43, 1.160, -0.165,  0.57   ),
         distribution="binomial",
         nest=patient,
         random="rand",
         mixture="cloglog-bridge-delta-eq-alpha",#"normal-var",#"cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method='nlminb',
         # abs.tol.nlminb=1e-8,#0,#1e-20, ## 1e-20,
         #  xf.tol.nlminb=2.2e-5, ##2.2e-14,
         #  x.tol.nlminb=1.5e-5, ##1.5e-8,
         #  rel.tol.nlminb=1e-5
  )
finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))
LogPS_alpha_0_mean[1:6]
LogPS_alpha_eq_delta_cond[1:6]

########################################
## 2021-11-15                         ##
## END: Senn Fertility data example   ##
########################################




########################################
## 2021-11-11                         ##
## START: cloglog / LogPS  free delta ##
########################################
# Derived expression for gamma
g <- function(a) {
  iu <- complex(real=0, imaginary=1)
  return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}

bridgecloglog_rstable <- function(n, alpha, delta){
  mult <- (delta/alpha)^(1/alpha)
  X <- stabledist::rstable(n , alpha, beta=1, gamma=g(alpha), delta=0, pm=1)
  Z <- log(mult * X)
  Z
}

sim_mrim_data <- function(n1, n2, J, a0, a1, v, mrim="PPN", alpha=1.89, gamma=1.2, delta=1){

  if(mrim=="PPN"){
    G <- function(n,v){rnorm(n, s=sqrt(v))}
    H <- function(x) pnorm(x)
  }

  if(mrim=="LLB"){
    G <- function(n,v){rbridge(n, scale=1/sqrt(1+3/pi^2*v))}
    H <- function(x) plogis(x)
  }

  if(mrim=="SSS"){
    G <- function(n,v){rstable(n, alpha, 0, v, 0, 0)}
    H <- function(x) pstable(x, alpha, 0, gamma, 0, 0)
  }

  if(mrim=="CLOGLOG-CLOGLOG-LPS" & alpha!=delta){
    #library(evd)
    #source("functions.R")
    G <- function(n, v){
      ## this combo should give variance 1 mean 0
      (aa <- 1/sqrt(1+6*pi^-2*v))
      (dd <- aa*exp(-digamma(1)*(aa-1)))
      bridgecloglog_rstable(n, aa, dd)
    }
    H <- function(x){
      #1-pgumbel(x, loc=0, scale=1)
      1-exp(-exp(x))
    }
  }
  if(mrim=="CLOGLOG-CLOGLOG-LPS" & alpha==delta){
    #library(evd)
    #source("functions.R")
    G <- function(n, v){
      ## this combo should give variance 1 mean 0
      (aa <- 1/sqrt(1+6*pi^-2*v))
      (dd <- aa)
      bridgecloglog_rstable(n, aa, dd)
    }
    H <- function(x){
      #1-pgumbel(x, loc=0, scale=1)
      1-exp(-exp(x))
    }
  }

  n <- n1 + n2
  u <- round(rep(G(n,v), each=J),2)

  x <- c(rep(1, n1*J), rep(0, n2*J))

  eta <- round(a0 + a1*x,2)

  eta_i <- round(eta + u,2)
  py1 <- round(H(eta_i),2)
  y <- rbinom(length(eta_i), 1, prob=py1 )

  data.frame(id=rep(1:n, each=J),
             j = rep(1:J),
             group = x,
             eta = eta,
             u_i = u,
             eta_i = eta_i,
             py1 = py1,
             y=y
  )

}

bc0 = -2.0; bc1 = 3.30; var=1.0;
## this combo should give variance 1 mean 0
(aa <- 1/sqrt(1+6*pi^-2*var))
(dd <- aa*exp(-digamma(1)*(aa-1)))
## this combo should give variance 1 mean 1
#(aa <- 1/sqrt(1+6*pi^-2*var))
#(dd <- 1.535924)

#detach(sim_data)
set.seed(102)
sim_data <- sim_mrim_data(n1=800,n2=800,J=80,a0=bc0,a1=bc1,v=var, mrim="CLOGLOG-CLOGLOG-LPS", alpha=aa, delta=dd)
summary(sim_data)
attach(sim_data)

y2 <- cbind(sim_data$y,(1-sim_data$y))

marg_fit <-
  glm(y~#-1 + I(1-group) + group,
        1+group,
      data=sim_data,
      family=binomial("cloglog"))

(bm0 <- marg_fit$coefficients[1])
(bm1 <- marg_fit$coefficients[2])

(bc0_hat <- -log(dd/aa)/aa + bm0/aa); bc0

(bc1_hat <-  bm1/aa); bc1

print(paste0("******************"))
print(paste0("******************"))
print(paste0("bc0_hat: ", round(bc0_hat,3)))
print(paste0("bc0    : ", round(bc0    ,3)))
print(paste0("******************"))
print(paste0("bm0    : ", round(bm0    ,3)))
print(paste0("******************"))
print(paste0("bc1_hat: ", round(bc1_hat,3)))
print(paste0("bc1    : ", round(bc1    ,3)))
print(paste0("******************"))
print(paste0("bm1    : ", round(bm1    ,3)))
print(paste0("******************"))
print(paste0("******************"))

alp <- function(var){1/sqrt(1+6*pi^(-2)*var)}
#bmean <- 1/aa*(log(dd/aa)+ -digamma(1)*(1-aa))
del <- function(bmean,aa){ aa*exp(aa*bmean-digamma(1)*(aa-1))}

start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_free_delta_free <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
         pmu=c(beta0_c=-2, beta1_c=3),
         pmix=c(alpha=0.80, scl=0.69),
         p_uppb = c(  10,   10,  0.95,  Inf   ),
         p_lowb = c( -10,  -10,  0.05,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="cloglog-bridge-delta-free",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb" #,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

# > start.time <- Sys.time()
# > print(paste0("Entering gnlrim at: ", start.time))
# [1] "Entering gnlrim at: 2021-11-11 14:58:02"
# > LogPS_alpha_free_delta_free <-
#   +   gnlrim(y=y2,
#              +          mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
#              +          pmu=c(beta0_c=-2, beta1_c=3),
#              +          pmix=c(alpha=0.80, scl=0.69),
#              +          p_uppb = c(  10,   10,  0.95,  Inf   ),
#              +          p_lowb = c( -10,  -10,  0.05,  0+1e-5),
#              +          distribution="binomial",
#              +          nest=id,
#              +          random="rand",
#              +          mixture="cloglog-bridge-delta-free",
#              +          ooo=TRUE,
#              +          compute_hessian = FALSE,
#              +          compute_kkt = FALSE,
#              +          trace=1,
#              +          method="nlminb" #,
#              +          # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
#                +          # xf.tol.nlminb=2.2e-7, ##2.2e-14,
#                +          # x.tol.nlminb=1.5e-7, ##1.5e-8,
#                +          # rel.tol.nlminb=1e-7
#                +   )
# fn is  fn
# Looking for method =  nlminb
# Function has  4  arguments
# par[ 1 ]:  -10   <? -2   <? 10     In Bounds
# par[ 2 ]:  -10   <? 3   <? 10     In Bounds
# par[ 3 ]:  0.05   <? 0.8   <? 0.95     In Bounds
# par[ 4 ]:  1e-05   <? 0.69   <? Inf     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6382722   log bounds ratio= 1.346787
# Method:  nlminb
# 0:     40100.997: -2.00000  3.00000 0.800000 0.690000
# 1:     40010.425: -1.95628  3.04045 0.786593 0.769199
# 2:     40005.706: -1.95624  3.04964 0.800184 0.769282
# 3:     39995.817: -1.96778  3.11038 0.811630 0.750536
# 4:     39994.488: -1.97104  3.11876 0.789264 0.745191
# 5:     39990.241: -1.97106  3.13999 0.801705 0.743053
# 6:     39987.512: -1.98223  3.18156 0.789231 0.722325
# 7:     39985.428: -1.98196  3.22603 0.803275 0.722835
# 8:     39983.219: -2.00469  3.26200 0.792262 0.707217
# 9:     39983.005: -2.00896  3.26518 0.793542 0.714551
# 10:     39982.986: -2.01089  3.26379 0.794580 0.713355
# 11:     39982.968: -2.01109  3.26180 0.794439 0.715388
# 12:     39982.968: -2.01112  3.26161 0.794238 0.715341
# 13:     39982.967: -2.01109  3.26144 0.794399 0.715498
# 14:     39982.966: -2.01121  3.26091 0.794300 0.715631
# 15:     39982.965: -2.01126  3.25990 0.794572 0.716106
# 16:     39982.964: -2.01022  3.25964 0.794367 0.715793
# 17:     39982.964: -2.00935  3.25902 0.794414 0.715374
# 18:     39982.964: -2.00831  3.25887 0.794527 0.714947
# 19:     39982.964: -2.00733  3.25884 0.794484 0.714365
# 20:     39982.964: -2.00477  3.25883 0.794483 0.712908
# 21:     39982.964: -1.99451  3.25882 0.794481 0.707088
# 22:     39982.964: -1.99431  3.25882 0.794481 0.706997
# 23:     39982.964: -1.98438  3.25882 0.794481 0.701443
# 24:     39982.964: -1.97401  3.25885 0.794484 0.695661
# 25:     39982.964: -1.96794  3.25889 0.794490 0.692337
# 26:     39982.964: -1.96360  3.25888 0.794490 0.689946
# 27:     39982.964: -1.95926  3.25886 0.794488 0.687564
# 28:     39982.964: -1.95353  3.25880 0.794482 0.684459
# 29:     39982.964: -1.95325  3.25879 0.794478 0.684329
# 30:     39982.964: -1.94734  3.25883 0.794478 0.681115
# 31:     39982.964: -1.93615  3.25884 0.794485 0.675046
# 32:     39982.964: -1.93987  3.25886 0.794489 0.677084
# 33:     39982.964: -1.93688  3.25885 0.794486 0.675475
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# beta0_c    beta1_c      alpha        scl
# -1.9368846  3.2588463  0.7944856  0.6754750
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 39982.96
#
# $fevals
# function
# 46
#
# $gevals
# gradient
# 225
#
# $nitns
# [1] 33
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 8692.615
#
# Assemble the answers
# Warning message:
#   In nlminb(start = par, objective = ufn, gradient = ugr, lower = lower,  :
#               NA/NaN function evaluation
#             >
#               > finish.time <- Sys.time()
#             > total.time <- finish.time - start.time
#             > total.time
#             Time difference of 55.55566 mins
#             > print(paste0("gnlrim took: ", total.time))
#             [1] "gnlrim took: 55.5556649168332"


# ## -aa*bc0_hat - log(dd/aa) ==  bm0
# > -0.7944856 * -1.9368846 - log(0.6754750/0.7944856)
# [1] 1.701106
#
# > ## beta1_c times alpha_hat
#   > 3.2588463  * 0.7944856
# [1] 2.589106



## now do the mean-0 parameterization.
## there is no alpha or delta; just the variance.
start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_0_mean <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
         pmu=c(beta0_c=-2, beta1_c=3),
         pmix=c(scl=pi^2*(0.80^-2-1)/6), ## alpha=0.80
         p_uppb = c(  10,   10,  Inf   ),
         p_lowb = c( -10,  -10,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="cloglog-bridge-0-mean",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb" #,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

# [1] "Entering gnlrim at: 2021-11-11 18:07:06"
# > LogPS_alpha_0_mean <-
#   +   gnlrim(y=y2,
#              +          mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
#              +          pmu=c(beta0_c=-2, beta1_c=3),
#              +          pmix=c(scl=pi^2*(0.80^-2-1)/6), ## alpha=0.80
#              +          p_uppb = c(  10,   10,  Inf   ),
#              +          p_lowb = c( -10,  -10,  0+1e-5),
#              +          distribution="binomial",
#              +          nest=id,
#              +          random="rand",
#              +          mixture="cloglog-bridge-0-mean",
#              +          ooo=TRUE,
#              +          compute_hessian = FALSE,
#              +          compute_kkt = FALSE,
#              +          trace=1,
#              +          method="nlminb" #,
#              +          # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
#                +          # xf.tol.nlminb=2.2e-7, ##2.2e-14,
#                +          # x.tol.nlminb=1.5e-7, ##1.5e-8,
#                +          # rel.tol.nlminb=1e-7
#                +   )
# fn is  fn
# Looking for method =  nlminb
# Function has  3  arguments
# par[ 1 ]:  -10   <? -2   <? 10     In Bounds
# par[ 2 ]:  -10   <? 3   <? 10     In Bounds
# par[ 3 ]:  1e-05   <? 0.9252754   <? Inf     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.5108502   log bounds ratio= 0
# Method:  nlminb
# 0:     40072.639: -2.00000  3.00000 0.925275
# 1:     40005.516: -1.93126  3.06911 0.902950
# 2:     39993.509: -1.91198  3.16665 0.892288
# 3:     39985.026: -1.98031  3.20920 0.951620
# 4:     39983.526: -1.97748  3.25743 0.957456
# 5:     39983.516: -2.02014  3.28073 0.955241
# 6:     39983.189: -2.00257  3.28050 0.972062
# 7:     39983.126: -2.00163  3.26995 0.970668
# 8:     39983.033: -1.99719  3.26940 0.960971
# 9:     39983.009: -2.00007  3.26222 0.953609
# 10:     39982.968: -1.99360  3.26158 0.962087
# 11:     39982.965: -1.99211  3.26007 0.961439
# 12:     39982.965: -1.99251  3.25956 0.961618
# 13:     39982.964: -1.99189  3.25937 0.961810
# 14:     39982.964: -1.99173  3.25895 0.961310
# 15:     39982.964: -1.99190  3.25888 0.960995
# 16:     39982.964: -1.99183  3.25881 0.961049
# 17:     39982.964: -1.99183  3.25881 0.961049
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# beta0_c    beta1_c        scl
# -1.9918274  3.2588144  0.9610486
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 39982.96
#
# $fevals
# function
# 24
#
# $gevals
# gradient
# 70
#
# $nitns
# [1] 17
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 3039.317
#
# Assemble the answers
# Warning message:
#   In nlminb(start = par, objective = ufn, gradient = ugr, lower = lower,  :
#               NA/NaN function evaluation
#             >
#               > finish.time <- Sys.time()
#             > total.time <- finish.time - start.time
#             > total.time
#             Time difference of 19.53819 mins
#             > print(paste0("gnlrim took: ", total.time))
#             [1] "gnlrim took: 19.5381904681524"

# > (aa <- 1/sqrt(1+6*pi^-2*0.9611163 ))
# [1] 0.7944798

# which means delta is del <- function(bmean,aa){ aa*exp(aa*bmean-digamma(1)*(aa-1))}
# R> 0.7944798*exp(0.7944798*0-digamma(1)*(0.7944798-1))
# [1] 0.7056068

# ## -aa*bc0_hat - log(dd/aa) ==  bm0
# > -0.7944798 * -1.9917919 - log(0.7056068/0.7944798)
# > -0.7944798 * -1.9917919 - log(0.7056068/0.7944798)
# [1] 1.701068
#
# > ## beta1_c times alpha_hat
#   > 3.2588396   * 0.7944798
# [1] 2.589082



## now do delta == alpha param
start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_eq_delta <-
  gnlrim(y=y2,
         mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
         pmu=c(beta0_c=-2, beta1_c=3),
         pmix=c(scl=0.80),
         p_uppb = c(  10,   10,  1-1e-5),
         p_lowb = c( -10,  -10,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="cloglog-bridge-delta-eq-alpha",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb" #,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))

# fn is  fn
# Looking for method =  nlminb
# Function has  3  arguments
# par[ 1 ]:  -10   <? -2   <? 10     In Bounds
# par[ 2 ]:  -10   <? 3   <? 10     In Bounds
# par[ 3 ]:  1e-05   <? 0.8   <? 0.99999     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.5740313   log bounds ratio= 1.301039
# Method:  nlminb
# 0:     40018.722: -2.00000  3.00000 0.800000
# 1:     40014.883: -1.99692  3.01256 0.807203
# 2:     40009.836: -1.99999  3.04061 0.789127
# 3:     40001.696: -2.01098  3.10437 0.806594
# 4:     39990.817: -2.06636  3.13962 0.793089
# 5:     39987.376: -2.11448  3.18342 0.809120
# 6:     39984.835: -2.09857  3.24394 0.785132
# 7:     39984.496: -2.16266  3.26239 0.791643
# 8:     39983.368: -2.16044  3.26394 0.797087
# 9:     39983.208: -2.15461  3.26263 0.798228
# 10:     39983.096: -2.14598  3.25478 0.794790
# 11:     39983.013: -2.13396  3.25288 0.794762
# 12:     39983.006: -2.13947  3.26357 0.792925
# 13:     39982.993: -2.13980  3.25751 0.793340
# 14:     39982.971: -2.14175  3.25843 0.795188
# 15:     39982.965: -2.14099  3.25884 0.794457
# 16:     39982.965: -2.14117  3.25885 0.794485
# 17:     39982.965: -2.14113  3.25884 0.794483
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# beta0_c   beta1_c       scl
# -2.141126  3.258838  0.794483
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 39982.96
#
# $fevals
# function
# 23
#
# $gevals
# gradient
# 63
#
# $nitns
# [1] 17
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 2877.638
#
# Assemble the answers
# Warning message:
#   In nlminb(start = par, objective = ufn, gradient = ugr, lower = lower,  :
#               NA/NaN function evaluation
#             >
#               > finish.time <- Sys.time()
#             > total.time <- finish.time - start.time
#             > total.time
#             Time difference of 19.82804 mins
#             > print(paste0("gnlrim took: ", total.time))
#             [1] "gnlrim took: 19.8280424316724"
#             >

library(data.table)
setDT(sim_data)
sim_data2 <-
      sim_data[,{j=list(y=sum(y),n_y=sum(y==0))}, by=c("id","group")]
detach(sim_data)
attach(sim_data2)

y3 <- cbind(sim_data2$y, sim_data2$n_y)

## now do delta == alpha param
start.time <- Sys.time()
print(paste0("Entering gnlrim at: ", start.time))
LogPS_alpha_eq_delta <-
  gnlrim(y=y3,
         mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
         pmu=c(beta0_c=-2, beta1_c=3),
         pmix=c(scl=0.80),
         p_uppb = c(  10,   10,  1-1e-5),
         p_lowb = c( -10,  -10,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="cloglog-bridge-delta-eq-alpha",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb" #,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

finish.time <- Sys.time()
total.time <- finish.time - start.time
total.time
print(paste0("gnlrim took: ", total.time))
#
# [1] "Entering gnlrim at: 2021-11-17 20:08:34"
# > LogPS_alpha_eq_delta <-
#   +   gnlrim(y=y3,
#              +          mu = ~1-exp(-exp(beta0_c + beta1_c*group + rand)),
#              +          pmu=c(beta0_c=-2, beta1_c=3),
#              +          pmix=c(scl=0.80),
#              +          p_uppb = c(  10,   10,  1-1e-5),
#              +          p_lowb = c( -10,  -10,  0+1e-5),
#              +          distribution="binomial",
#              +          nest=id,
#              +          random="rand",
#              +          mixture="cloglog-bridge-delta-eq-alpha",
#              +          ooo=TRUE,
#              +          compute_hessian = FALSE,
#              +          compute_kkt = FALSE,
#              +          trace=1,
#              +          method="nlminb" #,
#              +          # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
#                +          # xf.tol.nlminb=2.2e-7, ##2.2e-14,
#                +          # x.tol.nlminb=1.5e-7, ##1.5e-8,
#                +          # rel.tol.nlminb=1e-7
#                +   )
# fn is  fn
# Looking for method =  nlminb
# Function has  3  arguments
# par[ 1 ]:  -10   <? -2   <? 10     In Bounds
# par[ 2 ]:  -10   <? 3   <? 10     In Bounds
# par[ 3 ]:  1e-05   <? 0.8   <? 0.99999     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.5740313   log bounds ratio= 1.301039
# Method:  nlminb
# 0:     4981.8824: -2.00000  3.00000 0.800000
# 1:     4978.0429: -1.99692  3.01256 0.807203
# 2:     4973.1256: -1.99989  3.03943 0.789607
# 3:     4965.0543: -2.00959  3.10095 0.806407
# 4:     4954.4833: -2.06224  3.13568 0.792867
# 5:     4950.8574: -2.10893  3.17724 0.808780
# 6:     4948.0474: -2.09477  3.23584 0.785844
# 7:     4947.0386: -2.15558  3.25611 0.793065
# 8:     4946.4060: -2.15313  3.25794 0.798747
# 9:     4946.2627: -2.14807  3.25489 0.796160
# 10:     4946.2140: -2.13519  3.25530 0.795507
# 11:     4946.1678: -2.14202  3.26595 0.793002
# 12:     4946.1650: -2.13682  3.25415 0.793324
# 13:     4946.1498: -2.13840  3.25496 0.795482
# 14:     4946.1268: -2.13995  3.25715 0.794683
# 15:     4946.1249: -2.14116  3.25900 0.794472
# 16:     4946.1249: -2.14114  3.25882 0.794486
# 17:     4946.1249: -2.14112  3.25884 0.794483
# 18:     4946.1249: -2.14113  3.25884 0.794483
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# beta0_c    beta1_c        scl
# -2.1411262  3.2588387  0.7944831
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 4946.125
#
# $fevals
# function
# 24
#
# $gevals
# gradient
# 69
#
# $nitns
# [1] 18
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 2336.767
#
# Assemble the answers
# Warning message:
#   In nlminb(start = par, objective = ufn, gradient = ugr, lower = lower,  :
#               NA/NaN function evaluation
#             > finish.time <- Sys.time()
#             > total.time <- finish.time - start.time
#             > total.time
#             Time difference of 8.238013 mins
#             > print(paste0("gnlrim took: ", total.time))
#             [1] "gnlrim took: 8.23801318407059"



########################################
## 2021-11-11                         ##
## ENDxx: cloglog / LogPS  free delta ##
########################################

####################
## Pre 2021-11-11 ##
####################

sim_mrim_data2 <- function(n1, n2, J, a0, a1, v= 6/sqrt(1.69), mrim="SSS", alpha=1.69, gamma= 1/sqrt(1.69)){

  if(mrim=="PPN"){
    G <- function(n,v){rnorm(n, s=sqrt(v))}
    H <- function(x) pnorm(x)
  }

  if(mrim=="LLB"){
    G <- function(n,v){rbridge(n, scale=1/sqrt(1+3/pi^2*v))}
    H <- function(x) plogis(x)
  }

  if(mrim=="SSS"){
    G <- function(n,v){stabledist::rstable(n, alpha, 0, v, 0, 0)}
    H <- function(x) stabledist::pstable(x, alpha, 0, gamma, 0, 0)
  }


  n <- n1 + n2
  u <- round(rep(G(n,v), each=J),2)

  x <- c(rep(1, n1*J), rep(0, n2*J))

  eta <- round(a0 + a1*x,2)

  eta_i <- round(eta + u,2)
  py1 <- round(H(eta_i),2)
  y <- rbinom(length(eta_i), 1, prob=py1 )

  data.frame(id=rep(1:n, each=J),
             j = rep(1:J),
             group = x,
             eta = eta,
             u_i = u,
             eta_i = eta_i,
             py1 = py1,
             y=y
  )

}


#############3
JJ <- 60
a1.true <-  0.8
a0.true <- -1
alp.true= 1.69
gam.g.true <- 4

phi.true <- (1/sqrt(alp.true)) /( (1/sqrt(alp.true))^alp.true + (gam.g.true/sqrt(alp.true))^alp.true)^(1/alp.true)
bm0.true <- a0.true * phi.true
bm1.true <- a1.true * phi.true
(p0.true <- stable_cdf2( bm0.true           , c(alp.true, 0, 1/sqrt(alp.true), 0)))
(p1.true <- stable_cdf2( bm0.true + bm1.true, c(alp.true, 0, 1/sqrt(alp.true), 0)))

(tba.true <- 100*(1-p0.true/p1.true))

dputted_dat<-
  dput(data.table::data.table(sim_mrim_data2(n1=25,n2=100,J=JJ,a0=a0.true,a1=a1.true, v = gam.g.true / sqrt(alp.true), mrim="SSS", alpha=1.69, gam=1/sqrt(1.69)))[,sum(y), by=c("id","group")])

robby <- data.table::data.table(id=dputted_dat$id, group=dputted_dat$group, y1=dputted_dat$V1, y0=JJ-dputted_dat$V1); data.table::setorder(robby, group, y1, y0); robby

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(ggplot2::aes(color=factor(group)))

dose <- robby$group

y_01 <- robby$y1

id <- robby$id

y_cbind = cbind(y_01, JJ-y_01)

cbind(id, y_cbind, dose)

gapr_marg_alpha_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=0.50),
         p_uppb = c(  Inf,  Inf,  1.90,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1.45,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb" #,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

gapr_marg_alpha_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb" #,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )


## refit with confidence intervals:
fu_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_2.00_ci <-
  fu_2.00$coefficients[2] + 1.96*c(-1,1)*fu_2.00$se[2]
fu_2.00_ci

if(fu_2.00_ci[1] <0){print("Potential WINNER:  2.00 has confidence interval that includes 0:"); fu_2.00_ci
}


fu_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=5.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  0.1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_free_ci <-
  fu_free$coefficients[2] + 1.96*c(-1,1)*fu_free$se[2]
fu_free_ci

if(fu_free_ci[1] >0){print("Potential WINNER:  free has confidence interval that EXcludes 0:"); fu_free_ci
}

###########################################################################
## Trying to work in the gamma_G/sqrt(A) settings for mixture. See      ###
## if this section provides same estimates as the one immediately below.###
##  Trying to get 11 in each quartile.
##               ###
##                ###
###########################################################################

library(data.table)
JJ <- 100
dputted_dat <-
  structure(list(id = 1:88,
                 group=c(c(    0, #max, group 0

                               0,
                             0,0,0,
                           0,0,0,0,0,

                               0, #75+, group 0
                               0, #75-, group 0

                               0,
                             0,0,0,
                           0,0,0,0,0,

                               0, #50+, group 0
                               0, #50-, group 0

                               0,
                             0,0,0,
                           0,0,0,0,0,

                               0, #25+, group 0
                               0, #25-, group 0

                               0,
                             0,0,0,
                           0,0,0,0,0,

                               0 #min, group 0
                 ),
                         c(    1, #max, group 1

                               1,
                             1,1,1,
                           1,1,1,1,1,

                               1, #75+, group 1
                               1, #75-, group 1

                               1,
                             1,1,1,
                           1,1,1,1,1,

                               1, #51+, group 1
                               1, #51-, group 1

                               1,
                             1,1,1,
                           1,1,1,1,1,

                               1, #25+, group 1
                               1, #25-, group 1

                               1,
                             1,1,1,
                           1,1,1,1,1,

                               1 #min, group 1
                 )),
              quartile=c(c(    4, #max, group 0

                               4,
                             4,4,4,
                           4,4,4,4,4,

                               4, #75+, group 0
                               3, #75-, group 0

                               3,
                             3,3,3,
                           3,3,3,3,3,

                               3, #50+, group 0
                               2, #50-, group 0

                               2,
                             2,2,2,
                           2,2,2,2,2,

                               2, #25+, group 0
                               1, #25-, group 0

                               1,
                             1,1,1,
                           1,1,1,1,1,

                               1 #min, group 0
                 ),
                   c(    4, #max, group 1

                               4,
                             4,4,4,
                           4,4,4,4,4,

                               4, #75+, group 1
                               3, #75-, group 1

                               3,
                             3,3,3,
                           3,3,3,3,3,

                               3, #50+, group 1
                               2, #50-, group 1

                               2,
                             2,2,2,
                           2,2,2,2,2,

                               2, #25+, group 1
                               1, #25-, group 1

                               1,
                             1,1,1,
                           1,1,1,1,1,

                               1 #min, group 1
                 )),
                 V1   =c(c(      99, #max, group 0

                                 87,
                              85,85,85,
                           83,83,83,83,83,

                                 76, #75+, group 0
                                 74, #75-, group 0

                                 67,
                              65,65,65,
                           63,63,63,63,63,

                                 21, #50+, group 0
                                 19, #50-, group 0

                                 17,
                              15,15,15,
                           13,13,13,13,13,

                                  11, #25+, group 0
                                  9, #25-, group 0

                                1,
                              1,1,1,
                           1,1,1,1,1,

                                  0 #min, group 0
                 ),
                         c(      99, #max, group 1

                           97,97,97,97,97,
                              95,95,95,
                                 93,

                                 76, #75+, group 1
                                 74, #75-, group 1

                           c(c(67,67,67,67,67,
                                65,65,65,
                                   63,

                               61, #50+, group 1
                               59, #50-, group 1

                           57,57,57,57,57,
                             55,55,55,
                               53) -7),

                               11, #25+, group 1
                                9, #25-, group 1

                           7,7,7,7,7,
                             5,5,5,
                               3,

                               0 #min, group 1
                 ))/100))

robby <- data.table::data.table(id=dputted_dat$id,
                                group=dputted_dat$group,
                                y1=JJ*dputted_dat$V1,
                                y0=JJ-JJ*dputted_dat$V1,
                                quartile=dputted_dat$quartile);
data.table::setorder(robby, group, y1, y0); robby

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) +
  ggplot2::geom_boxplot() +
  ggbeeswarm::geom_beeswarm(ggplot2::aes(color=factor(quartile)))

dose <- dputted_dat$group

y_01 <- JJ*dputted_dat$V1

id <- dputted_dat$id

y_cbind = cbind(y_01, JJ-y_01)

cbind(id, y_cbind, dose)

gapr_marg_alpha_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.50),
         pmix=c(alpha=1.5, scl=5.50),
         p_uppb = c(  Inf,  Inf,  2.0,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1.0,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"#,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

gapr_marg_alpha_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_1.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1),
         pmix=c(alpha=1, scl=0.50),
         p_uppb = c(  Inf,  Inf,  1,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

roger <-
  rbind(gapr_marg_alpha_2.00,gapr_marg_alpha_free,gapr_marg_alpha_1.00)
rownames(roger) <- c("alpha set at 2.00", "alpha free range", "alpha set at 1.00")
roger

a_free <- gapr_marg_alpha_free$alpha
phi_free <- (1/sqrt(a_free)) /( (1/sqrt(a_free))^a_free + (gapr_marg_alpha_free$scl/sqrt(a_free))^a_free)^(1/a_free)
a_2.00 <- gapr_marg_alpha_2.00$alpha
phi_2.00 <- (1/sqrt(a_2.00)) /( (1/sqrt(a_2.00))^a_2.00 + (gapr_marg_alpha_2.00$scl/sqrt(a_2.00))^a_2.00)^(1/a_2.00)


bm0_free <- gapr_marg_alpha_free$a0 * phi_free
bm1_free <- gapr_marg_alpha_free$a1 * phi_free
p0_free <- stable_cdf2( bm0_free           , c(a_free, 0, 1/sqrt(a_free), 0))
p1_free <- stable_cdf2( bm0_free + bm1_free, c(a_free, 0, 1/sqrt(a_free), 0))

bm0_2.00 <- gapr_marg_alpha_2.00$a0 * phi_2.00
bm1_2.00 <- gapr_marg_alpha_2.00$a1 * phi_2.00
p0_2.00 <- stable_cdf2( bm0_2.00           , c(a_2.00, 0, 1/sqrt(a_2.00), 0))
p1_2.00 <- stable_cdf2( bm0_2.00 + bm1_2.00, c(a_2.00, 0, 1/sqrt(a_2.00), 0))

## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)

### mse: a_2.00 performs better
mean(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
mean(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
## median square error: a_free performs better
median(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
median(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
### mean absolute error: a_free performs better
mean(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
mean(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))
## median absolute error: : a_free performs better
median(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
median(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))

## refit with confidence intervals:
fu_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_2.00_ci <-
  fu_2.00$coefficients[2] + 1.96*c(-1,1)*fu_2.00$se[2]
fu_2.00_ci

## these are conditional probabilities!
stable_cdf2(fu_2.00$coefficients[1], c(2, 0, 1/sqrt(2), 0))
stable_cdf2(fu_2.00$coefficients[1] +fu_2.00$coefficients[2], c(2, 0, 1/sqrt(2), 0))
## fitted values
sort(unique(fu_2.00$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
## mean(c((robby[group==1]$y1 - 38.41672)^2, (robby[group==0]$y1 - 26.96960)^2 ))

if(fu_2.00_ci[1] <0){print("Potential WINNER:  2.00 has confidence interval that includes 0:"); fu_2.00_ci
}


fu_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=5.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  0.1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
         )

fu_free_ci <-
  fu_free$coefficients[2] + 1.96*c(-1,1)*fu_free$se[2]
fu_free_ci

## these are conditional probabilities!
stable_cdf2(fu_free$coefficients[1], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
stable_cdf2(fu_free$coefficients[1] +fu_free$coefficients[2], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
## fitted values
sort(unique(fu_free$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
#mean(c((robby[group==1]$y1 - 34.653820)^2, (robby[group==0]$y1 - 4.591598)^2 ))


if(fu_free_ci[1] >0){print("Potential WINNER:  free has confidence interval that EXcludes 0:"); fu_free_ci
}



## next two lines same?: no.  these are the conditional probabilities * 100
sort(unique(fu_2.00$fitted.values))
sort(unique(fu_free$fitted.values))
## next two lines same?: no.  these are the marginal probabilities * 100
JJ*c(p0_2.00,p1_2.00)
JJ*c(p0_free,p1_free)
## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)


#
#
#   > sort(unique(fu_2.00$fitted.values))
# [1] 26.96960 38.41672
#   > JJ*c(p0_2.00,p1_2.00)
# [1] 37.52979 43.93764
# >
#   >
#   > sort(unique(fu_free$fitted.values))
# [1]  4.591598 34.653820
# > JJ*c(p0_free,p1_free)
# [1] 20.67212 45.95814

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(color="blue")


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
###########################################################################
## Trying to work in the gamma_G/sqrt(A) settings for mixture. See      ###
## if this section provides same estimates as the one immediately below.###
## A:  No, it doesn't for alpha_free                                    ###
## A': Oh, wait, yes it does if you change starting value for scl!      ###
###########################################################################

library(data.table)
JJ <- 100
dputted_dat <-
  structure(list(id = 1:48,
                 group = c(1, 1, 1, 1, 1, 1,
                           1,
                           1, 1, 1, 1,
                           1,
                           1,
                           1, 1, 1, 1,
                           1,
                           1, 1, 1, 1, 1, 1,
                           0, 0, 0, 0, 0, 0,
                           0,
                           0, 0, 0, 0,
                           0,
                           0,
                           0, 0, 0, 0,
                           0,
                           0, 0, 0, 0, 0, 0),
                 V1 = c(1, 1, 1, 3, 3, 3,
                        5,
                        c(c(7, 7, 7, 9) + 30),
                        55,
                        59,
                        65, 68, 78, 88,
                        90,
                        97, 97, 97, 99, 99, 99,
                        1, 1, 1, 3, 3, 3,
                        5,
                        7, 9, 9, 9,
                        11,
                        13,
                        c(c(65, 68, 78, 88) - 25),
                        90,
                        97, 97, 97, 99, 99, 99)/100))

robby <- data.table::data.table(id=dputted_dat$id, group=dputted_dat$group, y1=JJ*dputted_dat$V1, y0=JJ-JJ*dputted_dat$V1); data.table::setorder(robby, group, y1, y0); robby

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(color="blue")

dose <- dputted_dat$group

y_01 <- JJ*dputted_dat$V1

id <- dputted_dat$id

y_cbind = cbind(y_01, JJ-y_01)

cbind(id, y_cbind, dose)

gapr_marg_alpha_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.50),
         pmix=c(alpha=1.5, scl=5.50),
         p_uppb = c(  Inf,  Inf,  2.0,  Inf   ),
         p_lowb = c( -Inf, -Inf,  0.1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"#,
         # abs.tol.nlminb=1e-7,#0,#1e-20, ## 1e-20,
         # xf.tol.nlminb=2.2e-7, ##2.2e-14,
         # x.tol.nlminb=1.5e-7, ##1.5e-8,
         # rel.tol.nlminb=1e-7
  )

gapr_marg_alpha_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )
## should match _2.00 -- it does b/c _2.00's sqrt(scl) == norm_norm_norm's scl
norm_norm_norm <-
  gnlrim(y=y_cbind,
         mu=~ pnorm(a0 + a1*dose + rand),

         pmu=c( a0=0, a1=0),
         pmix=c(scl=0.50),
         p_uppb = c(  Inf,  Inf,  Inf   ),
         p_lowb = c( -Inf, -Inf,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="normal-var",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_1.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1),
         pmix=c(alpha=1, scl=0.50),
         p_uppb = c(  Inf,  Inf,  1,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )
## should match _1.00
cauchy_cauchy_cauchy <-
  gnlrim(y=y_cbind,
         mu=~ pcauchy(a0 + a1*dose + rand),

         pmu=c( a0=0, a1=0),
         pmix=c(scl=0.50),
         p_uppb = c(  Inf,  Inf,  Inf   ),
         p_lowb = c( -Inf, -Inf,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="Cauchy-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

roger <-
  rbind(gapr_marg_alpha_2.00,gapr_marg_alpha_free,gapr_marg_alpha_1.00)
rownames(roger) <- c("alpha set at 2.00", "alpha free range", "alpha set at 1.00")
roger

a_free <- gapr_marg_alpha_free$alpha
phi_free <- (1/sqrt(a_free)) /( (1/sqrt(a_free))^a_free + (gapr_marg_alpha_free$scl/sqrt(a_free))^a_free)^(1/a_free)
a_2.00 <- gapr_marg_alpha_2.00$alpha
phi_2.00 <- (1/sqrt(a_2.00)) /( (1/sqrt(a_2.00))^a_2.00 + (gapr_marg_alpha_2.00$scl/sqrt(a_2.00))^a_2.00)^(1/a_2.00)


bm0_free <- gapr_marg_alpha_free$a0 * phi_free
bm1_free <- gapr_marg_alpha_free$a1 * phi_free
p0_free <- stable_cdf2( bm0_free           , c(a_free, 0, 1/sqrt(a_free), 0))
p1_free <- stable_cdf2( bm0_free + bm1_free, c(a_free, 0, 1/sqrt(a_free), 0))

bm0_2.00 <- gapr_marg_alpha_2.00$a0 * phi_2.00
bm1_2.00 <- gapr_marg_alpha_2.00$a1 * phi_2.00
p0_2.00 <- stable_cdf2( bm0_2.00           , c(a_2.00, 0, 1/sqrt(a_2.00), 0))
p1_2.00 <- stable_cdf2( bm0_2.00 + bm1_2.00, c(a_2.00, 0, 1/sqrt(a_2.00), 0))

## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)

### mse: a_2.00 performs better
mean(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
mean(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
## median square error: a_free performs better
median(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
median(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
### mean absolute error: a_free performs better
mean(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
mean(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))
## median absolute error: : a_free performs better
median(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
median(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))

## refit with confidence intervals:
fu_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_2.00_ci <-
  fu_2.00$coefficients[2] + 1.96*c(-1,1)*fu_2.00$se[2]
fu_2.00_ci

## these are conditional probabilities!
stable_cdf2(fu_2.00$coefficients[1], c(2, 0, 1/sqrt(2), 0))
stable_cdf2(fu_2.00$coefficients[1] +fu_2.00$coefficients[2], c(2, 0, 1/sqrt(2), 0))
## fitted values
sort(unique(fu_2.00$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
## mean(c((robby[group==1]$y1 - 38.41672)^2, (robby[group==0]$y1 - 26.96960)^2 ))

if(fu_2.00_ci[1] <0){print("Potential WINNER:  2.00 has confidence interval that includes 0:"); fu_2.00_ci
}


fu_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=5.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  0.1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl-over-sqrt",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
         )

fu_free_ci <-
  fu_free$coefficients[2] + 1.96*c(-1,1)*fu_free$se[2]
fu_free_ci

## these are conditional probabilities!
stable_cdf2(fu_free$coefficients[1], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
stable_cdf2(fu_free$coefficients[1] +fu_free$coefficients[2], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
## fitted values
sort(unique(fu_free$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
#mean(c((robby[group==1]$y1 - 34.653820)^2, (robby[group==0]$y1 - 4.591598)^2 ))


if(fu_free_ci[1] >0){print("Potential WINNER:  free has confidence interval that EXcludes 0:"); fu_free_ci
}



## next two lines same?: no.  these are the conditional probabilities * 100
sort(unique(fu_2.00$fitted.values))
sort(unique(fu_free$fitted.values))
## next two lines same?: no.  these are the marginal probabilities * 100
JJ*c(p0_2.00,p1_2.00)
JJ*c(p0_free,p1_free)
## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)


#
#
#   > sort(unique(fu_2.00$fitted.values))
# [1] 26.96960 38.41672
#   > JJ*c(p0_2.00,p1_2.00)
# [1] 37.52979 43.93764
# >
#   >
#   > sort(unique(fu_free$fitted.values))
# [1]  4.591598 34.653820
# > JJ*c(p0_free,p1_free)
# [1] 20.67212 45.95814

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(color="blue")


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
library(data.table)
JJ <- 100
dputted_dat <-
  structure(list(id = 1:48,
                 group = c(1, 1, 1, 1, 1, 1,
                           1,
                           1, 1, 1, 1,
                           1,
                           1,
                           1, 1, 1, 1,
                           1,
                           1, 1, 1, 1, 1, 1,
                           0, 0, 0, 0, 0, 0,
                           0,
                           0, 0, 0, 0,
                           0,
                           0,
                           0, 0, 0, 0,
                           0,
                           0, 0, 0, 0, 0, 0),
                 V1 = c(1, 1, 1, 3, 3, 3,
                        5,
                        c(c(7, 7, 7, 9) + 30),
                        55,
                        59,
                        65, 68, 78, 88,
                        90,
                        97, 97, 97, 99, 99, 99,
                        1, 1, 1, 3, 3, 3,
                        5,
                        5,5,5,5,#7, 9, 9, 9,
                        5,#11,
                        6,#13,
                        8,13,15,66,##c(c(65, 68, 78, 88) - 25),
                        90,
                        97, 97, 97, 99, 99, 99)))

robby <- data.table::data.table(id=dputted_dat$id, group=dputted_dat$group, y1=dputted_dat$V1, y0=JJ-dputted_dat$V1); data.table::setorder(robby, group, y1, y0); robby

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(color="blue")

dose <- dputted_dat$group

y_01 <- dputted_dat$V1

id <- dputted_dat$id

y_cbind = cbind(y_01, JJ-y_01)

cbind(id, y_cbind, dose)

gapr_marg_alpha_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_1.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1),
         pmix=c(alpha=1, scl=0.50),
         p_uppb = c(  Inf,  Inf,  1,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

roger <-
  rbind(gapr_marg_alpha_2.00,gapr_marg_alpha_free,gapr_marg_alpha_1.00)
rownames(roger) <- c("alpha set at 2.00", "alpha free range", "alpha set at 1.00")
roger

a_free <- gapr_marg_alpha_free$alpha
phi_free <- (1/sqrt(a_free)) /( (1/sqrt(a_free))^a_free + gapr_marg_alpha_free$scl^a_free)^(1/a_free)
a_2.00 <- gapr_marg_alpha_2.00$alpha
phi_2.00 <- (1/sqrt(a_2.00)) /( (1/sqrt(a_2.00))^a_2.00 + gapr_marg_alpha_2.00$scl^a_2.00)^(1/a_2.00)


bm0_free <- gapr_marg_alpha_free$a0 * phi_free
bm1_free <- gapr_marg_alpha_free$a1 * phi_free
p0_free <- stable_cdf2( bm0_free           , c(a_free, 0, 1/sqrt(a_free), 0))
p1_free <- stable_cdf2( bm0_free + bm1_free, c(a_free, 0, 1/sqrt(a_free), 0))

bm0_2.00 <- gapr_marg_alpha_2.00$a0 * phi_2.00
bm1_2.00 <- gapr_marg_alpha_2.00$a1 * phi_2.00
p0_2.00 <- stable_cdf2( bm0_2.00           , c(a_2.00, 0, 1/sqrt(a_2.00), 0))
p1_2.00 <- stable_cdf2( bm0_2.00 + bm1_2.00, c(a_2.00, 0, 1/sqrt(a_2.00), 0))

## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)

### mse: a_2.00 performs better
mean(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
mean(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
## median square error: a_free performs better
median(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
median(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
### mean absolute error: a_free performs better
mean(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
mean(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))
## median absolute error: : a_free performs better
median(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
median(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))

## refit with confidence intervals:
fu_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_2.00_ci <-
  fu_2.00$coefficients[2] + 1.96*c(-1,1)*fu_2.00$se[2]
fu_2.00_ci

## these are conditional probabilities!
stable_cdf2(fu_2.00$coefficients[1], c(2, 0, 1/sqrt(2), 0))
stable_cdf2(fu_2.00$coefficients[1] +fu_2.00$coefficients[2], c(2, 0, 1/sqrt(2), 0))
## fitted values
sort(unique(fu_2.00$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
## mean(c((robby[group==1]$y1 - 38.41672)^2, (robby[group==0]$y1 - 26.96960)^2 ))

if(fu_2.00_ci[1] <0){print("Potential WINNER:  2.00 has confidence interval that includes 0:"); fu_2.00_ci
}


fu_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_free_ci <-
  fu_free$coefficients[2] + 1.96*c(-1,1)*fu_free$se[2]
fu_free_ci

## these are conditional probabilities!
stable_cdf2(fu_free$coefficients[1], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
stable_cdf2(fu_free$coefficients[1] +fu_free$coefficients[2], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
## fitted values
sort(unique(fu_free$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
#mean(c((robby[group==1]$y1 - 34.653820)^2, (robby[group==0]$y1 - 4.591598)^2 ))


if(fu_free_ci[1] >0){print("Potential WINNER:  free has confidence interval that EXcludes 0:"); fu_free_ci
}



## next two lines same?: no.  these are the conditional probabilities * 100
sort(unique(fu_2.00$fitted.values))
sort(unique(fu_free$fitted.values))
## next two lines same?: no.  these are the marginal probabilities * 100
JJ*c(p0_2.00,p1_2.00)
JJ*c(p0_free,p1_free)
## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)


#
#
#   > sort(unique(fu_2.00$fitted.values))
# [1] 26.96960 38.41672
#   > JJ*c(p0_2.00,p1_2.00)
# [1] 37.52979 43.93764
# >
#   >
#   > sort(unique(fu_free$fitted.values))
# [1]  4.591598 34.653820
# > JJ*c(p0_free,p1_free)
# [1] 20.67212 45.95814

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(color="blue")


##############################################################
##############################################################
##############################################################
##############################################################
JJ <- 100
#
#   dput(data.table::data.table(sim_mrim_data(n1=10,n2=10,J=JJ,a0=-4,a1=3.65,v=3.42,
#                                             mrim="SSS", alpha=1, gam=1/sqrt(1.71)))[,sum(y), by=c("id","group")])
#
#
# dputted_dat <-
# structure(list(id = 1:20, group = c(1, 1, 1, 1, 1, 1, 1, 1, 1,
#                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), V1 = c(92L, 91L, 42L, 70L,
#                                                                              3L, 99L, 59L, 62L, 2L, 20L, 5L, 2L, 1L, 9L, 6L, 1L, 1L, 2L, 5L,
#                                                                              7L)))
#  dataset above gave this, with significant a1 95%CI for alpha=2 and not-significant on alpha=free
# > roger
#                          a0       a1  alpha       scl    value fevals gevals niter convcode kkt1 kkt2
# alpha set at 2.00 -1.907547 2.020999 2.0000 0.6818016 75.25485     17     51    13        0   NA   NA
# alpha free range  -4.786794 5.152611 1.2611 1.2892898 69.84298     19     81    15        0   NA   NA
# alpha set at 1.00 -8.537192 8.898735 1.0000 1.6536421 70.12491     23     70    19        0   NA   NA

## try #2:
# dput(data.table::data.table(sim_mrim_data(n1=24,n2=24,J=JJ,a0=-4,a1=3.65,v=3.42,
#                                           mrim="SSS", alpha=1, gam=1/sqrt(1.71)))[,sum(y), by=c("id","group")])
#
#
dputted_dat <-
  structure(list(id = 1:48,
                 group = c(1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 V1 = c(5L,
                        67L, 95L, 19L, 1L, 5L, 6L, 66L, 3L, 96L, 90L, 20L, 86L, 20L,
                        28L, 96L, 1L, 37L, 64L, 89L, 44L, 4L, 11L, 85L, 100L, 3L, 95L,
                        0L, 1L, 31L, 3L, 99L, 99L, 3L, 13L, 7L, 2L, 5L, 5L, 100L, 78L,
                        4L, 2L, 61L, 3L, 96L, 3L, 1L)))


# dputted_dat$V1[dputted_dat$group==1 & dputted_dat$V1 > 86] <- 99
# dputted_dat$V1[dputted_dat$group==1 & dputted_dat$V1 == 4] <- 3
# dputted_dat$V1[dputted_dat$group==1 & dputted_dat$V1 == 5] <- 3
# dputted_dat$V1[dputted_dat$group==1 & dputted_dat$V1 == 6] <- 4
#
# dputted_dat$V1[dputted_dat$group==0 & dputted_dat$V1 == 78] <- 85
# dputted_dat$V1[dputted_dat$group==0 & dputted_dat$V1 == 95] <- 86

## don't do this one makes the ses bonk: #dputted_dat$V1[dputted_dat$group==0 & dputted_dat$V1 >  95] <- 99


robby <- data.table(id=dputted_dat$id, group=dputted_dat$group, y1=dputted_dat$V1, y0=JJ-dputted_dat$V1);setorder(robby, group, y1, y0); robby
# dataset above gives this: and the 95% CI do what I want them to!
# [1] "Potential WINNER:  free has confidence interval that EXcludes 0:"
# [1] 2.512907 5.875752
# > if(fu_2.00_ci[1] <0){print("Potential WINNER:  2.00 has confidence interval that includes 0:"); fu_2.00_ci
#   + }
# [1] "Potential WINNER:  2.00 has confidence interval that includes 0:"
# [1] -0.630262  1.268619
#                           a0        a1    alpha      scl    value fevals gevals niter convcode kkt1
# alpha set at 2.00 -0.6137327 0.3191785 2.000000 1.167955 202.3304     14     46    11        0   NA
# alpha free range  -4.6950302 4.1943299 1.154963 3.074436 189.2001     43    163    29        0   NA
# alpha set at 1.00 -6.9827573 6.4302173 1.000000 4.328357 189.5539    131    175    52        1   NA

dose <- dputted_dat$group

y_01 <- dputted_dat$V1

id <- dputted_dat$id

y_cbind = cbind(y_01, JJ-y_01)

cbind(id, y_cbind, dose)

gapr_marg_alpha_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_1.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1),
         pmix=c(alpha=1, scl=0.50),
         p_uppb = c(  Inf,  Inf,  1,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

roger <-
  rbind(gapr_marg_alpha_2.00,gapr_marg_alpha_free,gapr_marg_alpha_1.00)
rownames(roger) <- c("alpha set at 2.00", "alpha free range", "alpha set at 1.00")
roger

a_free <- gapr_marg_alpha_free$alpha
phi_free <- (1/sqrt(a_free)) /( (1/sqrt(a_free))^a_free + gapr_marg_alpha_free$scl^a_free)^(1/a_free)
a_2.00 <- gapr_marg_alpha_2.00$alpha
phi_2.00 <- (1/sqrt(a_2.00)) /( (1/sqrt(a_2.00))^a_2.00 + gapr_marg_alpha_2.00$scl^a_2.00)^(1/a_2.00)


bm0_free <- gapr_marg_alpha_free$a0 * phi_free
bm1_free <- gapr_marg_alpha_free$a1 * phi_free
p0_free <- stable_cdf2( bm0_free           , c(a_free, 0, 1/sqrt(a_free), 0))
p1_free <- stable_cdf2( bm0_free + bm1_free, c(a_free, 0, 1/sqrt(a_free), 0))

bm0_2.00 <- gapr_marg_alpha_2.00$a0 * phi_2.00
bm1_2.00 <- gapr_marg_alpha_2.00$a1 * phi_2.00
p0_2.00 <- stable_cdf2( bm0_2.00           , c(a_2.00, 0, 1/sqrt(a_2.00), 0))
p1_2.00 <- stable_cdf2( bm0_2.00 + bm1_2.00, c(a_2.00, 0, 1/sqrt(a_2.00), 0))


### mse: a_2.00 performs better
mean(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
mean(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
## median square error: a_free performs better
median(c( (robby[group==1]$y1 - p1_free*JJ)^2, (robby[group==0]$y1 - p0_free*JJ)^2 ))
median(c( (robby[group==1]$y1 - p1_2.00*JJ)^2, (robby[group==0]$y1 - p0_2.00*JJ)^2 ))
### mean absolute error: a_free performs better
mean(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
mean(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))
## median absolute error: : a_free performs better
median(c( abs(robby[group==1]$y1 - p1_free*JJ), abs(robby[group==0]$y1 - p0_free*JJ) ))
median(c( abs(robby[group==1]$y1 - p1_2.00*JJ), abs(robby[group==0]$y1 - p0_2.00*JJ) ))

## refit with confidence intervals:
fu_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_2.00_ci <-
  fu_2.00$coefficients[2] + 1.96*c(-1,1)*fu_2.00$se[2]
fu_2.00_ci

## these are conditional probabilities!
stable_cdf2(fu_2.00$coefficients[1], c(2, 0, 1/sqrt(2), 0))
stable_cdf2(fu_2.00$coefficients[1] +fu_2.00$coefficients[2], c(2, 0, 1/sqrt(2), 0))
## fitted values
sort(unique(fu_2.00$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
mean(c((robby[group==1]$y1 - 38.41672)^2, (robby[group==0]$y1 - 26.96960)^2 ))

if(fu_2.00_ci[1] <0){print("Potential WINNER:  2.00 has confidence interval that includes 0:"); fu_2.00_ci
}


fu_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=FALSE,
         compute_hessian = TRUE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

fu_free_ci <-
  fu_free$coefficients[2] + 1.96*c(-1,1)*fu_free$se[2]
fu_free_ci

## these are conditional probabilities!
stable_cdf2(fu_free$coefficients[1], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
stable_cdf2(fu_free$coefficients[1] +fu_free$coefficients[2], c(fu_free$coefficients[3], 0, 1/sqrt(fu_free$coefficients[3]), 0))
## fitted values
sort(unique(fu_free$fitted.values)) ## on binomial, out of 100 scale
## MSE (conditional)
mean(c((robby[group==1]$y1 - 34.653820)^2, (robby[group==0]$y1 - 4.591598)^2 ))


if(fu_free_ci[1] >0){print("Potential WINNER:  free has confidence interval that EXcludes 0:"); fu_free_ci
}



## next two lines same?: no.  these are the conditional probabilities * 100
sort(unique(fu_2.00$fitted.values))
sort(unique(fu_free$fitted.values))
## next two lines same?: no.  these are the marginal probabilities * 100
JJ*c(p0_2.00,p1_2.00)
JJ*c(p0_free,p1_free)
## hypothetical TBAs
100*(1-p0_free/p1_free)
100*(1-p0_2.00/p1_2.00)


#
#
#   > sort(unique(fu_2.00$fitted.values))
# [1] 26.96960 38.41672
#   > JJ*c(p0_2.00,p1_2.00)
# [1] 37.52979 43.93764
# >
#   >
#   > sort(unique(fu_free$fitted.values))
# [1]  4.591598 34.653820
# > JJ*c(p0_free,p1_free)
# [1] 20.67212 45.95814

#library(ggplot2)
#library(ggbeeswarm)
#ggplot(robby) + geom_violin(aes(x=factor(group), y=y1))

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]

ggplot2::ggplot(robby,ggplot2::aes(x=factor(group), y=y1)) + ggplot2::geom_boxplot() + ggbeeswarm::geom_beeswarm(color="blue")


##############################################################
##############################################################
##############################################################
## this section shows that let alpha go free it will
## go to alpha=2 and give same results as setting alpha=2
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

gapr_marg_alpha_free <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5),
         pmix=c(alpha=1.5, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )




gapr_marg_alpha_2.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=2),
         pmix=c(alpha=2, scl=0.50),
         p_uppb = c(  Inf,  Inf,  2,  Inf   ),
         p_lowb = c( -Inf, -Inf,  2,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

gapr_marg_alpha_1.00 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2(a0 + a1*dose + rand, c(alpha, 0, 1/sqrt(alpha), 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1),
         pmix=c(alpha=1, scl=0.50),
         p_uppb = c(  Inf,  Inf,  1,  Inf   ),
         p_lowb = c( -Inf, -Inf,  1,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-scl",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )

roger <-
  rbind(gapr_marg_alpha_2.00,gapr_marg_alpha_free,gapr_marg_alpha_1.00)
rownames(roger) <- c("alpha set at 2.00", "alpha free range", "alpha set at 1.00")
roger

##############################################################
##############################################################
##############################################################
## this section gets into using nlminb vs bobyqa vs L-BFGS-B
## it uses the TSJYO (random identifier) dataset
##############################################################
##############################################################
##############################################################

sim_data_binomial <- readRDS("TSJYO_sim_data_binomial.RDS")
attach(sim_data_binomial)
y_cbind <- cbind(num_y1, num_y0)

start.time<-Sys.time()
gapr_marg_lbfgs <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2((a0 + a1*group)*(1^alpha)/phi + rand, c(alpha, 0, 1, 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5, phi=0.50),
         pmix=c(alpha=1.5, phi=0.50),
         p_uppb = c(  Inf ,  Inf, 2-1e-5,   1-1e-5),
         p_lowb = c( -Inf, -Inf,  0+1e-5,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="L-BFGS-B"
  )
end.time <- Sys.time()
time.taken <- end.time - start.time
lbfgs_time.taken <- time.taken
time.taken

start.time<-Sys.time()
gapr_marg_bobyqa <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2((a0 + a1*group)*(1^alpha)/phi + rand, c(alpha, 0, 1, 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5, phi=0.50),
         pmix=c(alpha=1.5, phi=0.50),
         p_uppb = c(  Inf ,  Inf, 2-1e-5,   1-1e-5),
         p_lowb = c( -Inf, -Inf,  0+1e-5,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="bobyqa"
  )
end.time <- Sys.time()
time.taken <- end.time - start.time
bobyqa_time.taken <- time.taken
time.taken
# Warning messages:
#   1: In (function (npt = min(n + 2L, 2L * n), rhobeg = NA, rhoend = NA,  :
#                      unused control arguments ignored
#                    2: In commonArgs(par, fn, control, environment()) :
#                      maxfun < 10 * length(par)^2 is not recommended.

start.time<-Sys.time()
gapr_marg_nlminb <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2((a0 + a1*group)*(1^alpha)/phi + rand, c(alpha, 0, 1, 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5, phi=0.50),
         pmix=c(alpha=1.5, phi=0.50),
         p_uppb = c(  Inf ,  Inf, 2-1e-5,   1-1e-5),
         p_lowb = c( -Inf, -Inf,  0+1e-5,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb"
  )
end.time <- Sys.time()
time.taken <- end.time - start.time
nlminb_time.taken <- time.taken
time.taken

rbind(gapr_marg_nlminb,gapr_marg_bobyqa,gapr_marg_lbfgs)

readRDS("TSJYO_summary_table1.RDS")[5:7,]

## now rerun bobyqa since we added new options:
start.time<-Sys.time()
gapr_marg_bobyqa2 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2((a0 + a1*group)*(1^alpha)/phi + rand, c(alpha, 0, 1, 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5, phi=0.50),
         pmix=c(alpha=1.5, phi=0.50),
         p_uppb = c(  Inf ,  Inf, 2-1e-5,   1-1e-5),
         p_lowb = c( -Inf, -Inf,  0+1e-5,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="bobyqa",
         maxfun.bobyqa=1e5,
         npt.bobyqa=6
  )
end.time <- Sys.time()
time.taken <- end.time - start.time
bobyqa2_time.taken <- time.taken
time.taken

## now rerun nlminb since we added new options
start.time<-Sys.time()
gapr_marg_nlminb2 <-
  gnlrim(y=y_cbind,
         mu=~ stable_cdf2((a0 + a1*group)*(1^alpha)/phi + rand, c(alpha, 0, 1, 0)),
         ##pmu=c( a0=bm0_start, a1=bm1_start, alpha=1.5, phi=0.50),
         pmu=c( a0=0, a1=0, alpha=1.5, phi=0.50),
         pmix=c(alpha=1.5, phi=0.50),
         p_uppb = c(  Inf ,  Inf, 2-1e-5,   1-1e-5),
         p_lowb = c( -Inf, -Inf,  1.15,  0+1e-5),
         distribution="binomial",
         nest=id,
         random="rand",
         mixture="libstableR-subgauss-phi",
         ooo=TRUE,
         compute_hessian = FALSE,
         compute_kkt = FALSE,
         trace=1,
         method="nlminb",
         abs.tol.nlminb=1e-20,#0,#1e-20, ## 1e-20,
         xf.tol.nlminb=2.2e-14, ##2.2e-14,
         x.tol.nlminb=1.5e-8, ##1.5e-8,
         rel.tol.nlminb=1e-10, ##1e-10,
  )
end.time <- Sys.time()
time.taken <- end.time - start.time
nlminb_time.taken <- time.taken
time.taken

##############################################################
##############################################################
##############################################################
## this section tests conditional vs. marginal inference
## hint:  definitely need to name pmu and pmix starting values:
## more specifically trying to get this to work
## when using libstableR-subgauss-phi
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)



## this example shows it is the exact same model -- just different parameterization
## one is condition coefficients
## the other is marginal (but assumes the scales are 1 when you specify phi -- this is hardcoded currently)
cond_PPN_phi <- gnlrim(y=y_cbind,
                       mu=~pnorm( a_cond+b_cond*dose+rand),
                       pmu=c(a_cond=0,b_cond=0),
                       pmix=c(phi=0.5),
                       p_uppb = c(Inf ,  Inf, 1-1e-5),
                       p_lowb = c(-Inf, -Inf, 0+1e-5),
                       distribution="binomial",
                       nest=id,
                       random="rand",
                       mixture="normal-phi")




marg_PPN_phi <- gnlrim(y=y_cbind,
                       mu=~pnorm( (a_marg+b_marg*dose)/phi+rand),
                       pmu=c(a_marg=0,b_marg=0,phi=0.5),
                       pmix=c(phi=0.5),
                       p_uppb = c(Inf ,  Inf,   1-1e-5),
                       p_lowb = c(-Inf, -Inf,   0+1e-5),
                       distribution="binomial",
                       nest=id,
                       random="rand",
                       mixture="normal-phi")


cond_PPN_phi$coefficients
marg_PPN_phi$coefficients

## next two lines same?
cond_PPN_phi$coefficients[1:2]*cond_PPN_phi$coefficients[3]
marg_PPN_phi$coefficients[1:2]

## next two lines same?
marg_PPN_phi$coefficients[1:2]/marg_PPN_phi$coefficients[3]
cond_PPN_phi$coefficients[1:2]


## next two lines same?
marg_PPN_phi$maxlike
cond_PPN_phi$maxlike

## next two lines same?
sort(unique(marg_PPN_phi$fitted.values))
sort(unique(cond_PPN_phi$fitted.values))


## so now we start to work in stable_cdf2
## keep in mind that the `-phi` parameterization assumes a scale of 1
## but as a Nolan parameterized stable distribution a standard normal has a scale
## of 1/sqrt(2)

cond_PPN_phi_libstableR <- gnlrim(y=y_cbind,
                                  mu=~stable_cdf2( a_cond+b_cond*dose+rand, c(alpha,0,1/sqrt(2),0)),
                                  pmu=c(a_cond=0,b_cond=0,alpha=2),
                                  pmix=c(alpha=2, phi=0.5),
                                  p_uppb = c(Inf ,  Inf, 2, 1-1e-5),
                                  p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                  distribution="binomial",
                                  nest=id,
                                  random="rand",
                                  mixture="libstableR-subgauss-phi")



## next two lines the same? No, because phi assumes scale=1 (this is hardcoded).
## and to get same beta coeffs you must assume the 1/sqrt(2)
## so not surprising that betas are the same but phi is different.
cond_PPN_phi$coefficients
cond_PPN_phi_libstableR$coefficients

cond_PPN_phi$maxlike
cond_PPN_phi_libstableR$maxlike



marg_PPN_phi_libstableR1 <- gnlrim(y=y_cbind,
                                   mu=~stable_cdf2( (a_marg+b_marg*dose)/phi+rand, c(alpha,0,1,0)),
                                   pmu=c(a_marg=0,b_marg=0,phi=0.5,alpha=2),
                                   pmix=c(alpha=2, phi=0.5),
                                   p_uppb = c(Inf ,  Inf, 2, 1-1e-5),
                                   p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                   distribution="binomial",
                                   nest=id,
                                   random="rand",
                                   mixture="libstableR-subgauss-phi")



## next two lines the same? No, because phi assumes scale=1
## and to get same beta coeffs you must assume the 1/sqrt(2)
cond_PPN_phi_libstableR1$coefficients
marg_PPN_phi_libstableR1$coefficients
## no agreement.

##If alpha is 1?
cond_CCC_phi_libstableR1 <- gnlrim(y=y_cbind,
                                   mu=~stable_cdf2( a_cond+b_cond*dose+rand, c(alpha,0,1,0)),
                                   pmu=c(a_cond=0,b_cond=0,alpha=1),
                                   pmix=c(alpha=1, phi=0.5),
                                   p_uppb = c(Inf ,  Inf, 1, 1-1e-5),
                                   p_lowb = c(-Inf, -Inf, 1, 0+1e-5),
                                   distribution="binomial",
                                   nest=id,
                                   random="rand",
                                   mixture="libstableR-subgauss-phi")

marg_CCC_phi_libstableR1 <- gnlrim(y=y_cbind,
                                   mu=~stable_cdf2( (a_marg+b_marg*dose)*((1^alpha)/phi)+rand, c(alpha,0,1,0)),
                                   pmu=c(a_marg=0,b_marg=0,alpha=1,phi=0.5),
                                   pmix=c(alpha=1, phi=0.5),
                                   p_uppb = c(Inf ,  Inf, 1, 1-1e-5),
                                   p_lowb = c(-Inf, -Inf, 1, 0+1e-5),
                                   distribution="binomial",
                                   nest=id,
                                   random="rand",
                                   mixture="libstableR-subgauss-phi",
                                   compute_hessian=TRUE,
                                   compute_kkt=FALSE,
                                   trace=1)

cond_CCC_phi_libstableR1$coefficients
marg_CCC_phi_libstableR1$coefficients
marg_CCC_phi_libstableR1$se

cond_CCC_phi_libstableR1$maxlike
marg_CCC_phi_libstableR1$maxlike



## no agreement.  I'm coming to the realization that `mixture=blah-blah-phi` should only be used when
## marginal coefficients and phi (explicitly) appear in the `mu` statement.  Furthermore, the
## parameters have to be listed correctly in `pmu` and `pmix`.  That is, the order of pmu must
## match the order the parameters appear in `mu`.  This means you might have to do a `1^alpha/phi` trick
## as opposed to `1/phi` because phi has to be last in pmix and pmu.  Note:  even the "trick" doesn't work for
## conditional models -- I put in alpha and phi explicitly:
# cond_CCC_phi_libstableR1_trick <- gnlrim(y=y_cbind,
#                                          mu=~stable_cdf2( (a_cond+b_cond*dose)*1^alpha*1^phi +rand, c(alpha,0,1,0)),
#                                          pmu=c(a_cond=0,b_cond=0,alpha=1,phi=0.5),
#                                          pmix=c(alpha=1, phi=0.5),
#                                          p_uppb = c(Inf ,  Inf, 1, 1-1e-5),
#                                          p_lowb = c(-Inf, -Inf, 1, 0+1e-5),
#                                          distribution="binomial",
#                                          nest=id,
#                                          random="rand",
#                                          mixture="libstableR-subgauss-phi")
# > cond_CCC_phi_libstableR1_trick$maxlike
# [1] 9.989502 ## should be 10.03177
# > cond_CCC_phi_libstableR1$maxlike
# [1] 9.989502 ## should be 10.03177


## observe these models that use `-scl` in a phi configuration in `mu`:

cond_PPN_var <- gnlrim(y=y_cbind,
                       mu=~pnorm( a_cond+b_cond*dose+rand),
                       pmu=c(a_cond=0,b_cond=0),
                       pmix=c(var=0.5),
                       p_uppb = c(Inf ,  Inf, Inf),
                       p_lowb = c(-Inf, -Inf, 0+1e-5),
                       distribution="binomial",
                       nest=id,
                       random="rand",
                       mixture="normal-var")




marg_PPN_var <- gnlrim(y=y_cbind,
                       mu=~pnorm( (a_marg+b_marg*dose)*sqrt(1+var)+rand),
                       pmu=c(a_marg=0,b_marg=0,var=0.5),
                       pmix=c(var=0.5),
                       p_uppb = c(Inf ,  Inf,   Inf),
                       p_lowb = c(-Inf, -Inf,   0+1e-5),
                       distribution="binomial",
                       nest=id,
                       random="rand",
                       mixture="normal-var")


cond_PPN_var$coefficients
marg_PPN_var$coefficients

## next two lines same?
cond_PPN_var$coefficients[1:2]*1/sqrt(1+cond_PPN_var$coefficients[3])
marg_PPN_var$coefficients[1:2]

## next two lines same?
marg_PPN_phi$coefficients[1:2]*sqrt(1+marg_PPN_var$coefficients[3])
cond_PPN_phi$coefficients[1:2]


## next two lines same?
marg_PPN_var$maxlike
cond_PPN_var$maxlike

## next two lines same?
sort(unique(marg_PPN_var$fitted.values))
sort(unique(cond_PPN_var$fitted.values))


## so now we start to work in stable_cdf2
## keep in mind that the `-phi` parameterization assumes a scale of 1
## but as a Nolan parameterized stable distribution a standard normal has a scale
## of 1/sqrt(2)

cond_PPN_scl_libstableR <- gnlrim(y=y_cbind,
                                  mu=~stable_cdf2( a_cond+b_cond*dose+rand, c(alpha,0,1/sqrt(2),0)),
                                  pmu=c(a_cond=0,b_cond=0,alpha=2),
                                  pmix=c(alpha=2, scl=0.5),
                                  p_uppb = c(Inf ,  Inf, 2, Inf),
                                  p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                  distribution="binomial",
                                  nest=id,
                                  random="rand",
                                  mixture="libstableR-subgauss-scl")



## next two lines the same? Yes!
c(cond_PPN_var$coefficients[1:2], 2, cond_PPN_var$coefficients[3],
  sqrt(cond_PPN_var$coefficients[3]/2))

c(cond_PPN_scl_libstableR$coefficients[1:3],
  2*cond_PPN_scl_libstableR$coefficients[4]^2 ,
  cond_PPN_scl_libstableR$coefficients[4])

## can we do marginal with libstableR?
## YES!, if we are sure to make sure to have alpha appear before scl in the `mu` statement
marg_PPN_scl_libstableR <- gnlrim(y=y_cbind,
                                  mu=~stable_cdf2( (a_marg+b_marg*dose)*sqrt((1/sqrt(2))^alpha+scl^2)/(1/sqrt(2)) +rand, c(alpha,0,1/sqrt(2),0)),
                                  pmu=c(a_marg=0,b_marg=0,alpha=2, scl=0.5),
                                  pmix=c(alpha=2, scl=0.5),
                                  p_uppb = c(Inf ,  Inf, 2, Inf),
                                  p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                  distribution="binomial",
                                  nest=id,
                                  random="rand",
                                  mixture="libstableR-subgauss-scl")

## next two lines the same? Yes!
marg_PPN_scl_libstableR$coefficients[1:2] * sqrt((1/sqrt(2))^2+marg_PPN_scl_libstableR$coefficients[4]^2)/(1/sqrt(2))

c(cond_PPN_var$coefficients[1:2], 2, cond_PPN_var$coefficients[3],
  sqrt(cond_PPN_var$coefficients[3]/2))


cond_PPN_var$maxlike
marg_PPN_var$maxlike
marg_PPN_scl_libstableR$maxlike
cond_PPN_scl_libstableR$maxlike



## see some models fit for Cauchy-Cauchy-Cauchy:

## and now for phi:
cond_fit_CCC_phi <- gnlrim(y=y_cbind,
                           mu=~pcauchy(a_cond+b_cond*dose+rand),
                           pmu=c(a_cond=0,b_cond=0),
                           pmix=c(phi=0.5),
                           p_uppb = c(Inf ,  Inf,   1-1e-5),
                           p_lowb = c(-Inf, -Inf,   0+1e-5),
                           distribution="binomial",
                           nest=id,
                           random="rand",
                           mixture="Cauchy-phi")
cond_fit_CCC_phi$coefficients
cond_fit_CCC_phi$se


marg_fit_CCC_phi <- gnlrim(y=y_cbind,
                           mu=~pcauchy((a_marg+b_marg*dose)/phi+rand),
                           pmu=c(a_marg=0,b_marg=0,phi=0.5),
                           pmix=c(phi=0.5),
                           p_uppb = c(Inf ,  Inf,   1-1e-5),
                           p_lowb = c(-Inf, -Inf,   0+1e-5),
                           distribution="binomial",
                           nest=id,
                           random="rand",
                           mixture="Cauchy-phi")
marg_fit_CCC_phi$coefficients
marg_fit_CCC_phi$se


## next two lines same? NO.  `phi` can only be used with marginal parms and must explicitly appear in mu
cond_fit_CCC_phi$coefficients[1:2] * cond_fit_CCC_phi$coefficients[3]
marg_fit_CCC_phi$coefficients[1:2]


## next two lines same? NO.  `phi` can only be used with marginal parms and must explicitly appear in mu
marg_fit_CCC_phi$coefficients[1:2] / marg_fit_CCC_phi$coefficients[3]
cond_fit_CCC_phi$coefficients[1:2]

marg_fit_CCC_phi$maxlike
cond_fit_CCC_phi$maxlike

## note the couplets above are different.
## now see `-scl` in a phi configuration (1/phi = 1+scl)

## and now for scl:
cond_fit_CCC_scl <- gnlrim(y=y_cbind,
                           mu=~pcauchy(a_cond+b_cond*dose+rand),
                           pmu=c(a_cond=0,b_cond=0),
                           pmix=c(scl=0.5),
                           p_uppb = c(Inf ,  Inf,   Inf),
                           p_lowb = c(-Inf, -Inf,   0+1e-5),
                           distribution="binomial",
                           nest=id,
                           random="rand",
                           mixture="Cauchy-scl")
cond_fit_CCC_scl$coefficients
cond_fit_CCC_scl$se


marg_fit_CCC_scl <- gnlrim(y=y_cbind,
                           mu=~pcauchy((a_marg+b_marg*dose)*(1+scl)+rand),
                           pmu=c(a_marg=0,b_marg=0,scl=0.5),
                           pmix=c(scl=0.5),
                           p_uppb = c(Inf ,  Inf,   Inf),
                           p_lowb = c(-Inf, -Inf,   0+1e-5),
                           distribution="binomial",
                           nest=id,
                           random="rand",
                           mixture="Cauchy-scl")
marg_fit_CCC_scl$coefficients
marg_fit_CCC_scl$se


## next two lines same?
cond_fit_CCC_scl$coefficients[1:2] * (1/(1+cond_fit_CCC_scl$coefficients[3] ))
marg_fit_CCC_scl$coefficients[1:2]


## next two lines same?
marg_fit_CCC_scl$coefficients[1:2] * (1+marg_fit_CCC_scl$coefficients[3])
cond_fit_CCC_scl$coefficients[1:2]


marg_fit_CCC_scl$maxlike
cond_fit_CCC_scl$maxlike



##############################################################
##############################################################
##############################################################
## this section
## get estimates but don't calculate hessian:
## To optimize with method=="nlminb" without computing the hessian,
## do so with the following settings:
## `ooo=TRUE, compute_hessian=FALSE, compute_kkt=FALSE.`
## mainly a testing feature
## speeds things up considerably
## to get SEs though, must get Hessian
## this will error if ooo=FALSE because ooo=FALSE assumes
## hessian was calculated and tries to invert to get covariance/SEs
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)


fit_PPN_libstableR_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf2(a+b*dose+rand, c(alpha, 0, 1/sqrt(2), 0)),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstableR-subgauss-scl",
                                 ooo=TRUE,
                                 compute_hessian = FALSE,
                                 compute_kkt = FALSE,
                                 trace=1
)

fit_PPN_libstableR_scl

##

dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
y <- c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
       1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
       20.972267, 17.178012)
id <- rep(1:4, each=5)

## fit with lmer
lmer_fit <- lme4::lmer(y~dose + (1|id), REML=FALSE)

## fit with gnlrim
gnlrim_fit <-
gnlrim(y,
       mu=~a+b*dose+rand,
       random="rand",
       nest=id,
       pmu=c(a=8.7,b=0.25),
       pshape = c(shape=1),
       pmix=c(var=3.0938^2),
       p_uppb = c(10,  1, 5, 3.0938^3),
       p_lowb = c( 5, -1, 0, 0)
       )

## show fits are the same:
## intercept (a) slope (b)
summary(lmer_fit)$coeff[,1]
gnlrim_fit$coeff[1:2]

## Residuals standard deviation
## sigma_epsilon = 5.58
summary(lmer_fit)$varcor
sqrt(exp(gnlrim_fit$coeff[3]))

## random effects standard deviation
## sigma_id = 3.0938
summary(lmer_fit)$varcor
sqrt(gnlrim_fit$coeff[4])

## likelihood
summary(lmer_fit)$logLik
-gnlrim_fit$maxlike

## Take same model but hold
## random effects standard deviation
## sigma_id   :=  9
## sigma^2_id := 81
gnlrim_fit2 <-
  gnlrim(y,
         mu=~a+b*dose+rand,
         random="rand",
         nest=id,
         pmu=c(a=8.7,b=0.25),
         pshape = c(shape=1),
         pmix=c(var=9^2),
         p_uppb = c(10,  1, 5, 9^2),
         p_lowb = c( 5, -1, 0, 9^2)
  )

gnlrim_fit2$coeff
gnlrim_fit2$se


##############################################################
##############################################################
##############################################################
## this section tests the libstableR implementation of
## subgauss and compares it to stabledist implementation of
## subgauss and compares to PPN with [dp]norm
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

fit_PPN_scl <- gnlrim(y=y_cbind,
                      mu=~pnorm(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(scl=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="normal-var")
fit_PPN_scl$coefficients
fit_PPN_scl$se


#Now try reproducing PPN with *!*stabledist*!* subgauss:

fit_PPN_stabledist_scl <- gnlrim(y=y_cbind,
                                 mu=~pstable2(a+b*dose+rand, alpha, 0, 1/sqrt(2), 0),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="stabledist-subgauss-scl"#,
                                 #ooo=TRUE,
                                 #trace=1
                                 )
fit_PPN_stabledist_scl$coefficients
fit_PPN_stabledist_scl$se

#Now try reproducing PPN with **libstableR** subgauss:

fit_PPN_libstableR_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf2(a+b*dose+rand, c(alpha, 0, 1/sqrt(2), 0)),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstableR-subgauss-scl"#,
                                 # ooo=TRUE,
                                 # compute_hessian = FALSE,
                                 # compute_kkt = FALSE,
                                 # trace=1
)

## should be 0s, not NAs.
attr(fit_PPN_libstableR_scl, "details")[1,]$nhatend
bonk <- attr(fit_PPN_libstableR_scl, "details")[1,]$nhatend
bonk_no_na <- matrix(bonk[!is.na(bonk)], nrow=3)
bonk_no_na_cov <- solve(bonk_no_na)
bonk_se <- sqrt(diag(bonk_no_na_cov))
## compared these two:
bonk_se
fit_PPN_scl$se



fit_PPN_libstableR_scl$coefficients
fit_PPN_libstableR_scl$se


# next two function are now in the package gnlrim
# ## what if we made a wrapper for libstableR
# ## that would take out of bound pars and force them in bounds
# ## ***internally***
# stable_cdf2 <- function(x, pars, parameterization=0L, tol=1e-12){
#   libstableR::stable_cdf(x,
#                          c(min(max(pars[1],0+1e-20),2),
#                            min(max(pars[2],-1),1),
#                            min(max(pars[3],0+1e-20),Inf),
#                            pars[4]),
#                          parameterization,
#                          tol
#   )
# }
#
# stable_pdf2 <- function(x, pars, parameterization=0L, tol=1e-12){
#   libstableR::stable_pdf(x,
#                          c(min(max(pars[1],0+1e-20),2),
#                            min(max(pars[2],-1),1),
#                            min(max(pars[3],0+1e-20),Inf),
#                            pars[4]),
#                          parameterization,
#                          tol
#   )
# }
#

## alpha <= 2 is forced internally
stable_cdf2(1, c( 2.0,0,1,0))
stable_cdf2(1, c( 2.1,0,1,0))
stable_cdf2(1, c(22.1,0,1,0))

## alpha >= 0+1e-20 is forced internally
stable_cdf2(1, c( 0+1e-20,0,1,0))
stable_cdf2(1, c( 0,0,1,0))
stable_cdf2(1, c(22.1,0,1,0))

## beta >= -1 is forced internally
stable_cdf2(1, c( 1,-1.0,1,0))
stable_cdf2(1, c( 1,-1.1,1,0))
stable_cdf2(1, c( 1,-2,1,0))


fit2_PPN_libstableR_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf2(a+b*dose+rand, c(alpha, 0, 1/sqrt(2), 0)),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstableR-subgauss-scl"#,
                                 #ooo=TRUE,
                                 #trace=1
)

## should be 0s, not NAs.
attr(fit2_PPN_libstableR_scl, "details")[1,]$nhatend
bonk <- attr(fit2_PPN_libstableR_scl, "details")[1,]$nhatend
bonk_no_na <- matrix(bonk[!is.na(bonk)], nrow=3)
bonk_no_na_cov <- solve(bonk_no_na)
bonk_se <- sqrt(diag(bonk_no_na_cov))
## compared these two:
bonk_se
fit_PPN_scl$se

fit2_PPN_libstableR_scl$coefficients
fit2_PPN_libstableR_scl$se

fit_PPN_scl$coefficients
fit_PPN_scl$se




## all the same?
fit_PPN_scl$coefficients
fit_PPN_stabledist_scl$coefficients
fit2_PPN_libstableR_scl$coefficients

fit_PPN_scl$se
fit_PPN_stabledist_scl$se
fit2_PPN_libstableR_scl$se


##############################################################
##############################################################
##############################################################
## this section tests the libstableR implementation of
## subgauss and compares it to stabledist implementation of
## subgauss and compares to CCC with [dp]cauchy
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

fit_CCC_scl <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(scl=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-scl")
fit_CCC_scl$coefficients
fit_CCC_scl$se


#Now try reproducing CCC with *!*stabledist*!* subgauss:

fit_CCC_stabledist_scl <- gnlrim(y=y_cbind,
                                 mu=~pstable(a+b*dose+rand, alpha, 0, 1, 0),
                                 pmu=c(a=0,b=0, alpha=1),
                                 pmix=c(alpha=1, scl=1),
                                 p_uppb = c(Inf ,  Inf, 1, Inf),
                                 p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="stabledist-subgauss-scl")
fit_CCC_stabledist_scl$coefficients
fit_CCC_stabledist_scl$se

#Now try reproducing CCC with **libstableR** subgauss:

fit_CCC_libstableR_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf(a+b*dose+rand, c(alpha, 0, 1, 0)),
                                 pmu=c(a=0,b=0, alpha=1),
                                 pmix=c(alpha=1, scl=1),
                                 p_uppb = c(Inf ,  Inf, 1, Inf),
                                 p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstableR-subgauss-scl")
fit_CCC_libstableR_scl$coefficients
fit_CCC_libstableR_scl$se


## all the same?
fit_CCC_scl$coefficients
fit_CCC_stabledist_scl$coefficients
fit_CCC_libstableR_scl$coefficients

fit_CCC_scl$se
fit_CCC_stabledist_scl$se
fit_CCC_libstableR_scl$se



## and now for phi:
fit_CCC_phi <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(phi=0.5),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-phi")
fit_CCC_phi$coefficients
fit_CCC_phi$se


#Now try reproducing CCC with *!*stabledist*!* subgauss:

fit_CCC_stabledist_phi <- gnlrim(y=y_cbind,
                                 mu=~pstable(a+b*dose+rand, alpha, 0, 1, 0),
                                 pmu=c(a=0,b=0, alpha=1),
                                 pmix=c(alpha=1, phi=0.5),
                                 p_uppb = c(Inf ,  Inf, 1,   1-1e-5),
                                 p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="stabledist-subgauss-phi")
fit_CCC_stabledist_phi$coefficients
fit_CCC_stabledist_phi$se

#Now try reproducing CCC with **libstableR** subgauss:

fit_CCC_libstableR_phi <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf(a+b*dose+rand, c(alpha, 0, 1, 0)),
                                 pmu=c(a=0,b=0, alpha=1),
                                 pmix=c(alpha=1, phi=0.5),
                                 p_uppb = c(Inf ,  Inf, 1,   1-1e-5),
                                 p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstableR-subgauss-phi")
fit_CCC_libstableR_phi$coefficients
fit_CCC_libstableR_phi$se


## all the same?
fit_CCC_phi$coefficients
fit_CCC_stabledist_phi$coefficients
fit_CCC_libstableR_phi$coefficients

fit_CCC_phi$se
fit_CCC_stabledist_phi$se
fit_CCC_libstableR_phi$se



##############################################################
##############################################################
##############################################################
## this section tests the stabledist implementation of subgauss
## and compares to CCC
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

fit_CCC_scl <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(scl=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-scl")

fit_CCC_phi <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=c(a=fit_CCC_scl$coefficients[1],
                            b=fit_CCC_scl$coefficients[2]),
                      pmix=c(phi=0.25),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-phi")

fit_CCC_scl$coefficients
fit_CCC_phi$coefficients

fit_CCC_scl$se
fit_CCC_phi$se

1  /( 1 + fit_CCC_scl$coefficients[3])

#Now try reproducing CCC with subguass:

fit_CCC_subg_scl <- gnlrim(y=y_cbind,
                      mu=~pstable(a+b*dose+rand, alpha, 0, 1, 0),
                      pmu=c(a=0,b=0, alpha=1),
                      pmix=c(alpha=1, scl=1),
                      p_uppb = c(Inf ,  Inf, 1, Inf),
                      p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="stabledist-subgauss-scl")

fit_CCC_subg_phi <- gnlrim(y=y_cbind,
                      mu=~pstable(a+b*dose+rand, alpha, 0, 1, 0),
                      pmu=c(a=fit_CCC_scl$coefficients[1],
                            b=fit_CCC_scl$coefficients[2],
                            alpha=1),
                      pmix=c(alpha=1,phi=0.25),
                      p_uppb = c(Inf ,  Inf,  1, 1-1e-5),
                      p_lowb = c(-Inf, -Inf,  1, 0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="stabledist-subgauss-phi")

fit_CCC_subg_scl$coefficients
fit_CCC_subg_phi$coefficients

fit_CCC_subg_scl$se
fit_CCC_subg_phi$se

1  /( 1 + fit_CCC_subg_scl$coefficients[3])


fit_CCC_scl$coefficients
fit_CCC_phi$coefficients

fit_CCC_scl$se
fit_CCC_phi$se

## hessian?
fit_CCC_subg_scl_hess <- gnlrim(y=y_cbind,
                           mu=~pstable(a+b*dose+rand, alpha, 0, 1, 0),
                           pmu=c(a=0,b=0, alpha=1),
                           pmix=c(alpha=1, scl=1),
                           p_uppb = c(Inf ,  Inf, 1, Inf),
                           p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                           distribution="binomial",
                           nest=id,
                           random="rand",
                           mixture="stabledist-subgauss-scl",
                           ooo=TRUE)

## with `ooo=TRUE` hessian extraction (rows or by name):
hesser <- attr(fit_CCC_subg_scl_hess, "details")[1,]$nhatend
hesser
## use number p_uppb == p_lowb to identify how many held constant
## also need to identify which col(s) and row(s) are constant
## because we need to rebuild after inverting
hesser_no_0 <- matrix(hesser[hesser != 0], ncol=3,nrow=3)
hesser_no_0
cov_no_0 <- solve(hesser_no_0)
cov_no_0
se_no_0 <- sqrt(diag(cov_no_0))
se_no_0
fit_CCC_scl$se

## below is future LLB vs subgauss...

## now try that a=1.89 gam_d = 1.2 = gam_h but wait phi assumes 1

fit_LLB_var <- gnlrim(y=y_cbind,
                      mu=~plogis(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(var=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="logit-bridge-var")

fit_LLB_phi <- gnlrim(y=y_cbind,
                      mu=~plogis(a+b*dose+rand),
                      pmu=c(a=fit_LLB_var$coefficients[1],
                            b=fit_LLB_var$coefficients[2]),
                      pmix=c(phi=0.25),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="logit-bridge-phi")

fit_LLB_var$coefficients
fit_LLB_phi$coefficients

fit_LLB_var$se
fit_LLB_phi$se

sqrt( pi^2/3  /( pi^2/3 + fit_LLB_var$coefficients[3]))





##############################################################
##############################################################
##############################################################
## this section tests conditional vs. marginal inference
## hint:  definitely need to name pmu and pmix starting values:
##############################################################
##############################################################
##############################################################

dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
                         1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
                         11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
                         20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

cond_PPN_phi <- gnlrim(y=y_cbind,
                       mu=~pnorm( a_cond+b_cond*dose+rand),
                       pmu=c(a_cond=0,b_cond=0),
                       pmix=c(phi=0.5),
                       p_uppb = c(Inf ,  Inf,   1-1e-5),
                       p_lowb = c(-Inf, -Inf,   0+1e-5),
                       distribution="binomial",
                       nest=id,
                       random="rand",
                       mixture="normal-phi")




marg_PPN_phi <- gnlrim(y=y_cbind,
                      mu=~pnorm( (a_marg+b_marg*dose)/phi+rand),
                      pmu=c(a_marg=0,b_marg=0,phi=0.5),
                      pmix=c(phi=0.5),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="normal-phi")


cond_PPN_phi$coefficients
marg_PPN_phi$coefficients

## next two lines same?
cond_PPN_phi$coefficients[1:2]*cond_PPN_phi$coefficients[3]
marg_PPN_phi$coefficients[1:2]

## next two lines same?
marg_PPN_phi$coefficients[1:2]/marg_PPN_phi$coefficients[3]
cond_PPN_phi$coefficients[1:2]


## next two lines same?
marg_PPN_phi$maxlike
cond_PPN_phi$maxlike

## next two lines same?
sort(unique(marg_PPN_phi$fitted.values))
sort(unique(cond_PPN_phi$fitted.values))







dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
       1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
       20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

fit_PPN_var <- gnlrim(y=y_cbind,
                      mu=~pnorm(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(var=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="normal-var")

fit_PPN_phi <- gnlrim(y=y_cbind,
                      mu=~pnorm(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(phi=0.5),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="normal-phi")

fit_PPN_var$coefficients
fit_PPN_phi$coefficients

fit_PPN_var$se
fit_PPN_phi$se

tau2 <- fit_PPN_var$coefficients[3]
1/sqrt(1+tau2)


fit_LLB_var <- gnlrim(y=y_cbind,
                      mu=~plogis(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(var=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="logit-bridge-var")

fit_LLB_phi <- gnlrim(y=y_cbind,
                      mu=~plogis(a+b*dose+rand),
                      pmu=c(a=fit_LLB_var$coefficients[1],
                            b=fit_LLB_var$coefficients[2]),
                      pmix=c(phi=0.25),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="logit-bridge-phi")

fit_LLB_var$coefficients
fit_LLB_phi$coefficients

fit_LLB_var$se
fit_LLB_phi$se

sqrt( pi^2/3  /( pi^2/3 + fit_LLB_var$coefficients[3]))



fit_CCC_scl <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(scl=1),
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-scl")

fit_CCC_phi <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=c(a=fit_CCC_scl$coefficients[1],
                            b=fit_CCC_scl$coefficients[2]),
                      pmix=c(phi=0.25),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-phi")

fit_CCC_scl$coefficients
fit_CCC_phi$coefficients

fit_CCC_scl$se
fit_CCC_phi$se

1  /( 1 + fit_CCC_scl$coefficients[3])




















dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
y <- c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
       1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
       20.972267, 17.178012)
id <- rep(1:4, each=5)


fit_nlminb <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(a=8.7,b=0.25),
                            pshape=3.44, pmix=2.3, trace=1)
fit_nlminb$coefficients

fit_nlminb_bounds <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(a=8.7,b=0.25),
                     pshape=3.44, pmix=8, trace=1, p_lowb=c(8,.2,3,8), p_uppb=c(9,.3,4,10))
fit_nlminb_bounds$coefficients



fit_nlm <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(a=8.7,b=0.25),
                  pshape=3.44, pmix=2.3, method="nlm")

fit_both <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(a=8.7,b=0.25),
                  pshape=3.44, pmix=2.3, method=c("nlminb","nlm"), ooo=TRUE)
fit_both
## with `ooo=TRUE` hessian extraction (rows or by name):
attr(fit_both, "details")[1,]$nhatend
attr(fit_both, "details")[2,]$nhatend
attr(fit_both, "details")["nlminb",]$nhatend
attr(fit_both, "details")["nlm",]$nhatend


fit_both[1, 1:attr(fit_both, "npar")]



fit_both <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(a=8.7,b=0.25),
                   pshape=3.44, pmix=c(var=2.3), method=c("nlminb","nlm"))
fit_both
rmutil::print.gnlm(fit_both)
fit_both


fit_gnlmix <- gnlmix(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(8.7,0.25),
                            pshape=3.44, pmix=2.3)
## hessian extraction (1st row or by name):
rmutil::print.gnlm(fit_gnlmix)

fit_nlm$fitted.values
fit_gnlmix$fitted.values

fit_nlm$residuals
fit_gnlmix$residuals
