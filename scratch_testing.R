########################################$$$$#########
## 2023-10-11A                                     ##
## START two random parameters: run on cluster    ###
## bivariate cauchy marg/cond                      ##
#####################################################
# (rand.int.rand.slopes.nonzero.corr.CUBA.cauchy <-
#    gnlrim::gnlrem(y=ybind,
#                   mu = ~ pcauchy(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
#
#                   pmu = c(Intercept=-1.2, b_p=1),
#                   pmix=c(var1=1, var2=1.3, corr12= 0.30),
#
#                   p_uppb = c(  0,   4, 1.0, 4.00, 4.00, 0.90),
#                   p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
#                   distribution="binomial",
#                   nest=id,
#                   random=c("rand1", "rand2"),
#                   mixture="bivariate-cauchy-corr",
#                   ooo=TRUE,
#                   compute_hessian = FALSE,
#                   compute_kkt = FALSE,
#                   trace=1,
#                   method='nlminb',
#                   int2dmethod="cuba",
#                   tol.pcubature = 0.1,
#                   abs.tol.nlminb = 1e-2,
#                   xf.tol.nlminb =  1e-2,
#                   x.tol.nlminb =   1e-2,
#                   rel.tol.nlminb = 1e-2
#    )
# )
# ## ran this and recorded the results locally (see far below).  change tol.pcubature
# ## and run a few more times on cluster:
#
# (rand.int.rand.slopes.nonzero.corr.CUBA.cauchy <-
#     gnlrim::gnlrem(y=ybind,
#                    mu = ~ pcauchy(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
#
#                    pmu = c(Intercept=-1.2, b_p=1),
#                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#
#                    p_uppb = c(  0,   4, 4.00, 4.00, 0.90),
#                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                    distribution="binomial",
#                    nest=id,
#                    random=c("rand1", "rand2"),
#                    mixture="bivariate-cauchy-corr",
#                    ooo=TRUE,
#                    compute_hessian = FALSE,
#                    compute_kkt = FALSE,
#                    trace=1,
#                    method='nlminb',
#                    int2dmethod="cuba",
#                    tol.pcubature = 0.5,
#                    abs.tol.nlminb = 1e-2,
#                    xf.tol.nlminb =  1e-2,
#                    x.tol.nlminb =   1e-2,
#                    rel.tol.nlminb = 1e-2
#     )
# )
# Assemble the answers
#        Intercept      b_p      var1     var2    corr12    value fevals gevals
# nlminb -1.970935 0.995794 0.7449743 1.677083 0.6275397 3671.545      7     31
# niter convcode kkt1 kkt2    xtime
# nlminb     5        0   NA   NA 9763.134

## change tol.pcubature=0.5 -> 0.09
# (rand.int.rand.slopes.nonzero.corr.CUBA.cauchy <-
#     gnlrim::gnlrem(y=ybind,
#                    mu = ~ pcauchy(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
#
#                    pmu = c(Intercept=-1.2, b_p=1),
#                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#
#                    p_uppb = c(  0,   4, 4.00, 4.00, 0.90),
#                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                    distribution="binomial",
#                    nest=id,
#                    random=c("rand1", "rand2"),
#                    mixture="bivariate-cauchy-corr",
#                    ooo=TRUE,
#                    compute_hessian = FALSE,
#                    compute_kkt = FALSE,
#                    trace=1,
#                    method='nlminb',
#                    int2dmethod="cuba",
#                    tol.pcubature = 0.09,
#                    abs.tol.nlminb = 1e-2,
#                    xf.tol.nlminb =  1e-2,
#                    x.tol.nlminb =   1e-2,
#                    rel.tol.nlminb = 1e-2
#     )
# )
# Assemble the answers
#        Intercept      b_p      var1     var2    corr12    value fevals gevals
# nlminb -2.063432 1.032551 0.7800393 1.453341 0.5415967 3672.277      4     15
# niter convcode kkt1 kkt2    xtime
# nlminb     3        0   NA   NA 2457.624

## now do marginal parms for the same 0.5 and 0.09 tol.pcubatures, get parameters for graph
## run on cluster
# (marginal.CUBA.cauchy <-
#     gnlrim::gnlrem(y=ybind,
#                    mu = ~ pcauchy( (Intercept + period_numeric*b_p) *
#                                      (1.000^1.00 + (var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric)^(1.00/2) )^(1/1.00) / (1.000)  +
#                                      rand1 + period_numeric*rand2),
#
#                    pmu = c(Intercept=-1.2, b_p=1, var1=1, var2=1.3, corr12= 0.30),
#                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#
#                    p_uppb = c(  0,   4, 4.00, 4.00, 0.90),
#                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                    distribution="binomial",
#                    nest=id,
#                    random=c("rand1", "rand2"),
#                    mixture="bivariate-cauchy-corr",
#                    ooo=TRUE,
#                    compute_hessian = FALSE,
#                    compute_kkt = FALSE,
#                    trace=1,
#                    method='nlminb',
#                    int2dmethod="cuba",
#                    tol.pcubature = 0.5,
#                    abs.tol.nlminb = 1e-2,
#                    xf.tol.nlminb =  1e-2,
#                    x.tol.nlminb =   1e-2,
#                    rel.tol.nlminb = 1e-2
#     )
# )
## output:
#                           > (marginal.CUBA.cauchy <-
#                                +     gnlrim::gnlrem(y=ybind,
#                                                     +                    mu = ~ pcauchy( (Intercept + period_numeric*b_p) *
#                                                                                            +                                     (1.000^1.00 + (var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric)^(1.00/2) )^(1/1.00) / (1.000)  +
#                                                                                            +                                     rand1 + period_numeric*rand2),
#                                                     +
#                                                       +                    pmu = c(Intercept=-1.2, b_p=1, var1=1, var2=1.3, corr12= 0.30),
#                                                     +                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#                                                     +
#                                                       +                    p_uppb = c(  0,   4, 4.00, 4.00, 0.90),
#                                                     +                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                                                     +                    distribution="binomial",
#                                                     +                    nest=id,
#                                                     +                    random=c("rand1", "rand2"),
#                                                     +                    mixture="bivariate-cauchy-corr",
#                                                     +                    ooo=TRUE,
#                                                     +                    compute_hessian = FALSE,
#                                                     +                    compute_kkt = FALSE,
#                                                     +                    trace=1,
#                                                     +                    method='nlminb',
#                                                     +                    int2dmethod="cuba",
#                                                     +                    tol.pcubature = 0.5,
#                                                     +                    abs.tol.nlminb = 1e-2,
#                                                     +                    xf.tol.nlminb =  1e-2,
#                                                     +                    x.tol.nlminb =   1e-2,
#                                                     +                    rel.tol.nlminb = 1e-2
#                                                     +     )
#                              + )
#                           fn is  fn1
#                           Looking for method =  nlminb
#                           Function has  5  arguments
#                           par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
#                           par[ 2 ]:  -2   <? 1   <? 4     In Bounds
#                           par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
#                           par[ 4 ]:  0.05   <? 1.3   <? 4     In Bounds
#                           par[ 5 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
#                           [1] "2023-10-11 09:59:35.818407 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:01:22.196008 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           Analytic gradient not made available.
#                           Analytic Hessian not made available.
#                           Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.5228787
#                           Method:  nlminb
#                           [1] "2023-10-11 10:01:22.31376 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:03:08.29424 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           [1] "2023-10-11 10:03:08.294467 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:04:54.088164 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           [1] "2023-10-11 10:04:54.088384 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:06:40.152978 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           [1] "2023-10-11 10:06:40.153198 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:08:26.031108 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           [1] "2023-10-11 10:08:26.031319 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:10:11.746651 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           [1] "2023-10-11 10:10:11.746875 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:11:57.636014 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.16513"
#                           0:     3681.1651: -1.20000  1.00000  1.00000  1.30000 0.300000
#                           [1] "2023-10-11 10:11:57.636257 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:13:44.210635 ... ending   pcubature -- tol=0.5 -- ret.val is: 3715.29924"
#                           [1] "2023-10-11 10:13:44.210861 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:14:09.245176 ... ending   pcubature -- tol=0.5 -- ret.val is: 3679.93606"
#                           [1] "2023-10-11 10:14:09.245383 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:14:34.419124 ... ending   pcubature -- tol=0.5 -- ret.val is: 3679.93572"
#                           [1] "2023-10-11 10:14:34.419329 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:14:59.608576 ... ending   pcubature -- tol=0.5 -- ret.val is: 3679.93507"
#                           [1] "2023-10-11 10:14:59.608788 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:15:24.700667 ... ending   pcubature -- tol=0.5 -- ret.val is: 3679.9363"
#                           [1] "2023-10-11 10:15:24.700876 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:15:49.68237 ... ending   pcubature -- tol=0.5 -- ret.val is: 3679.9363"
#                           [1] "2023-10-11 10:15:49.682591 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:16:14.851552 ... ending   pcubature -- tol=0.5 -- ret.val is: 3679.9368"
#                           1:     3679.9361: -1.07296  1.00420 0.968245  1.30025 0.308327
#                           [1] "2023-10-11 10:16:14.851787 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:17:59.945832 ... ending   pcubature -- tol=0.5 -- ret.val is: 3674.96344"
#                           [1] "2023-10-11 10:17:59.946042 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:25:10.437264 ... ending   pcubature -- tol=0.5 -- ret.val is: 3676.03472"
#                           [1] "2023-10-11 10:25:10.437527 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:26:54.219828 ... ending   pcubature -- tol=0.5 -- ret.val is: 3674.96336"
#                           [1] "2023-10-11 10:26:54.22003 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:28:37.979184 ... ending   pcubature -- tol=0.5 -- ret.val is: 3674.96342"
#                           [1] "2023-10-11 10:28:37.979389 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:30:21.665556 ... ending   pcubature -- tol=0.5 -- ret.val is: 3674.96357"
#                           [1] "2023-10-11 10:30:21.665747 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:32:05.539072 ... ending   pcubature -- tol=0.5 -- ret.val is: 3674.96351"
#                           [1] "2023-10-11 10:32:05.539284 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:33:49.656699 ... ending   pcubature -- tol=0.5 -- ret.val is: 3674.96359"
#                           2:     3674.9634: -1.10601 0.906439 0.943690  1.32457 0.381745
#                           [1] "2023-10-11 10:33:49.656923 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:34:16.118491 ... ending   pcubature -- tol=0.5 -- ret.val is: 3681.85143"
#                           [1] "2023-10-11 10:34:16.1187 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:34:40.850864 ... ending   pcubature -- tol=0.5 -- ret.val is: 3677.27722"
#                           [1] "2023-10-11 10:34:40.851061 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:36:25.69031 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.75514"
#                           [1] "2023-10-11 10:36:25.690536 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:38:10.768688 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.75515"
#                           [1] "2023-10-11 10:38:10.768896 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:39:55.711218 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.75511"
#                           [1] "2023-10-11 10:39:55.711422 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:41:40.587735 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.75522"
#                           [1] "2023-10-11 10:41:40.587938 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:43:25.550888 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.75525"
#                           [1] "2023-10-11 10:43:25.5511 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:45:10.411189 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.75529"
#                           3:     3673.7551: -1.05387 0.909284 0.928322  1.32763 0.392125
#                           [1] "2023-10-11 10:45:10.411428 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:46:55.513955 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.15816"
#                           [1] "2023-10-11 10:46:55.514173 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:47:20.292275 ... ending   pcubature -- tol=0.5 -- ret.val is: 3676.84165"
#                           [1] "2023-10-11 10:47:20.29248 ... starting pcubature for bivariate normal or cauchy"
#                           [1] "2023-10-11 10:49:05.050862 ... ending   pcubature -- tol=0.5 --[1] "2023-10-11 10:49:05.051068 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 10:50:49.948548 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.1582"
# [1] "2023-10-11 10:50:49.948755 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 10:52:34.816716 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.1582"
# [1] "2023-10-11 10:52:34.816932 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 10:54:19.350091 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.15831"
# [1] "2023-10-11 10:54:19.350304 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 10:56:04.111473 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.1584"
# 4:     3673.1582: -1.03591 0.885063 0.907852  1.34176 0.431520
# [1] "2023-10-11 10:56:04.111714 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 10:57:50.245139 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.32919"
# 5:     3673.1582: -1.03591 0.885063 0.907852  1.34176 0.431520
# [1] "2023-10-11 10:57:50.245365 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 10:59:35.159869 ... ending   pcubature -- tol=0.5 -- ret.val is: 3673.15816"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -1.0359081  0.8850626  0.9078518  1.3417596  0.4315196
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 3673.158
#
# $fevals
# function
# 11
#
# $gevals
# gradient
# 25
#
# $nitns
# [1] 5
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 3471.885
#
# Assemble the answers
#        Intercept       b_p      var1    var2    corr12    value fevals gevals
# nlminb -1.035908 0.8850626 0.9078518 1.34176 0.4315196 3673.158     11     25
# niter convcode kkt1 kkt2    xtime
# nlminb     5        0   NA   NA 3471.885 ret.val is: 3673.15815"
#

## ran this farther below.  see farther below for output.
##glm(cbind(r, n_r) ~ x1, summed_binom_dat, family=binomial(link="cauchit"))

## then *Locally* run this part after inputting the values
## some things are commented out but it should make the graph
## at the end
phi_x_deluxe <- function(x, var1=0.9078518, var2=1.34176, corr12=0.4315196, gamma_d = 1, gamma_h=1, a=1.00 ){

  (gamma_d)  / (gamma_h^a + (var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x)^(a/2) )^(1/a)

}

datavalueset <- c(0.1, 0.12, 0.2, 0.22, 0.3, 0.32, 0.4, 0.42)
xvalueset <- c(0, datavalueset, 1, 10)
xdomain <- seq(-2,4,0.05)

plot(NA, NA,col="grey", ylim=c(0,1),xlim=range(xdomain))
abline(h=0.50, lty=3, col="lightblue")
abline(v=min(datavalueset), lty=2, col="lightgray")
abline(v=max(datavalueset), lty=2, col="lightgray")

lines(xdomain, pcauchy( -1.970935+0.995794*xdomain), ylim=c(0,1), lty=3)
text(3.6,0.68, "(b0c+b1c*x)       ", col="black")

lines(xdomain, pcauchy( -1.035908+0.8850626*xdomain), ylim=c(0,1),col="orange",lty=3)
text(3.6,0.95, "(b0m+b1m*x)       ", col="orange")

lines(xdomain, phi_x_deluxe(xdomain), ylim=c(0,1), lty=2, col='cyan')
text(3.6,0.19, "            phi(x)", col="cyan")

lines(xdomain, pcauchy((-1.970935+0.995794*xdomain)*phi_x_deluxe(xdomain)), col="red",lty=3)
text(3.6,0.52, "(b0c+b1c*x)*phi(x)", col="red")


# lines(xdomain, pcauchy( -0.72507+0.54384*xdomain), ylim=c(0,1), col="darkgreen",lty=3)
# text(3.6,0.78, "glm(b'0m+b'1m*x)       ", col="darkgreen")


## phi ranges over the dataset:
phi_x_deluxe(min(datavalueset))
phi_x_deluxe(max(datavalueset))
abline(h=phi_x_deluxe(min(datavalueset)), lty=2, col="lightpink")
abline(h=phi_x_deluxe(max(datavalueset)), lty=2, col="lightpink")




## change 0.5 to 0.09; run on cluster; then make graphic for the 0.09 case marg/cond.
# (marginal.CUBA.cauchy <-
#     gnlrim::gnlrem(y=ybind,
#                    mu = ~ pcauchy( (Intercept + period_numeric*b_p) *
#                                      (1.000^1.00 + (var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric)^(1.00/2) )^(1/1.00) / (1.000)  +
#                                      rand1 + period_numeric*rand2),
#
#                    pmu = c(Intercept=-1.2, b_p=1, var1=1, var2=1.3, corr12= 0.30),
#                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#
#                    p_uppb = c(  0,   4, 4.00, 4.00, 0.90),
#                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                    distribution="binomial",
#                    nest=id,
#                    random=c("rand1", "rand2"),
#                    mixture="bivariate-cauchy-corr",
#                    ooo=TRUE,
#                    compute_hessian = FALSE,
#                    compute_kkt = FALSE,
#                    trace=1,
#                    method='nlminb',
#                    int2dmethod="cuba",
#                    tol.pcubature = 0.09,
#                    abs.tol.nlminb = 1e-2,
#                    xf.tol.nlminb =  1e-2,
#                    x.tol.nlminb =   1e-2,
#                    rel.tol.nlminb = 1e-2
#     )
# )
## output:
# > (marginal.CUBA.cauchy <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pcauchy( (Intercept + period_numeric*b_p) *
#                                                                  +                                      (1.000^1.00 + (var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric)^(1.00/2) )^(1/1.00) / (1.000)  +
#                                                                  +                                      rand1 + period_numeric*rand2),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, var1=1, var2=1.3, corr12= 0.30),
#                           +                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   4, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-cauchy-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.09,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 4     In Bounds
# par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 4 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-10-11 11:12:07.634038 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:13:51.242385 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.5228787
# Method:  nlminb
# [1] "2023-10-11 11:13:51.356963 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:15:35.11645 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# [1] "2023-10-11 11:15:35.116672 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:17:19.210105 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# [1] "2023-10-11 11:17:19.210313 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:19:03.211856 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# [1] "2023-10-11 11:19:03.212062 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:20:46.897993 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# [1] "2023-10-11 11:20:46.898191 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:22:31.60603 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# [1] "2023-10-11 11:22:31.606238 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:24:15.931335 ... ending   pcubature -- tol=0.09 -- ret.val is: 3681.16513"
# 0:     3681.1651: -1.20000  1.00000  1.00000  1.30000 0.300000
# [1] "2023-10-11 11:24:15.931557 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:26:02.226034 ... ending   pcubature -- tol=0.09 -- ret.val is: 3715.29924"
# [1] "2023-10-11 11:26:02.226264 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:27:48.340861 ... ending   pcubature -- tol=0.09 -- ret.val is: 3675.93527"
# [1] "2023-10-11 11:27:48.341069 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:29:34.63026 ... ending   pcubature -- tol=0.09 -- ret.val is: 3675.93535"
# [1] "2023-10-11 11:29:34.630479 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:31:20.785568 ... ending   pcubature -- tol=0.09 -- ret.val is: 3675.93438"
# [1] "2023-10-11 11:31:20.785767 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:33:06.95763 ... ending   pcubature -- tol=0.09 -- ret.val is: 3675.93559"
# [1] "2023-10-11 11:33:06.957835 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:34:53.189072 ... ending   pcubature -- tol=0.09 -- ret.val is: 3675.93545"
# [1] "2023-10-11 11:34:53.189271 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:36:39.733939 ... ending   pcubature -- tol=0.09 -- ret.val is: 3675.93591"
# 1:     3675.9353: -1.07296  1.00420 0.968245  1.30025 0.308327
# [1] "2023-10-11 11:36:39.734164 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:38:26.339045 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.94682"
# [1] "2023-10-11 11:38:26.339276 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:40:13.055671 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.94684"
# [1] "2023-10-11 11:40:13.055882 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:41:59.723995 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.94685"
# [1] "2023-10-11 11:41:59.724213 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:43:46.18868 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.94689"
# [1] "2023-10-11 11:43:46.188895 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:45:32.585884 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.94694"
# [1] "2023-10-11 11:45:32.586107 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:47:19.368291 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.94704"
# 2:     3673.9468: -1.06407 0.903002 0.932116  1.32084 0.380330
# [1] "2023-10-11 11:47:19.368532 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:49:06.008158 ... ending   pcubature -- tol=0.09 -- ret.val is: 3674.0785"
# [1] "2023-10-11 11:49:06.008394 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:50:52.727401 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.46444"
# [1] "2023-10-11 11:50:52.727631 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:52:39.447943 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.46438"
# [1] "2023-10-11 11:52:39.448164 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:54:26.282312 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.46448"
# [1] "2023-10-11 11:54:26.282541 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:56:13.098666 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.46449"
# [1] "2023-10-11 11:56:13.098888 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:57:59.907894 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.46466"
# [1] "2023-10-11 11:57:59.908113 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 11:59:46.568775 ... ending   pcubature -- tol=0.09 -- ret.val is: 3673.4647"
# 3:     3673.4644: -1.02243 0.889744 0.910138  1.33277 0.415877
# [1] "2023-10-11 11:59:46.569013 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:01:33.626545 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.89978"
# [1] "2023-10-11 12:01:33.626777 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:03:20.704967 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.89975"
# [1] "2023-10-11 12:03:20.705183 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:05:07.638127 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.89978"
# [1] "2023-10-11 12:05:07.638347 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:06:54.692795 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.89985"
# [1] "2023-10-11 12:06:54.693015 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:08:41.499194 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.89991"
# [1] "2023-10-11 12:08:41.49941 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:10:28.325247 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.89995"
# 4:     3672.8998: -1.05423 0.876339 0.900124  1.35249 0.461906
# [1] "2023-10-11 12:10:28.325492 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:12:15.264472 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.49933"
# [1] "2023-10-11 12:12:15.264692 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:14:02.190624 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.4993"
# [1] "2023-10-11 12:14:02.19083 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:15:49.113758 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.49936"
# [1] "2023-10-11 12:15:49.113969 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:17:36.156631 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.49935"
# [1] "2023-10-11 12:17:36.156839 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:19:22.935065 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.4995"
# [1] "2023-10-11 12:19:22.935282 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:21:09.984468 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.49942"
# 5:     3672.4993: -1.02260 0.878296 0.871742  1.37424 0.500836
# [1] "2023-10-11 12:21:09.984694 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-11 12:22:56.510018 ... ending   pcubature -- tol=0.09 -- ret.val is: 3672.49933"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -1.0225960  0.8782959  0.8717420  1.3742424  0.5008356
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 3672.499
#
# $fevals
# function
# 8
#
# $gevals
# gradient
# 30
#
# $nitns
# [1] 5
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 4117.675
#
# Assemble the answers
# Intercept       b_p     var1     var2    corr12    value fevals gevals
# nlminb -1.022596 0.8782959 0.871742 1.374242 0.5008356 3672.499      8     30
# niter convcode kkt1 kkt2    xtime
# nlminb     5        0   NA   NA 4117.675
# >
#   > glm(cbind(r, n_r) ~ x1, summed_binom_dat, family=binomial(link="cauchit"))
#
# Call:  glm(formula = cbind(r, n_r) ~ x1, family = binomial(link = "cauchit"),
#            data = summed_binom_dat)
#
# Coefficients:
#   (Intercept)           x1
# -1.0873       0.9787
#
# Degrees of Freedom: 1599 Total (i.e. Null);  1598 Residual
# Null Deviance:      16340
# Residual Deviance: 16250        AIC: 20530


## then *Locally* run this part after inputting the values
## some things are commented out but it should make the graph
## at the end
phi_x_deluxe <- function(x, var1=0.871, var2=1.37, corr12=0.50, gamma_d = 1, gamma_h=1, a=1.00 ){

  (gamma_d)  / (gamma_h^a + (var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x)^(a/2) )^(1/a)

}

datavalueset <- c(0.1, 0.12, 0.2, 0.22, 0.3, 0.32, 0.4, 0.42)
xvalueset <- c(0, datavalueset, 1, 10)
xdomain <- seq(-2,4,0.05)

plot(NA, NA,col="grey", ylim=c(0,1),xlim=range(xdomain))
abline(h=0.50, lty=3, col="lightblue")
abline(v=min(datavalueset), lty=2, col="lightgray")
abline(v=max(datavalueset), lty=2, col="lightgray")

lines(xdomain, pcauchy(-2.063432+1.032551*xdomain), ylim=c(0,1), lty=3)
text(3.6,0.68, "(b0c+b1c*x)       ", col="black")

lines(xdomain, pcauchy( -1.022596+0.8782959*xdomain), ylim=c(0,1),col="orange",lty=3)
text(3.6,0.95, "(b0m+b1m*x)       ", col="orange")

lines(xdomain, phi_x_deluxe(xdomain), ylim=c(0,1), lty=2, col='cyan')
text(3.6,0.19, "            phi(x)", col="cyan")

lines(xdomain, pcauchy((-2.063432+1.032551*xdomain)*phi_x_deluxe(xdomain)), col="red",lty=3)
text(3.6,0.52, "(b0c+b1c*x)*phi(x)", col="red")


lines(xdomain, pcauchy(-1.0873+0.9787*xdomain), ylim=c(0,1), col="darkgreen",lty=3)
text(3.6,0.78, "glm(b'0m+b'1m*x)       ", col="darkgreen")


## phi ranges over the dataset:
phi_x_deluxe(min(datavalueset))
phi_x_deluxe(max(datavalueset))
abline(h=phi_x_deluxe(min(datavalueset)), lty=2, col="lightpink")
abline(h=phi_x_deluxe(max(datavalueset)), lty=2, col="lightpink")




########################################$$$$#########
## 2023-10-11A                                     ##
## _END_ two random parameters: run on cluster    ###
## bivariate cauchy marg/cond                      ##
#####################################################

########################################$$$$#########
## 2023-10-10A                                     ##
## START two random parameters: run on cluster    ###
## bivariate normal marg/cond                      ##
#####################################################


#
# > library(gnlrim)
# > library(mvpd)
# > library(libstable4u)
# > library(data.table)
# data.table 1.14.8 using 1 threads (see ?getDTthreads).  Latest news: r-datatable.com
# > # Derived expression for gamma
#   > g <- function(a) {
#     +   iu <- complex(real=0, imaginary=1)
#     +   return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
#     + }
# >
#   > bridgecloglog_rstable <- function(n, alpha, delta){
#     +   mult <- (delta/alpha)^(1/alpha)
#     +   X <- stabledist::rstable(n , alpha, beta=1, gamma=g(alpha), delta=0, pm=1)
#     +   Z <- log(mult * X)
#     +   Z
#     + }
# > ## add in beta random effect
#   > sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){
#     +
#       +   if(mrim=="Logistic-BIVARIATE_NORM"){
#         +     print("HERE")
#         +     G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
#         +     H <- function(x) plogis(x)
#         +   }
#     +   if(mrim=="Probit-BIVARIATE_NORM"){
#       +     print("Probit for the win")
#       +     G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
#       +     H <- function(x) pnorm(x)
#       +   }
#     +   if(mrim=="Stabit-BIVARIATE_STABLE"){
#       +     print("STABLE for the win")
#       +     print(paste0("alpha set to ", alpha))
#       +     G <- function(n, v1, v2, rho){mvpd::rmvss(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
#       +     H <- function(x) stable_cdf(x, c(alpha,0,1,0))
#       +   }
#     +
#       +   n <- n1 + n2
#       +   u <- round(apply(G(n,v1,v2,rho),2, function(W) rep(W,each=J)),2)
#       +
#         +   ## x <- c(rep(1, n1*J), rep(0, n2*J))
#         +   ##  x <-c(runif(n1*J, 0.5,1.5), runif(n2*J, 0,0.60))
#         +   x <-c(sample(c(1,2,3,4)/10,n1*J,TRUE), sample(c(1.2,2.2,3.2,4.2)/10,n2*J,TRUE))
#         +   eta <- round(a0 + a1*x,2)
#         +
#           +   eta_i <- round( (a0 + u[,1]) + (a1+u[,2])*x, 2)
#           +   py1 <- round(H(eta_i),2)
#           +   y <- rbinom(length(eta_i), 1, prob=py1 )
#           +
#             +   data.frame(id=rep(1:n, each=J),
#                            +              j = rep(1:J),
#                            +              x1 = x,
#                            +              eta = eta,
#                            +              u_i = u,
#                            +              eta_i = eta_i,
#                            +              py1 = py1,
#                            +              y=y
#                            +   )
#           +
#             + }
# >
#   >
#   > detach(summed_binom_dat)
# Error in detach(summed_binom_dat) : invalid 'name' argument
# > set.seed(1709)
# > binom_dat <-
#   +   #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
#   +   #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
#   +   sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1)
# [1] "STABLE for the win"
# [1] "alpha set to 1"
# > data.table::setDT(binom_dat)
# >
#   > summed_binom_dat <-
#   +   binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
# > data.table::setkey(summed_binom_dat,id, x1)
# > summed_binom_dat
# id   x1  r n_r
# 1:   1 0.10  6  15
# 2:   1 0.20  5  24
# 3:   1 0.30  5  16
# 4:   1 0.40  3  26
# 5:   2 0.10  8  19
# ---
#   1596: 399 0.42 11  19
# 1597: 400 0.12  7  23
# 1598: 400 0.22  6  22
# 1599: 400 0.32  4  17
# 1600: 400 0.42  2  19
# > attach(summed_binom_dat)
# >
#   > ybind <- cbind(r,n_r)
# > period_numeric <- x1
# >
#   > ## glmer with correlation between random intercept and random slope
#   > lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
# Generalized linear mixed model fit by maximum likelihood (Adaptive
#                                                           Gauss-Hermite Quadrature, nAGQ = 0) [glmerMod]
# Family: binomial  ( probit )
# Formula: cbind(r, n_r) ~ x1 + (x1 | id)
# Data: summed_binom_dat
# AIC       BIC    logLik  deviance  df.resid
# 7583.362  7610.251 -3786.681  7573.362      1595
# Random effects:
#   Groups Name        Std.Dev. Corr
# id     (Intercept) 0.8478
# x1          1.3740   0.18
# Number of obs: 1600, groups:  id, 400
# Fixed Effects:
#   (Intercept)           x1
# -0.8324       0.5497
# > lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)
# Generalized linear mixed model fit by maximum likelihood (Laplace
#                                                           Approximation) [glmerMod]
# Family: binomial  ( probit )
# Formula: cbind(r, n_r) ~ x1 + (x1 | id)
# Data: summed_binom_dat
# AIC       BIC    logLik  deviance  df.resid
# 7583.292  7610.181 -3786.646  7573.292      1595
# Random effects:
#   Groups Name        Std.Dev. Corr
# id     (Intercept) 0.8463
# x1          1.3683   0.18
# Number of obs: 1600, groups:  id, 400
# Fixed Effects:
#   (Intercept)           x1
# -0.8440       0.5522
# >
#   >
#   > ## starting at 12:52pm 2023-10-10
#   > ## starting at 12:52pm 2023-10-10
#   > ## starting at 12:52pm 2023-10-10
#   >
#   >
#   > (rand.int.rand.slopes.nonzero.corr.CUBA.MARGINAL <-
#        +     gnlrim::gnlrem(y=ybind,
#                             +                    mu = ~ pnorm(
#                               +                      (Intercept + period_numeric*b_p)*
#                                 +                        sqrt(1 + var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric ) +
#                                 +                        rand1 + period_numeric*rand2
#                               +                    ),
#                             +                    pmu = c(Intercept=-0.95, b_p=0.55, var1=1, var2=1, corr12= 0.20),
#                             +                    pmix=c(var1=1, var2=1, corr12= 0.20),
#                             +                    p_uppb = c(  0,   2, 4.00, 4.00, 0.90),
#                             +                    p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                             +                    distribution="binomial",
#                             +                    nest=id,
#                             +                    random=c("rand1", "rand2"),
#                             +                    mixture="bivariate-normal-corr",
#                             +                    ooo=TRUE,
#                             +                    compute_hessian = FALSE,
#                             +                    compute_kkt = FALSE,
#                             +                    trace=1,
#                             +                    method='nlminb',
#                             +                    int2dmethod="cuba"
#                             +     )
#      + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -4   <? -0.95   <? 0     In Bounds
# par[ 2 ]:  -2   <? 0.55   <? 2     In Bounds
# par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? 0.2   <? 0.9     In Bounds
# [1] "2023-10-10 12:52:07.246674 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:52:32.692152 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.69897   log bounds ratio= 0.3467875
# Method:  nlminb
# [1] "2023-10-10 12:52:32.807845 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:52:58.244022 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# [1] "2023-10-10 12:52:58.244244 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:53:23.578015 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# [1] "2023-10-10 12:53:23.578226 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:53:49.220704 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# [1] "2023-10-10 12:53:49.220912 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:54:14.938711 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# [1] "2023-10-10 12:54:14.938924 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:54:40.390297 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# [1] "2023-10-10 12:54:40.390511 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:55:06.003416 ... ending   pcubature -- tol=0 -- ret.val is: 3837.79604"
# 0:     3837.7960: -0.950000 0.550000  1.00000  1.00000 0.200000
# [1] "2023-10-10 12:55:06.003663 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:55:18.68468 ... ending   pcubature -- tol=0 -- ret.val is: 3803.95171"
# [1] "2023-10-10 12:55:18.68491 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:55:31.352621 ... ending   pcubature -- tol=0 -- ret.val is: 3803.94405"
# [1] "2023-10-10 12:55:31.352834 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:55:44.082723 ... ending   pcubature -- tol=0 -- ret.val is: 3803.94856"
# [1] "2023-10-10 12:55:44.082931 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:55:56.676746 ... ending   pcubature -- tol=0 -- ret.val is: 3803.95219"
# [1] "2023-10-10 12:55:56.676961 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:56:09.383235 ... ending   pcubature -- tol=0 -- ret.val is: 3803.95294"
# [1] "2023-10-10 12:56:09.383454 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:56:22.104123 ... ending   pcubature -- tol=0 -- ret.val is: 3803.95095"
# 1:     3803.9517: -0.466261 0.609965 0.896501  1.03234 0.174653
# [1] "2023-10-10 12:56:22.104365 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:56:47.398168 ... ending   pcubature -- tol=0 -- ret.val is: 3795.79822"
# [1] "2023-10-10 12:56:47.398388 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:57:12.67701 ... ending   pcubature -- tol=0 -- ret.val is: 3795.79823"
# [1] "2023-10-10 12:57:12.677227 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:57:38.02516 ... ending   pcubature -- tol=0 -- ret.val is: 3795.79761"
# [1] "2023-10-10 12:57:38.02538 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:58:03.254251 ... ending   pcubature -- tol=0 -- ret.val is: 3795.79756"
# [1] "2023-10-10 12:58:03.254462 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 12:58:28.607691 ... ending   pcubature -- tol=0 -- ret.val is: 3795.79885"
# [1] "2023-10-10 12:58:28.607903 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 13:48:46.064348 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 13:48:58.837123 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11016"
# [1] "2023-10-10 13:48:58.837332 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 13:49:11.603476 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11016"
# [1] "2023-10-10 13:49:11.603691 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 13:49:24.27714 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11016"
# [1] "2023-10-10 13:49:24.277344 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 13:49:37.079996 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11016"
# 28:     3782.1102: -0.645891 0.568588 0.718746  1.92642 0.169825
# [1] "2023-10-10 13:49:37.080228 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-10-10 13:49:49.864315 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11016"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -0.6458908  0.5685879  0.7187462  1.9264188  0.1698254
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 3782.11
#
# $fevals
# function
# 37
#
# $gevals
# gradient
# 182
#
# $nitns
# [1] 28
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 3420.541
#
# Assemble the answers
# Intercept       b_p      var1     var2    corr12   value fevals gevals
# nlminb -0.6458908 0.5685879 0.7187462 1.926419 0.1698254 3782.11     37    182
# niter convcode kkt1 kkt2    xtime
# nlminb    28        0   NA   NA 3420.541
# > 0.8463^2
# [1] 0.7162237
# > 1.3683^2
# [1] 1.872245
# > 1.3740^2
# [1] 1.887876
# >
#   >
#   > ## now try conditional
#   > ## now try conditional
#   > ## now try conditional
#   >
#   >
#   > (rand.int.rand.slopes.nonzero.corr.CUBA.CONDITIONAL <-
#        +   gnlrim::gnlrem(y=ybind,
#                           +                  mu = ~ pnorm(
#                             +                    (Intercept + period_numeric*b_p) +
#                               +
#                               +                      rand1 + period_numeric*rand2
#                             +                  ),
#                           +                  pmu = c(Intercept=-0.95, b_p=0.55, var1=1, var2=1, corr12= 0.20),
#                           +                  pmix=c(var1=1, var2=1, corr12= 0.20),
#                           +                  p_uppb = c(  0,   2, 4.00, 4.00, 0.90),
#                           +                  p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                           +                  distribution="binomial",
#                           +                  nest=id,
#                           +                  random=c("rand1", "rand2"),
#                           +                  mixture="bivariate-normal-corr",
#                           +                  ooo=TRUE,
#                           +                  compute_hessian = FALSE,
#                           +                  compute_kkt = FALSE,
#                           +                  trace=1,
#                           +                  method='nlminb',
#                           +                  int2dmethod="cuba"
#                           +   )
#      + )
#
# Parameters are Intercept b_p
# Error in gnlrim::gnlrem(y = ybind, mu = ~pnorm((Intercept + period_numeric *  :
#                                                   pmu should have 2 estimates
#                                                 > (rand.int.rand.slopes.nonzero.corr.CUBA.CONDITIONAL <-
#                                                      +   gnlrim::gnlrem(y=ybind,
#                                                                         +                  mu = ~ pnorm(
#                                                                           +                    (Intercept + period_numeric*b_p) +
#                                                                             +
#                                                                             +                      rand1 + period_numeric*rand2
#                                                                           +                  ),
#                                                                         +                  pmu = c(Intercept=-0.95, b_p=0.55),
#                                                                         +                  pmix=c(var1=1, var2=1, corr12= 0.20),
#                                                                         +                  p_uppb = c(  0,   2, 4.00, 4.00, 0.90),
#                                                                         +                  p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#                                                                         +                  distribution="binomial",
#                                                                         +                  nest=id,
#                                                                         +                  random=c("rand1", "rand2"),
#                                                                         +                  mixture="bivariate-normal-corr",
#                                                                         +                  ooo=TRUE,
#                                                                         +                  compute_hessian = FALSE,
#                                                                         +                  compute_kkt = FALSE,
#                                                                         +                  trace=1,
#                                                                         +                  method='nlminb',
#                                                                         +                  int2dmethod="cuba"
#                                                                         +   )
#                                                    + )
#                                                 fn is  fn1
#                                                 Looking for method =  nlminb
#                                                 Function has  5  arguments
#                                                 par[ 1 ]:  -4   <? -0.95   <? 0     In Bounds
#                                                 par[ 2 ]:  -2   <? 0.55   <? 2     In Bounds
#                                                 par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
#                                                 par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
#                                                 par[ 5 ]:  -0.9   <? 0.2   <? 0.9     In Bounds
#                                                 [1] "2023-10-10 14:03:27.258739 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:03:52.162627 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 Analytic gradient not made available.
#                                                 Analytic Hessian not made available.
#                                                 Scale check -- log parameter ratio= 0.69897   log bounds ratio= 0.3467875
#                                                 Method:  nlminb
#                                                 [1] "2023-10-10 14:03:52.276128 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:04:17.225 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 [1] "2023-10-10 14:04:17.225217 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:04:42.138645 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 [1] "2023-10-10 14:04:42.138861 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:05:07.057225 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 [1] "2023-10-10 14:05:07.057452 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:05:31.96138 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 [1] "2023-10-10 14:05:31.961603 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:05:56.938058 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 [1] "2023-10-10 14:05:56.938272 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:06:21.890938 ... ending   pcubature -- tol=0 -- ret.val is: 3797.25428"
#                                                 0:     3797.2543: -0.950000 0.550000  1.00000  1.00000 0.200000
#                                                 [1] "2023-10-10 14:06:21.891178 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:06:34.45409 ... ending   pcubature -- tol=0 -- ret.val is: 3802.27595"
#                                                 [1] "2023-10-10 14:06:34.45431 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:06:59.321851 ... ending   pcubature -- tol=0 -- ret.val is: 3790.03651"
#                                                 [1] "2023-10-10 14:06:59.322073 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:07:24.308886 ... ending   pcubature -- tol=0 -- ret.val is: 3790.03562"
#                                                 [1] "2023-10-10 14:07:24.309094 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:07:49.317374 ... ending   pcubature -- tol=0 -- ret.val is: 3790.03638"
#                                                 [1] "2023-10-10 14:07:49.317593 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:08:14.283301 ... ending   pcubature -- tol=0 -- ret.val is: 3790.03774"
#                                                 [1] "2023-10-10 14:08:14.283527 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:08:39.248438 ... ending   pcubature -- tol=0 -- ret.val is: 3790.0375"
#                                                 [1] "2023-10-10 14:08:39.24866 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:09:04.122351 ... ending   pcubature -- tol=0 -- ret.val is: 3790.03647"
#                                                 1:     3790.0365: -0.809497 0.581045 0.862906  1.06795 0.171303
#                                                 [1] "2023-10-10 14:09:04.122596 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:09:29.130384 ... ending   pcubature -- tol=0 -- ret.val is: 3787.82137"
#                                                 [1] "2023-10-10 14:09:29.130611 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:09:54.070104 ... ending   pcubature -- tol=0 -- ret.val is: 3787.82126"
#                                                 [1] "2023-10-10 14:09:54.070323 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:10:19.008267 ... ending   pcubature -- tol=0 -- ret.val is: 3787.8213"
#                                                 [1] "2023-10-10 14:10:19.008483 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:10:43.940641 ... ending   pcubature -- tol=0 -- ret.val is: 3787.82126"
#                                                 [1] "2023-10-10 14:10:43.940855 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:11:08.967751 ... ending   pcubature -- tol=0 -- ret.val is: 3787.82224"
#                                                 [1] "2023-10-10 14:11:08.967969 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:11:33.905201 ... ending   pcubature -- tol=0 -- ret.val is: 3787.82155"
#                                                 2:     3787.8214: -0.913303 0.565213 0.719413  1.18328 0.175914
#                                                 [1] "2023-10-10 14:11:33.905436 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:11:58.777467 ... ending   pcubature -- tol=0 -- ret.val is: 3785.49372"
#                                                 [1] "2023-10-10 14:11:58.777699 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:12:23.760679 ... ending   pcubature -- tol=0 -- ret.val is: 3785.49363"
#                                                 [1] "2023-10-10 14:12:23.760897 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:12:48.611654 ... ending   pcubature -- tol=0 -- ret.val is: 3785.49365"
#                                                 [1] "2023-10-10 14:12:48.611873 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:13:13.474517 ... ending   pcubature -- tol=0 -- ret.val is: 3785.49367"
#                                                 [1] "2023-10-10 14:13:13.474731 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:13:38.435634 ... ending   pcubature -- tol=0 -- ret.val is: 3785.49415"
#                                                 [1] "2023-10-10 14:13:38.435848 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:14:03.361964 ... ending   pcubature -- tol=0 -- ret.val is: 3785.49358"
#                                                 3:     3785.4937: -0.791894 0.573469 0.751732  1.33583 0.252163
#                                                 [1] "2023-10-10 14:14:03.362204 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:14:28.270411 ... ending   pcubature -- tol=0 -- ret.val is: 3783.77302"
#                                                 [1] "2023-10-10 14:14:28.270635 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:14:53.280552 ... ending   pcubature -- tol=0 -- ret.val is: 3783.77297"
#                                                 [1] "2023-10-10 14:14:53.28076 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:15:18.204954 ... ending   pcubature -- tol=0 -- ret.val is: 3783.77299"
#                                                 [1] "2023-10-10 14:15:18.20517 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:15:43.110496 ... ending   pcubature -- tol=0 -- ret.val is: 3783.77311"
#                                                 [1] "2023-10-10 14:15:43.110717 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:16:07.993413 ... ending   pcubature -- tol=0 -- ret.val is: 3783.77316"
#                                                 [1] "2023-10-10 14:16:07.993634 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:16:32.880492 ... ending   pcubature -- tol=0 -- ret.val is: 3783.7731"
#                                                 4:     3783.7730: -0.881611 0.559697 0.758799  1.52520 0.224057
#                                                 [1] "2023-10-10 14:16:32.880732 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:16:45.432345 ... ending   pcubature -- tol=0 -- ret.val is: 3783.86111"
#                                                 [1] "2023-10-10 14:16:45.432579 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:17:10.269422 ... ending   pcubature -- tol=0 -- ret.val is: 3783.20673"
#                                                 [1] "2023-10-10 14:17:10.269888 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:17:35.299181 ... ending   pcubature -- tol=0 -- ret.val is: 3783.20669"
#                                                 [1] "2023-10-10 14:17:35.299394 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:18:00.208107 ... ending   pcubature -- tol=0 -- ret.val is: 3783.2067"
#                                                 [1] "2023-10-10 14:18:00.20832 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:18:25.20188 ... ending   pcubature -- tol=0 -- ret.val is: 3783.20668"
#                                                 [1] "2023-10-10 14:18:25.202095 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:18:50.21269 ... ending   pcubature -- tol=0 -- ret.val is: 3783.20687"
#                                                 [1] "2023-10-10 14:18:50.212903 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:19:15.044116 ... ending   pcubature -- tol=0 -- ret.val is: 3783.2066"
#                                                 5:     3783.2067: -0.810901 0.558003 0.721573  1.56415 0.171811
#                                                 [1] "2023-10-10 14:19:15.044357 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:19:27.613025 ... ending   pcubature -- tol=0 -- ret.val is: 3782.66917"
#                                                 [1] "2023-10-10 14:19:27.613241 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:19:40.153576 ... ending   pcubature -- tol=0 -- ret.val is: 3782.66915"
#                                                 [1] "2023-10-10 14:19:40.153786 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:19:52.715585 ... ending   pcubature -- tol=0 -- ret.val is: 3782.66909"
#                                                 [1] "2023-10-10 14:19:52.715798 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:20:05.1999 ... ending   pcubature -- tol=0 -- ret.val is: 3782.66916"
#                                                 [1] "2023-10-10 14:20:05.200106 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:20:17.757064 ... ending   pcubature -- tol=0 -- ret.val is: 3782.66931"
#                                                 [1] "2023-10-10 14:20:17.757272 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:20:30.288244 ... ending   pcubature -- tol=0 -- ret.val is: 3782.66921"
#                                                 6:     3782.6692: -0.851256 0.537617 0.740275  1.65168 0.147714
#                                                 [1] "2023-10-10 14:20:30.288478 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:20:55.176038 ... ending   pcubature -- tol=0 -- ret.val is: 3782.47647"
#                                                 [1] "2023-10-10 14:20:55.176261 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:21:20.147663 ... ending   pcubature -- tol=0 -- ret.val is: 3782.47645"
#                                                 [1] "2023-10-10 14:21:20.147873 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:21:45.074424 ... ending   pcubature -- tol=0 -- ret.val is: 3782.47642"
#                                                 [1] "2023-10-10 14:21:45.074657 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:22:10.095118 ... ending   pcubature -- tol=0 -- ret.val is: 3782.4765"
#                                                 [1] "2023-10-10 14:22:10.095335 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:22:35.022801 ... ending   pcubature -- tol=0 -- ret.val is: 3782.47653"
#                                                 [1] "2023-10-10 14:22:35.023011 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:22:59.917949 ... ending   pcubature -- tol=0 -- ret.val is: 3782.47634"
#                                                 7:     3782.4765: -0.831999 0.579978 0.725026  1.71807 0.209603
#                                                 [1] "2023-10-10 14:22:59.918188 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:23:12.465099 ... ending   pcubature -- tol=0 -- ret.val is: 3782.25986"
#                                                 [1] "2023-10-10 14:23:12.465314 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:23:25.038697 ... ending   pcubature -- tol=0 -- ret.val is: 3782.25985"
#                                                 [1] "2023-10-10 14:23:25.038913 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:23:37.593635 ... ending   pcubature -- tol=0 -- ret.val is: 3782.25982"
#                                                 [1] "2023-10-10 14:23:37.593845 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:23:50.190609 ... ending   pcubature -- tol=0 -- ret.val is: 3782.2598"
#                                                 [1] "2023-10-10 14:23:50.190814 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:24:02.641546 ... ending   pcubature -- tol=0 -- ret.val is: 3782.25989"
#                                                 [1] "2023-10-10 14:24:02.641755 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:24:15.181537 ... ending   pcubature -- tol=0 -- ret.val is: 3782.25986"
#                                                 8:     3782.2599: -0.839498 0.519580 0.702422  1.79593 0.190681
#                                                 [1] "2023-10-10 14:24:15.181767 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:24:27.731253 ... ending   pcubature -- tol=0 -- ret.val is: 3782.17427"
#                                                 [1] "2023-10-10 14:24:27.731475 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:24:40.193587 ... ending   pcubature -- tol=0 -- ret.val is: 3782.17426"
#                                                 [1] "2023-10-10 14:24:40.193809 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:24:52.750163 ... ending   pcubature -- tol=0 -- ret.val is: 3782.17426"
#                                                 [1] "2023-10-10 14:24:52.750376 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:25:05.194252 ... ending   pcubature -- tol=0 -- ret.val is: 3782.17424"
#                                                 [1] "2023-10-10 14:25:05.194457 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:25:17.729067 ... ending   pcubature -- tol=0 -- ret.val is: 3782.17429"
#                                                 [1] "2023-10-10 14:25:17.729277 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:25:30.260138 ... ending   pcubature -- tol=0 -- ret.val is: 3782.17425"
#                                                 9:     3782.1743: -0.841044 0.569436 0.747083  1.85610 0.140360
#                                                 [1] "2023-10-10 14:25:30.26038 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:25:42.69796 ... ending   pcubature -- tol=0 -- ret.val is: 3782.74231"
#                                                 [1] "2023-10-10 14:25:42.698185 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:25:55.232464 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11886"
#                                                 [1] "2023-10-10 14:25:55.232689 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:26:07.650964 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11885"
#                                                 [1] "2023-10-10 14:26:07.651172 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:26:20.174198 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11887"
#                                                 [1] "2023-10-10 14:26:20.174407 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:26:32.68154 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11887"
#                                                 [1] "2023-10-10 14:26:32.681751 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:26:45.19891 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11888"
#                                                 [1] "2023-10-10 14:26:45.199123 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:26:57.60809 ... ending   pcubature -- tol=0 -- ret.val is: 3782.11889"
#                                                 10:     3782.1189: -0.844241 0.563709 0.734348  1.85912 0.145594
#                                                 [1] "2023-10-10 14:26:57.608322 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:27:10.128474 ... ending   pcubature -- tol=0 -- ret.val is: 3782.08646"
#                                                 [1] "2023-10-10 14:27:10.128696 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:27:22.65741 ... ending   pcubature -- tol=0 -- ret.val is: 3782.08646"
#                                                 [1] "2023-10-10 14:27:22.657624 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:27:35.18873 ... ending   pcubature -- tol=0 -- ret.val is: 3782.08647"
#                                                 [1] "2023-10-10 14:27:35.188935 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:27:47.618965 ... ending   pcubature -- tol=0 -- ret.val is: 3782.08647"
#                                                 [1] "2023-10-10 14:27:47.61919 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:28:00.136037 ... ending   pcubature -- tol=0 -- ret.val is: 3782.08648"
#                                                 [1] "2023-10-10 14:28:00.136243 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:28:12.665865 ... ending   pcubature -- tol=0 -- ret.val is: 3782.08648"
#                                                 11:     3782.0865: -0.840023 0.557702 0.726910  1.86430 0.155878
#                                                 [1] "2023-10-10 14:28:12.666101 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:28:25.097805 ... ending   pcubature -- tol=0 -- ret.val is: 3782.07945"
#                                                 [1] "2023-10-10 14:28:25.098029 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:28:37.613702 ... ending   pcubature -- tol=0 -- ret.val is: 3782.07944"
#                                                 [1] "2023-10-10 14:28:37.613924 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:28:50.150319 ... ending   pcubature -- tol=0 -- ret.val is: 3782.07945"
#                                                 [1] "2023-10-10 14:28:50.150542 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:29:02.569545 ... ending   pcubature -- tol=0 -- ret.val is: 3782.07945"
#                                                 [1] "2023-10-10 14:29:02.569758 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:29:15.100914 ... ending   pcubature -- tol=0 -- ret.val is: 3782.07946"
#                                                 [1] "2023-10-10 14:29:15.101122 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:29:27.630985 ... ending   pcubature -- tol=0 -- ret.val is: 3782.07945"
#                                                 12:     3782.0794: -0.847022 0.551544 0.722916  1.87175 0.165008
#                                                 [1] "2023-10-10 14:29:27.631229 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:29:40.052481 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06449"
#                                                 [1] "2023-10-10 14:29:40.052706 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:29:52.577616 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06449"
#                                                 [1] "2023-10-10 14:29:52.577826 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:30:05.125895 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06449"
#                                                 [1] "2023-10-10 14:30:05.126102 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:30:17.67471 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06449"
#                                                 [1] "2023-10-10 14:30:17.674922 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:30:30.222351 ... ending   pcubature -- tol=0 -- ret.val is: 3782.0645"
#                                                 [1] "2023-10-10 14:30:30.222564 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:30:42.659299 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06446"
#                                                 [1] "2023-10-10 14:30:42.659515 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:30:55.199368 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06452"
#                                                 13:     3782.0645: -0.841104 0.546893 0.722105  1.88477 0.161132
#                                                 [1] "2023-10-10 14:30:55.199603 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:31:07.751467 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06536"
#                                                 [1] "2023-10-10 14:31:07.751697 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:31:20.204088 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06118"
#                                                 [1] "2023-10-10 14:31:20.204306 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:31:32.760091 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06119"
#                                                 [1] "2023-10-10 14:31:32.760303 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:31:45.325174 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06118"
#                                                 [1] "2023-10-10 14:31:45.325387 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:31:57.778189 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06118"
#                                                 [1] "2023-10-10 14:31:57.778408 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:32:10.337106 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06119"
#                                                 [1] "2023-10-10 14:32:10.337536 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:32:23.018048 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06118"
#                                                 [1] "2023-10-10 14:32:23.018261 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:32:35.468837 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06119"
#                                                 [1] "2023-10-10 14:32:35.469057 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:32:47.983391 ... ending   pcubature -- tol=0 -- ret.val is: 3782.06118"
#                                                 14:     3782.0612: -0.839486 0.550043 0.723051  1.88962 0.165212
#                                                 [1] "2023-10-10 14:32:47.983638 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:33:00.530707 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05788"
#                                                 [1] "2023-10-10 14:33:00.530922 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:33:13.064824 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05788"
#                                                 [1] "2023-10-10 14:33:13.065032 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:33:25.61622 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05789"
#                                                 [1] "2023-10-10 14:33:25.616437 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:33:38.058577 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05788"
#                                                 [1] "2023-10-10 14:33:38.058796 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:33:50.615139 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05788"
#                                                 [1] "2023-10-10 14:33:50.615349 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:34:03.064604 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05789"
#                                                 [1] "2023-10-10 14:34:03.064815 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:34:15.611222 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05788"
#                                                 [1] "2023-10-10 14:34:15.611443 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:34:28.160015 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05788"
#                                                 15:     3782.0579: -0.841889 0.551464 0.719984  1.89560 0.164415
#                                                 [1] "2023-10-10 14:34:28.160245 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:34:40.595241 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:34:40.595459 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:34:53.139882 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:34:53.140096 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:35:05.686802 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:35:05.687012 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:35:18.22707 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:35:18.227278 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:35:30.777441 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:35:30.777659 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:35:43.221007 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:35:43.22122 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:35:55.761081 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:35:55.761299 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:36:08.2996 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05497"
#                                                 [1] "2023-10-10 14:36:08.299808 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:36:20.731111 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 [1] "2023-10-10 14:36:20.731322 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:36:33.285393 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05496"
#                                                 16:     3782.0550: -0.840563 0.549603 0.721710  1.90233 0.164377
#                                                 [1] "2023-10-10 14:36:33.285631 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:36:45.823858 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05298"
#                                                 [1] "2023-10-10 14:36:45.824078 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:36:58.255229 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05297"
#                                                 [1] "2023-10-10 14:36:58.25544 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:37:10.783699 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05298"
#                                                 [1] "2023-10-10 14:37:10.783908 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:37:23.222308 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05297"
#                                                 [1] "2023-10-10 14:37:23.222528 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:37:35.803899 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05298"
#                                                 [1] "2023-10-10 14:37:35.805807 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:37:48.300796 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05297"
#                                                 [1] "2023-10-10 14:37:48.301009 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:38:00.833806 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05298"
#                                                 [1] "2023-10-10 14:38:00.834022 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:38:13.361488 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05298"
#                                                 [1] "2023-10-10 14:38:13.361717 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:38:25.804032 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05297"
#                                                 [1] "2023-10-10 14:38:25.80425 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:38:38.327992 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05298"
#                                                 17:     3782.0530: -0.841563 0.550007 0.720042  1.90930 0.163365
#                                                 [1] "2023-10-10 14:38:38.328326 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:38:50.854653 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:38:50.854877 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:39:03.383638 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:39:03.383857 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:39:15.786388 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:39:15.786617 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:39:28.302895 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05174"
#                                                 [1] "2023-10-10 14:39:28.303105 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:39:40.826038 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05172"
#                                                 [1] "2023-10-10 14:39:40.826246 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:39:53.24529 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:39:53.245497 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:40:05.773686 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:40:05.773897 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:40:18.305681 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:40:18.305907 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:40:30.726542 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 [1] "2023-10-10 14:40:30.72676 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:40:43.248691 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05173"
#                                                 18:     3782.0517: -0.840197 0.552004 0.719670  1.91620 0.163436
#                                                 [1] "2023-10-10 14:40:43.248935 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:40:55.782885 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:40:55.783098 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:41:08.186895 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:41:08.187101 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:41:20.703584 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:41:20.703806 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:41:33.107345 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:41:33.107558 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:41:45.620861 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:41:45.621068 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:41:58.140467 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:41:58.140683 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:42:10.548711 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:42:10.548919 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:42:23.066214 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:42:23.066428 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:42:35.581015 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 [1] "2023-10-10 14:42:35.581228 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:42:47.987924 ... ending   pcubature -- tol=0 -- ret.val is: 3782.05033"
#                                                 19:     3782.0503: -0.841303 0.549970 0.721372  1.92293 0.163111
#                                                 [1] "2023-10-10 14:42:47.988168 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:43:00.504992 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:43:00.505221 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:43:12.919988 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:43:12.920192 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:43:25.447716 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:43:25.44792 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:43:37.978465 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:43:37.978684 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:43:50.511601 ... ending   pcubature -- tol=0 -- ret.val is: 3782.0498"
#                                                 [1] "2023-10-10 14:43:50.511807 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:44:02.934576 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:44:02.934785 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:44:15.495779 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:44:15.495989 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:44:28.017251 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:44:28.017472 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:44:40.693481 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:44:40.693715 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:44:53.624917 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 [1] "2023-10-10 14:44:53.625129 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:45:06.54079 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04979"
#                                                 20:     3782.0498: -0.840529 0.549327 0.721052  1.93013 0.162350
#                                                 [1] "2023-10-10 14:45:06.541032 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:45:19.347847 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:45:19.348068 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:45:32.185448 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:45:32.185679 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:45:44.986706 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:45:44.986924 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:45:57.977231 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:45:57.977454 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:46:10.914158 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:46:10.914381 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:46:23.598474 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:46:23.598674 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:46:36.486903 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:46:36.487113 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:46:49.271968 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:46:49.272187 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:47:02.020411 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:47:02.020625 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:47:14.874708 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:47:14.874915 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:47:27.568192 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 21:     3782.0497: -0.840853 0.549731 0.720741  1.93217 0.161993
#                                                 [1] "2023-10-10 14:47:27.568427 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:47:40.424602 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:47:40.424818 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:47:53.278435 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:47:53.278675 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:48:06.124281 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:48:06.124552 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:48:18.69372 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:48:18.693931 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:48:31.444313 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:48:31.444537 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:48:44.186872 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:48:44.187067 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:48:56.946695 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:48:56.946914 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:49:09.568438 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:49:09.568667 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:49:22.326414 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:49:22.326635 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:49:35.105477 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:49:35.105699 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:49:47.827239 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 22:     3782.0497: -0.840792 0.549620 0.720750  1.93197 0.162080
#                                                 [1] "2023-10-10 14:49:47.82746 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:50:00.599865 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:50:00.600094 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:50:13.35555 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:50:13.355754 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:50:25.90817 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:50:25.908382 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:50:38.543511 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:50:38.543723 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:50:51.251332 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:50:51.251559 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:51:03.85216 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:51:03.852368 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:51:16.570397 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:51:16.570627 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:51:29.253925 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:51:29.254138 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:51:41.943352 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:51:41.943565 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:51:54.534525 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 [1] "2023-10-10 14:51:54.534733 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:52:07.14849 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 23:     3782.0497: -0.840791 0.549632 0.720761  1.93197 0.162068
#                                                 [1] "2023-10-10 14:52:07.148726 ... starting pcubature for bivariate normal or cauchy"
#                                                 [1] "2023-10-10 14:52:19.682764 ... ending   pcubature -- tol=0 -- ret.val is: 3782.04974"
#                                                 Post processing for method  nlminb
#                                                 Successful convergence!
#                                                   Save results from method  nlminb
#                                                 $par
#                                                 Intercept        b_p       var1       var2     corr12
#                                                 -0.8407910  0.5496316  0.7207606  1.9319714  0.1620683
#
#                                                 $message
#                                                 [1] "relative convergence (4)"
#
#                                                 $convcode
#                                                 [1] 0
#
#                                                 $value
#                                                 [1] 3782.05
#
#                                                 $fevals
#                                                 function
#                                                 28
#
#                                                 $gevals
#                                                 gradient
#                                                 161
#
#                                                 $nitns
#                                                 [1] 23
#
#                                                 $kkt1
#                                                 [1] NA
#
#                                                 $kkt2
#                                                 [1] NA
#
#                                                 $xtimes
#                                                 user.self
#                                                 2892.59
#
#                                                 Assemble the answers
#                                                 Intercept       b_p      var1     var2    corr12   value fevals gevals
#                                                 nlminb -0.840791 0.5496316 0.7207606 1.931971 0.1620683 3782.05     28    161
#                                                 niter convcode kkt1 kkt2   xtime
#                                                 nlminb    23        0   NA   NA 2892.59
#
# > glm.fit <-
#   +   glm(cbind(r, n_r) ~ x1, summed_binom_dat, family=binomial(link="probit"))
# > summary(glm.fit)
#
# Call:
#   glm(formula = cbind(r, n_r) ~ x1, family = binomial(link = "probit"),
#       data = summed_binom_dat)
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.72507    0.01695 -42.781   <2e-16 ***
#   x1           0.54384    0.05932   9.168   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 16336  on 1599  degrees of freedom
# Residual deviance: 16252  on 1598  degrees of freedom
# AIC: 20527
#
# Number of Fisher Scoring iterations: 4

## then *Locally* run this part after inputting the values
## some things are commented out but it should make the graph
## at the end
phi_x <- function(x, var1=0.72, var2=1.93, corr12=0.162){

  1/sqrt(1 + var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x )

}

##xvalueset <- c(0, sort(unique(period_numeric)), 1, 10)
datavalueset <- c(0.1, 0.12, 0.2, 0.22, 0.3, 0.32, 0.4, 0.42)
xvalueset <- c(0, datavalueset, 1, 10)

cbind(
  xvalue = xvalueset,

  "phi(xvalue)" = phi_x(x=xvalueset),

  beta0_c = -0.8407910,
  beta1_c =   0.5496316,

  beta0_m =  -0.8407910*  phi_x(x=xvalueset),
  beta1_m =    0.5496316 *  phi_x(x=xvalueset)
)

xdomain <- seq(-2,4,0.05)
plot(xdomain, phi_x(xdomain), ylim=c(0,1))

library(data.table)
#dat.interim <- data.table::copy(summed_binom_dat[, r/(r+n_r), by=c("id","period_numeric")])

## do it either way:
#dat.interim2 <- data.table::copy(dat.interim[order(period_numeric),mean(V1), by=period_numeric])
#summed_binom_dat[, mean(r/(r+n_r)), by=c("period_numeric")]

dat.interim2[,marginal.fit:=pnorm( -0.6458908+0.5685879*period_numeric)]
dat.interim2[,conditional.fit:=pnorm( -0.8407910+0.5496316*period_numeric)]
setnames(dat.interim2, "V1", "empirical.avg")
setcolorder(dat.interim2, c("period_numeric", "conditional.fit","empirical.avg","marginal.fit" ))

dat.interim2[, glm.fit:=pnorm(-1.42435 + 1.10648 *period_numeric)]

print(dat.interim2, digits=2)
## plot(x1, r/(r+n_r),col="grey", ylim=c(0,1),xlim=range(xdomain))
plot(NA, NA,col="grey", ylim=c(0,1),xlim=range(xdomain))
abline(h=0.50, lty=3, col="lightblue")
abline(v=min(datavalueset), lty=2, col="lightgray")
abline(v=max(datavalueset), lty=2, col="lightgray")

lines(xdomain, pnorm( -0.8407910+0.5496316*xdomain), ylim=c(0,1))
text(3.6,0.68, "(b0c+b1c*x)       ", col="black")

lines(xdomain, pnorm( -0.6458908+0.5685879*xdomain), ylim=c(0,1),col="orange")
text(3.6,0.95, "(b0m+b1m*x)       ", col="orange")

lines(xdomain, phi_x(xdomain), ylim=c(0,1), lty=2, col='cyan')
text(3.6,0.19, "            phi(x)", col="cyan")

lines(xdomain, pnorm((-0.8407910+0.5496316*xdomain)*phi_x(xdomain)), col="red")
text(3.6,0.52, "(b0c+b1c*x)*phi(x)", col="red")


lines(xdomain, pnorm( -0.72507+0.54384*xdomain), ylim=c(0,1), col="darkgreen")
text(3.6,0.78, "glm(b'0m+b'1m*x)       ", col="darkgreen")


## phi ranges over the dataset:
phi_x(min(datavalueset))
phi_x(max(datavalueset))
abline(h=phi_x(min(datavalueset)), lty=2, col="lightpink")
abline(h=phi_x(max(datavalueset)), lty=2, col="lightpink")



########################################$$$$#########
## 2023-10-10A                                     ##
## _END_ two random parameters: run on cluster    ###
## bivariate normal marg/cond                      ##
#####################################################


########################################$$$$#########
## 2023-10-03A                                     ##
## START two random param: QUICKER                ###
#####################################################

## I'm starting to think that one can't estimate df
## but could still set it and test on a grid (?)

## tl;dr: false convergence?!? let's try not baby-ing the tolerances in 2023-10-02A
## and we can let df range

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t_df_gt_2"){
    print("t-BIVARIATE-t where shape is in terms for variance")
    a <- alpha
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c((a-2)/a*v1,(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
##detach(summed_binom_dat)
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  #sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t", alpha=8, v1=1, v2=2, rho=0.5)
  sim_mrim_data(20,20, J=10, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t", alpha=8, v1=1, v2=2, rho=0.5)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=nu_set),

                   pmu = c(Intercept=-1.2, b_p=1, nu_set=8),
                   pmix=c(nu_set=8, v1=1, v2=1, corr= 0.30),

                   p_uppb = c(  0,   2,    8, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2,    8, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba"
    )
)
########################################$$$$#########
## 2023-10-03A                                     ##
## END two random param: QUICKER                  ###
#####################################################

########################################$$$$#########
## 2023-10-02A                                     ##
## START two random param: fix-df bivariate-t sim ###
#####################################################

## I'm starting to think that one can't estimate df
## but could still set it and test on a grid (?)

## tl;dr: false convergence?!? let's try not baby-ing the tolerances in 2023-10-02A
## and we can let df range

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t_df_gt_2"){
    print("t-BIVARIATE-t where shape is in terms for variance")
    a <- alpha
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c((a-2)/a*v1,(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
##detach(summed_binom_dat)
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t", alpha=8, v1=1, v2=2, rho=0.5)

data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=nu_set),

                   pmu = c(Intercept=-1.2, b_p=1, nu_set=8),
                   pmix=c(nu_set=8, v1=1, v2=1, corr= 0.30),

                   p_uppb = c(  0,   2,   20, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2,    4, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba"
    )
)
# > (rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=nu_set),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, nu_set=8),
#                           +                    pmix=c(nu_set=8, v1=1, v2=1, corr= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2,   20, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2,    4, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-t-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba"
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  4   <? 8   <? 20     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-10-02 18:46:34 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 18:50:37 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 1.425969   log bounds ratio= 0.9488475
# Method:  nlminb
# [1] "2023-10-02 18:50:37 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 18:54:39 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# [1] "2023-10-02 18:54:39 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 18:58:38 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# [1] "2023-10-02 18:58:38 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:02:35 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# [1] "2023-10-02 19:02:35 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:06:31 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# [1] "2023-10-02 19:06:31 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:10:27 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# [1] "2023-10-02 19:10:27 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:14:24 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# [1] "2023-10-02 19:14:24 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:18:21 ... ending   pcubature -- tol=0 -- ret.val is: 3087.95878"
# 0:     3087.9588: -1.20000  1.00000  8.00000  1.00000  1.00000 0.300000
# [1] "2023-10-02 19:18:21 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:26:11 ... ending   pcubature -- tol=0 -- ret.val is: 3019.32991"
# [1] "2023-10-02 19:26:11 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:34:01 ... ending   pcubature -- tol=0 -- ret.val is: 3019.3336"
# [1] "2023-10-02 19:34:01 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:41:52 ... ending   pcubature -- tol=0 -- ret.val is: 3019.3309"
# [1] "2023-10-02 19:41:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:49:42 ... ending   pcubature -- tol=0 -- ret.val is: 3019.33005"
# [1] "2023-10-02 19:49:42 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 19:57:30 ... ending   pcubature -- tol=0 -- ret.val is: 3019.32985"
# [1] "2023-10-02 19:57:30 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 20:37:40 ... ending   pcubature -- tol=0 -- ret.val is: 3019.33145"
# [1] "2023-10-02 20:37:40 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:04:20 ... ending   pcubature -- tol=0 -- ret.val is: 3019.33093"
# 1:     3019.3299: -1.67074  1.02083  8.01219  1.10456  1.09636 0.387224
# [1] "2023-10-02 21:04:20 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:12:10 ... ending   pcubature -- tol=0 -- ret.val is: 3000.73981"
# [1] "2023-10-02 21:12:10 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:20:18 ... ending   pcubature -- tol=0 -- ret.val is: 3000.73974"
# [1] "2023-10-02 21:20:18 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:28:30 ... ending   pcubature -- tol=0 -- ret.val is: 3000.74002"
# [1] "2023-10-02 21:28:30 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:36:21 ... ending   pcubature -- tol=0 -- ret.val is: 3000.73981"
# [1] "2023-10-02 21:36:21 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:44:13 ... ending   pcubature -- tol=0 -- ret.val is: 3000.73976"
# [1] "2023-10-02 21:44:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:52:05 ... ending   pcubature -- tol=0 -- ret.val is: 3000.74034"
# [1] "2023-10-02 21:52:05 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 21:59:56 ... ending   pcubature -- tol=0 -- ret.val is: 3000.73981"
# 2:     3000.7398: -2.06890  1.17818  8.03084  1.04505  1.31049 0.517466
# [1] "2023-10-02 21:59:56 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:07:47 ... ending   pcubature -- tol=0 -- ret.val is: 3000.57675"
# [1] "2023-10-02 22:07:47 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:15:34 ... ending   pcubature -- tol=0 -- ret.val is: 3000.5766"
# [1] "2023-10-02 22:15:34 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:23:22 ... ending   pcubature -- tol=0 -- ret.val is: 3000.57653"
# [1] "2023-10-02 22:23:22 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:31:09 ... ending   pcubature -- tol=0 -- ret.val is: 3000.5767"
# [1] "2023-10-02 22:31:09 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:38:57 ... ending   pcubature -- tol=0 -- ret.val is: 3000.57646"
# [1] "2023-10-02 22:38:57 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:46:46 ... ending   pcubature -- tol=0 -- ret.val is: 3000.57696"
# [1] "2023-10-02 22:46:46 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 22:54:33 ... ending   pcubature -- tol=0 -- ret.val is: 3000.57686"
# [1] "2023-10-02 22:54:33 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 23:02:22 ... ending   pcubature -- tol=0 -- ret.val is: 3000.57665"
# 3:     3000.5768: -1.78921  1.33393  8.03358  1.18573  1.66664 0.547087
# [1] "2023-10-02 23:02:22 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 23:17:42 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21755"
# [1] "2023-10-02 23:17:42 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-02 23:33:03 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21752"
# [1] "2023-10-02 23:33:03 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 06:51:07 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21761"
# [1] "2023-10-03 06:51:07 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 07:06:51 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21755"
# [1] "2023-10-03 07:06:51 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 07:22:15 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21759"
# [1] "2023-10-03 07:22:15 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 07:37:36 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21776"
# [1] "2023-10-03 07:37:36 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 07:52:56 ... ending   pcubature -- tol=0 -- ret.val is: 2995.21769"
# 4:     2995.2176: -2.01803  1.28400  8.03998  1.13921  1.71624 0.492481
# [1] "2023-10-03 07:52:56 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 08:00:38 ... ending   pcubature -- tol=0 -- ret.val is: 2995.42258"
# [1] "2023-10-03 08:00:38 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 09:46:13 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44042"
# [1] "2023-10-03 09:46:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 09:53:58 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44039"
# [1] "2023-10-03 09:53:58 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:01:40 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44047"
# [1] "2023-10-03 10:01:40 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:09:23 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44039"
# [1] "2023-10-03 10:09:23 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:17:45 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44044"
# [1] "2023-10-03 10:17:45 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:25:52 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44047"
# [1] "2023-10-03 10:25:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:34:09 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44068"
# [1] "2023-10-03 10:34:09 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:42:16 ... ending   pcubature -- tol=0 -- ret.val is: 2994.44034"
# 5:     2994.4404: -1.94494  1.24864  8.03981  1.12016  1.77956 0.438425
# [1] "2023-10-03 10:42:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 10:58:49 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92394"
# [1] "2023-10-03 10:58:49 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 11:14:55 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92389"
# [1] "2023-10-03 11:14:55 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 11:30:52 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92398"
# [1] "2023-10-03 11:30:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 11:46:49 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92394"
# [1] "2023-10-03 11:46:49 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 12:02:45 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92392"
# [1] "2023-10-03 12:02:45 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 12:19:19 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92415"
# [1] "2023-10-03 12:19:19 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 12:36:19 ... ending   pcubature -- tol=0 -- ret.val is: 2993.92393"
# 6:     2993.9239: -2.01876  1.22606  8.04691  1.09538  1.86088 0.463924
# [1] "2023-10-03 12:36:19 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 12:52:00 ... ending   pcubature -- tol=0 -- ret.val is: 2993.14489"
# [1] "2023-10-03 12:52:00 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 13:07:40 ... ending   pcubature -- tol=0 -- ret.val is: 2993.14487"
# [1] "2023-10-03 13:07:40 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 13:23:15 ... ending   pcubature -- tol=0 -- ret.val is: 2993.14494"
# [1] "2023-10-03 13:23:15 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 13:38:46 ... ending   pcubature -- tol=0 -- ret.val is: 2993.14492"
# [1] "2023-10-03 13:38:46 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 13:54:18 ... ending   pcubature -- tol=0 -- ret.val is: 2993.14484"
# [1] "2023-10-03 13:54:18 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 14:09:52 ... ending   pcubature -- tol=0 -- ret.val is: 2993.14507"
# [1] "2023-10-03 14:09:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 14:25:27 ... ending   pcubature -- tol=0 -- ret.val is: 2993.1449"
# 7:     2993.1449: -1.93671  1.18606  8.05122  1.09239  1.93474 0.472438
# [1] "2023-10-03 14:25:27 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 14:33:20 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93198"
# [1] "2023-10-03 14:33:20 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 14:41:11 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93195"
# [1] "2023-10-03 14:41:11 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 14:49:03 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93194"
# [1] "2023-10-03 14:49:03 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 14:56:53 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93199"
# [1] "2023-10-03 14:56:53 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 15:04:44 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93197"
# [1] "2023-10-03 15:04:44 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 15:12:35 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93214"
# [1] "2023-10-03 15:12:35 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 15:20:25 ... ending   pcubature -- tol=0 -- ret.val is: 2991.93203"
# 8:     2991.9320: -1.97204  1.08462  8.07997  1.03511  2.13344 0.492909
# [1] "2023-10-03 15:20:25 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 15:36:01 ... ending   pcubature -- tol=0 -- ret.val is: 2992.81366"
# [1] "2023-10-03 15:36:01 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 15:43:50 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41309"
# [1] "2023-10-03 15:43:50 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 15:53:54 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41307"
# [1] "2023-10-03 15:53:54 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:02:32 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41307"
# [1] "2023-10-03 16:02:32 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:11:02 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41313"
# [1] "2023-10-03 16:11:02 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:19:09 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41311"
# [1] "2023-10-03 16:19:09 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:27:47 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41321"
# [1] "2023-10-03 16:27:47 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:37:00 ... ending   pcubature -- tol=0 -- ret.val is: 2991.41293"
# 9:     2991.4131: -1.92123  1.13224  8.11508  1.03111  2.27346 0.417036
# [1] "2023-10-03 16:37:00 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:45:48 ... ending   pcubature -- tol=0 -- ret.val is: 2990.88038"
# [1] "2023-10-03 16:45:48 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 16:53:55 ... ending   pcubature -- tol=0 -- ret.val is: 2990.88039"
# [1] "2023-10-03 16:53:55 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 17:02:11 ... ending   pcubature -- tol=0 -- ret.val is: 2990.88043"
# [1] "2023-10-03 17:02:11 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 17:10:42 ... ending   pcubature -- tol=0 -- ret.val is: 2990.8804"
# [1] "2023-10-03 17:10:42 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 17:19:02 ... ending   pcubature -- tol=0 -- ret.val is: 2990.88024"
# [1] "2023-10-03 17:19:02 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 17:27:56 ... ending   pcubature -- tol=0 -- ret.val is: 2990.88046"
# [1] "2023-10-03 17:27:56 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-03 17:37:09 ... ending   pcubature -- tol=0 -- ret.val is: 2990.88034"
# 10:     2990.8804: -1.93180  1.11665  8.17100  1.08979  2.40599 0.500515
# [1] "2023-10-03 17:37:09 ... starting pcubature for bivariate-t-corr"

## ... -- > Umm.  It stalled for two hours.  Might need to runs these on biowulf.
## Also: maybe lock down df...I think print out `8` was the best and
## maybe got overfit in 9/10 and maybe the reason is that df
## isn't gaining traction and the algo is spinning its wheels (?)

########################################$$$$#########
## 2023-10-02A                                     ##
## _END_ two random param: fix-df bivariate-t sim ###
#####################################################

########################################$$$$#########
## 2023-10-01A                                     ##
## START two random param: fix-df bivariate-t sim ###
#####################################################

## I'm starting to think that one can't estimate df
## but could still set it and test on a grid (?)

## tl;dr: false convergence?!? let's try not baby-ing the tolerances in 2023-10-02A

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t_df_gt_2"){
    print("t-BIVARIATE-t where shape is in terms for variance")
    a <- alpha
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c((a-2)/a*v1,(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
##detach(summed_binom_dat)
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t", alpha=8, v1=1, v2=2, rho=0.5)

data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=nu_set),

                   pmu = c(Intercept=-1.2, b_p=1, nu_set=8),
                   pmix=c(nu_set=8, v1=1, v2=1, corr= 0.30),

                   p_uppb = c(  0,   2,    8, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2,    8, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
# > (rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=nu_set),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, nu_set=8),
#                           +                    pmix=c(nu_set=8, v1=1, v2=1, corr= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2,    8, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2,    8, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-t-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.1,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  8   <? 8   <? 8     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-10-01 22:51:17 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:51:45 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 1.425969   log bounds ratio= 0.3467875
# Method:  nlminb
# [1] "2023-10-01 22:51:45 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:52:12 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# [1] "2023-10-01 22:52:12 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:52:40 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# [1] "2023-10-01 22:52:40 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:53:07 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# [1] "2023-10-01 22:53:07 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:53:34 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# [1] "2023-10-01 22:53:34 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:54:01 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# [1] "2023-10-01 22:54:01 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:54:28 ... ending   pcubature -- tol=0.1 -- ret.val is: 3087.95747"
# 0:     3087.9575: -1.20000  1.00000  8.00000  1.00000  1.00000 0.300000
# [1] "2023-10-01 22:54:28 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:54:55 ... ending   pcubature -- tol=0.1 -- ret.val is: 3019.35227"
# [1] "2023-10-01 22:54:55 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:55:22 ... ending   pcubature -- tol=0.1 -- ret.val is: 3019.35595"
# [1] "2023-10-01 22:55:22 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:55:49 ... ending   pcubature -- tol=0.1 -- ret.val is: 3019.35326"
# [1] "2023-10-01 22:55:49 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:56:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 3019.35221"
# [1] "2023-10-01 22:56:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:56:43 ... ending   pcubature -- tol=0.1 -- ret.val is: 3019.35381"
# [1] "2023-10-01 22:56:43 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:57:11 ... ending   pcubature -- tol=0.1 -- ret.val is: 3019.35329"
# 1:     3019.3523: -1.67087  1.02085  8.00000  1.10459  1.09639 0.387244
# [1] "2023-10-01 22:57:11 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:58:07 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.70614"
# [1] "2023-10-01 22:58:07 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:59:03 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.70607"
# [1] "2023-10-01 22:59:03 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 22:59:59 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.70635"
# [1] "2023-10-01 22:59:59 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:00:55 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.70609"
# [1] "2023-10-01 23:00:55 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:01:51 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.70666"
# [1] "2023-10-01 23:01:51 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:02:47 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.70614"
# 2:     3000.7061: -2.06836  1.17931  8.00000  1.04474  1.31134 0.518033
# [1] "2023-10-01 23:02:47 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:03:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.66054"
# [1] "2023-10-01 23:03:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:03:41 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.66038"
# [1] "2023-10-01 23:03:41 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:04:08 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.66032"
# [1] "2023-10-01 23:04:08 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:04:35 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.66025"
# [1] "2023-10-01 23:04:35 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:05:02 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.66074"
# [1] "2023-10-01 23:05:02 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:05:29 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.6607"
# 3:     3000.6605: -1.78792  1.33309  8.00000  1.18654  1.66748 0.545661
# [1] "2023-10-01 23:05:29 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:07:23 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.21184"
# [1] "2023-10-01 23:07:23 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:09:18 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.2118"
# [1] "2023-10-01 23:09:18 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:11:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.2119"
# [1] "2023-10-01 23:11:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:13:07 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.21187"
# [1] "2023-10-01 23:13:07 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:15:02 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.21204"
# [1] "2023-10-01 23:15:02 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:16:57 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.21198"
# 4:     2995.2118: -2.01752  1.28412  8.00000  1.14040  1.71718 0.492862
# [1] "2023-10-01 23:16:57 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:18:51 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.43234"
# 5:     2995.2118: -2.01752  1.28412  8.00000  1.14040  1.71718 0.492862
# [1] "2023-10-01 23:18:51 ... starting pcubature for bivariate-t-corr"
# [1] "2023-10-01 23:20:46 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.21184"
# Post processing for method  nlminb
# Save results from method  nlminb
# $par
# Intercept        b_p     nu_set         v1         v2       corr
# -2.0175213  1.2841241  8.0000000  1.1404001  1.7171831  0.4928617
#
# $message
# [1] "false convergence (8)"
#
# $convcode
# [1] 1
#
# $value
# [1] 2995.212
#
# $fevals
# function
# 6
#
# $gevals
# gradient
# 25
#
# $nitns
# [1] 5
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 1713.869
#
# Assemble the answers
# Intercept      b_p nu_set     v1       v2      corr    value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb -2.017521 1.284124      8 1.1404 1.717183 0.4928617 2995.212      6     25     5        1   NA   NA 1713.869
########################################$$$$#########
## 2023-10-01A                                     ##
## _END_ two random param: fix-df bivariate-t sim ###
#####################################################

########################################$$$$#########
## 2023-09-30F                                     ##
## START two random parameters: bivariate-t sim   ###
#####################################################

## try framing in terms of variance.

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t_df_gt_2"){
    print("t-BIVARIATE-t where shape is in terms for variance")
    a <- alpha
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c((a-2)/a*v1,(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
##detach(summed_binom_dat)
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t_df_gt_2", alpha=8, v1=1, v2=2, rho=0.5)

data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=degfree),

                   pmu = c(Intercept=-1.2, b_p=1, degfree=2.1),
                   pmix=c(degfree=2.1, v1=1, v2=1, corr= 0.30),

                   p_uppb = c(  0,   2,    200, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2,   2.01, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-varcorr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)

########################################$$$$#########
## 2023-09-30F                                     ##
## START two random parameters: bivariate-t sim   ###
#####################################################


########################################$$$$#########
## 2023-09-30E                                     ##
## START two random parameters: bivariate-t sim   ###
#####################################################

## try framing in terms of variance.

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t_df_gt_2"){
    print("t-BIVARIATE-t where shape is in terms for variance")
    a <- alpha
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c((a-2)/a*v1,(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*rho*sqrt(v1*v2),(a-2)/a*v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t_df_gt_2", alpha=8, v1=1, v2=1, rho=0.5)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=degfree),

                   pmu = c(Intercept=-1.2, b_p=1, degfree=40),
                   pmix=c(degfree=40, v1=1, v2=1, corr= 0.30),

                   p_uppb = c(  0,   2, 200, 1.00, 1.00, 0.90),
                   p_lowb = c( -4,  -2,   3, 1.00, 1.00,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-varcorr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
# > (rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=degfree),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, degfree=40),
#                           +                    pmix=c(degfree=40, v1=1, v2=1, v12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2, 200, 1.00, 1.00, 0.90),
#                           +                    p_lowb = c( -4,  -2,   1, 1.00, 1.00,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-t-varcorr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.1,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  1   <? 40   <? 200     In Bounds
# par[ 4 ]:  1   <? 1   <? 1     In Bounds
# par[ 5 ]:  1   <? 1   <? 1     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-09-30 20:58:53 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 20:59:36 ... ending   pcubature -- tol=0.1 -- ret.val is: 2998.06627"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 2.124939   log bounds ratio= 2.043581
# Method:  nlminb
# [1] "2023-09-30 20:59:36 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:00:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 2998.06627"
# [1] "2023-09-30 21:00:16 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:00:55 ... ending   pcubature -- tol=0.1 -- ret.val is: 2998.06628"
# [1] "2023-09-30 21:00:55 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:01:34 ... ending   pcubature -- tol=0.1 -- ret.val is: 2998.06627"
# [1] "2023-09-30 21:01:34 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:02:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 2998.06627"
# [1] "2023-09-30 21:02:16 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:02:57 ... ending   pcubature -- tol=0.1 -- ret.val is: 2998.06627"
# 0:     2998.0663: -1.20000  1.00000  40.0000  1.00000  1.00000 0.300000
# [1] "2023-09-30 21:02:57 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:05:42 ... ending   pcubature -- tol=0.1 -- ret.val is: 2940.39655"
# [1] "2023-09-30 21:05:42 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:08:29 ... ending   pcubature -- tol=0.1 -- ret.val is: 2940.39827"
# [1] "2023-09-30 21:08:29 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:11:14 ... ending   pcubature -- tol=0.1 -- ret.val is: 2940.39601"
# [1] "2023-09-30 21:11:14 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:13:56 ... ending   pcubature -- tol=0.1 -- ret.val is: 2940.39655"
# [1] "2023-09-30 21:13:56 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:16:40 ... ending   pcubature -- tol=0.1 -- ret.val is: 2940.39614"
# 1:     2940.3965: -1.69877 0.965164  40.0004  1.00000  1.00000 0.304061
# [1] "2023-09-30 21:16:40 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:19:42 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41495"
# [1] "2023-09-30 21:19:42 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:22:41 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41493"
# [1] "2023-09-30 21:22:41 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:25:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41478"
# [1] "2023-09-30 21:25:37 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:28:32 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41496"
# [1] "2023-09-30 21:28:32 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:31:28 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41494"
# [1] "2023-09-30 21:31:28 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:34:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41512"
# 2:     2936.4150: -1.88723  1.35420  39.9997  1.00000  1.00000 0.0527960
# [1] "2023-09-30 21:34:13 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:36:53 ... ending   pcubature -- tol=0.1 -- ret.val is: 2951.90746"
# 3:     2936.4150: -1.88723  1.35420  39.9997  1.00000  1.00000 0.0527960
# [1] "2023-09-30 21:36:53 ... starting pcubature for bivariate-t-varcorr"
# [1] "2023-09-30 21:39:32 ... ending   pcubature -- tol=0.1 -- ret.val is: 2936.41495"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept         b_p     degfree          v1          v2         v12
# -1.88723376  1.35419633 39.99974394  1.00000000  1.00000000  0.05279598
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2936.415
#
# $fevals
# function
# 4
#
# $gevals
# gradient
# 13
#
# $nitns
# [1] 3
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 2355.741
#
# Assemble the answers
# Intercept      b_p  degfree v1 v2        v12    value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb -1.887234 1.354196 39.99974  1  1 0.05279598 2936.415      4     13     3        0   NA   NA 2355.741

########################################$$$$#########
## 2023-09-30E                                     ##
## _END_ two random parameters: bivariate-t sim   ###
#####################################################

########################################$$$$#########
## 2023-09-30DD                                    ##
## START two random parameters: bivariate-t sim   ###
#####################################################

##tl;dr: locking corners to 1 still did not get df mixing

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=,alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t", alpha=8, v1=1, v2=1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=degfree),

                   pmu = c(Intercept=-1.2, b_p=1, degfree=40),
                   pmix=c(degfree=40, var1=1, var2=1, corr12= 0.30),

                   p_uppb = c(  0,   2, 200, 1.00, 1.00, 0.90),
                   p_lowb = c( -4,  -2,   1, 1.00, 1.00,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
# > (rand.int.rand.slopes.nonzero.corr.estimate.degfree.corr.mat <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=degfree),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, degfree=40),
#                           +                    pmix=c(degfree=40, var1=1, var2=1, corr12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2, 200, 1.00, 1.00, 0.90),
#                           +                    p_lowb = c( -4,  -2,   1, 1.00, 1.00,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-t-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.1,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  1   <? 40   <? 200     In Bounds
# par[ 4 ]:  1   <? 1   <? 1     In Bounds
# par[ 5 ]:  1   <? 1   <? 1     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-09-30 18:58:19 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 18:58:58 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.73074"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 2.124939   log bounds ratio= 2.043581
# Method:  nlminb
# [1] "2023-09-30 18:58:58 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 18:59:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.73074"
# [1] "2023-09-30 18:59:37 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:00:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.73075"
# [1] "2023-09-30 19:00:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:00:54 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.73074"
# [1] "2023-09-30 19:00:54 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:01:33 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.73074"
# [1] "2023-09-30 19:01:33 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:02:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 3000.73074"
# 0:     3000.7307: -1.20000  1.00000  40.0000  1.00000  1.00000 0.300000
# [1] "2023-09-30 19:02:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:04:58 ... ending   pcubature -- tol=0.1 -- ret.val is: 2948.78076"
# [1] "2023-09-30 19:04:58 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:07:39 ... ending   pcubature -- tol=0.1 -- ret.val is: 2948.78211"
# [1] "2023-09-30 19:07:39 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:10:23 ... ending   pcubature -- tol=0.1 -- ret.val is: 2948.77992"
# [1] "2023-09-30 19:10:23 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:13:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 2948.78076"
# [1] "2023-09-30 19:13:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:16:01 ... ending   pcubature -- tol=0.1 -- ret.val is: 2948.78153"
# 1:     2948.7808: -1.69552 0.998636  40.0004  1.00000  1.00000 0.366791
# [1] "2023-09-30 19:16:01 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:18:44 ... ending   pcubature -- tol=0.1 -- ret.val is: 2950.15986"
# 2:     2948.7808: -1.69552 0.998636  40.0004  1.00000  1.00000 0.366791
# [1] "2023-09-30 19:18:44 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 19:21:24 ... ending   pcubature -- tol=0.1 -- ret.val is: 2948.78076"
# Post processing for method  nlminb
# Save results from method  nlminb
# $par
# Intercept        b_p    degfree       var1       var2     corr12
# -1.6955168  0.9986360 40.0004283  1.0000000  1.0000000  0.3667915
#
# $message
# [1] "false convergence (8)"
#
# $convcode
# [1] 1
#
# $value
# [1] 2948.781
#
# $fevals
# function
# 3
#
# $gevals
# gradient
# 8
#
# $nitns
# [1] 2
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 1322.054
#
# Assemble the answers
# Intercept      b_p  degfree var1 var2    corr12    value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb -1.695517 0.998636 40.00043    1    1 0.3667915 2948.781      3      8     2        1   NA   NA 1322.054
########################################$$$$#########
## 2023-09-30DD                                   ##
## _END_ two random parameters: bivariate-t sim   ###
#####################################################

########################################$$$$#########
## 2023-09-30A                                     ##
## START two random parameters: bivariate-t sim   ###
#####################################################

## tl;dr: can't have v1,v2, and df.  will have to do
## the trick where you lock v1,v2 and estimate
## df.  See your paper


library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="t-BIVARIATE_t"){
    print("t-BIVARIATE-t could be fun")
    G <- function(n, v1, v2, rho){mvtnorm::rmvt(n=n, df=,alpha,sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pt(x, df=alpha)
  }
  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
set.seed(166)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="t-BIVARIATE_t", alpha=8)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


(rand.int.rand.slopes.nonzero.corr.estimate.degfree <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=8),

                   pmu = c(Intercept=-1.2, b_p=1),#, degfree=40),
                   pmix=c(degfree=40, var1=1, var2=1.3, corr12= 0.30),

                   p_uppb = c(  0,   2, 200, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2,  1, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
# > (rand.int.rand.slopes.nonzero.corr.estimate.degfree <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=degfree),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, degfree=10),
#                           +                    pmix=c(degfree=10, var1=1, var2=1.3, corr12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2, 20, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2,  1, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-t-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.1,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  1   <? 10   <? 20     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-09-30 12:30:02 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:30:31 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 1.522879   log bounds ratio= 1.023481
# Method:  nlminb
# [1] "2023-09-30 12:30:32 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:31:00 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# [1] "2023-09-30 12:31:00 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:31:29 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# [1] "2023-09-30 12:31:29 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:31:58 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# [1] "2023-09-30 12:31:58 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:32:26 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# [1] "2023-09-30 12:32:26 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:32:55 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# [1] "2023-09-30 12:32:55 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:33:24 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# [1] "2023-09-30 12:33:24 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:33:52 ... ending   pcubature -- tol=0.1 -- ret.val is: 3070.54475"
# 0:     3070.5447: -1.20000  1.00000  10.0000  1.00000  1.30000 0.300000
# [1] "2023-09-30 12:33:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:34:21 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.19796"
# [1] "2023-09-30 12:34:21 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:34:49 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.20096"
# [1] "2023-09-30 12:34:49 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:35:18 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.19872"
# [1] "2023-09-30 12:35:18 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:35:47 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.19803"
# [1] "2023-09-30 12:35:47 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:36:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.19762"
# [1] "2023-09-30 12:36:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:36:44 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.19905"
# [1] "2023-09-30 12:36:44 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:37:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 3008.19877"
# 1:     3008.1980: -1.68066  1.00799  10.0075  1.08270  1.37132 0.383195
# [1] "2023-09-30 12:37:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:38:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.863"
# [1] "2023-09-30 12:38:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:39:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.86289"
# [1] "2023-09-30 12:39:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:40:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.86291"
# [1] "2023-09-30 12:40:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:41:12 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.86297"
# [1] "2023-09-30 12:41:12 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:42:12 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.86283"
# [1] "2023-09-30 12:42:12 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:43:11 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.86338"
# [1] "2023-09-30 12:43:11 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:44:11 ... ending   pcubature -- tol=0.1 -- ret.val is: 2997.86326"
# 2:     2997.8630: -2.01866  1.20247  10.0178 0.910451  1.59362 0.520091
# [1] "2023-09-30 12:44:11 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:46:19 ... ending   pcubature -- tol=0.1 -- ret.val is: 3004.73148"
# [1] "2023-09-30 12:46:19 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:48:22 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64787"
# [1] "2023-09-30 12:48:22 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:50:25 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64784"
# [1] "2023-09-30 12:50:25 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:52:28 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64795"
# [1] "2023-09-30 12:52:28 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:54:30 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64786"
# [1] "2023-09-30 12:54:30 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:56:32 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64786"
# [1] "2023-09-30 12:56:32 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:58:35 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64811"
# [1] "2023-09-30 12:58:35 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:00:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 2995.64734"
# 3:     2995.6479: -1.87396  1.18979  10.0157 0.979392  1.63565 0.570909
# [1] "2023-09-30 13:00:37 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:02:38 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.84132"
# [1] "2023-09-30 13:02:38 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:04:39 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.84131"
# [1] "2023-09-30 13:04:39 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:06:40 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.84128"
# [1] "2023-09-30 13:06:40 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:08:42 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.84133"
# [1] "2023-09-30 13:08:42 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:10:52 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.84131"
# [1] "2023-09-30 13:10:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:12:57 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.84177"
# [1] "2023-09-30 13:12:57 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:15:06 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.8413"
# 4:     2993.8413: -1.91147  1.14085  10.0170 0.990261  1.74902 0.455024
# [1] "2023-09-30 13:15:06 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:17:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69523"
# [1] "2023-09-30 13:17:13 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:19:21 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69519"
# [1] "2023-09-30 13:19:21 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:21:26 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69519"
# [1] "2023-09-30 13:21:26 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:23:29 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69523"
# [1] "2023-09-30 13:23:29 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:25:33 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69515"
# [1] "2023-09-30 13:25:33 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:27:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69555"
# [1] "2023-09-30 13:27:37 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:29:41 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.6952"
# 5:     2992.6952: -1.88295  1.15998  10.0200  1.01825  1.91039 0.501852
# [1] "2023-09-30 13:29:41 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:31:44 ... ending   pcubature -- tol=0.1 -- ret.val is: 2993.13413"
# 6:     2992.6952: -1.88295  1.15998  10.0200  1.01825  1.91039 0.501852
# [1] "2023-09-30 13:31:44 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 13:33:47 ... ending   pcubature -- tol=0.1 -- ret.val is: 2992.69523"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p    degfree       var1       var2     corr12
# -1.8829471  1.1599813 10.0200049  1.0182483  1.9103897  0.5018522
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2992.695
#
# $fevals
# function
# 8
#
# $gevals
# gradient
# 36
#
# $nitns
# [1] 6
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 3725.777
#
# Assemble the answers
# Intercept      b_p degfree     var1    var2    corr12    value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb -1.882947 1.159981   10.02 1.018248 1.91039 0.5018522 2992.695      8     36     6        0   NA   NA 3725.777
########################################$$$$#########
## 2023-09-30B                                     ##
## _END_ two random parameters: bivariate-t sim   ###
#####################################################

########################################$$$$#########
## 2023-09-30A                                     ##
## START two random parameters: glmer == gnlrim.  ###
#####################################################

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
set.seed(1709001)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1, mrim="Probit-BIVARIATE_NORM")
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat





attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)


## takes about Yhr XXmin to run:
(rand.int.rand.slopes.nonzero.corr.CUBA.CONDITIONAL <-
  gnlrim::gnlrem(y=ybind,
                 mu = ~ pnorm(
                   (Intercept + period_numeric*b_p)+

                     rand1 + period_numeric*rand2
                 ),
                 pmu = c(Intercept=-0.95, b_p=0.55),
                 pmix=c(var1=1, var2=1, corr12= 0.20),
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
                 method='nlminb',
                 int2dmethod="cuba"
  )
)



########################################$$$$#########
## 2023-09-30A                                     ##
## _END_ two random parameters: glmer == gnlrim.  ###
#####################################################





























####################################################
## BEGIN malcolm idea                      ##########
####################################################
ab_dat <-
structure(list(monkey_id = c("DC4C", "DC4C", "DC4C", "DC4C",
                             "DC4C", "DCG3", "DCG3", "DCG3", "DECD", "DECD", "DECD", "DECD",
                             "DEV2", "DEV2", "DEV2", "DEV2", "DEV2", "DEV2", "DEWH", "DEWH",
                             "DEWH", "DEWH", "DEWH", "DEWH", "DEXW", "DEXW", "DEXW", "DEXW",
                             "DEXi", "DEXi", "DEXi", "DEXi", "DF02", "DF02", "DF02", "DF02",
                             "DF02", "DF02", "DF02", "DF87", "DF87", "DF87", "DF87", "DFCL",
                             "DFCL", "DFCL", "DFCL", "DFCL", "J7ZA", "J7ZA", "J7ZA", "KZIA",
                             "KZIA", "KZIA", "KZIA", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M",
                             "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M",
                             "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M",
                             "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M", "DF3M",
                             "DF3M", "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z",
                             "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z", "DF6Z",
                             "DF6Z", "DF6Z", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V",
                             "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V",
                             "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V",
                             "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V",
                             "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V", "DF9V",
                             "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0",
                             "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0",
                             "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFB0", "DFC6", "DFC6",
                             "DFC6", "DFC6", "DFC6", "DFC6", "DFC6", "DFC6", "DFC6", "DFC6",
                             "DFC6", "DFC6", "DFC6", "DFC6", "DFC6", "DFC6", "DFC6", "DFC6",
                             "DFC6", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7",
                             "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7",
                             "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7", "DFD7",
                             "DFD7", "DFD7", "DFH3", "DFH3", "DFH3", "DFH3", "DFH3", "DFH3",
                             "DFH3", "DFH3", "DFH3", "DFH3", "DFH3", "DFH3", "DFJ4", "DFJ4",
                             "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4",
                             "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFJ4",
                             "DFJ4", "DFJ4", "DFJ4", "DFJ4", "DFL4", "DFL4", "DFL4", "DFL4",
                             "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4",
                             "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4",
                             "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4", "DFL4",
                             "DFL4", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5",
                             "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5",
                             "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5", "DFL5",
                             "DFL5", "DFL5", "DFL5", "DFL5", "DFM1", "DFM1", "DFM1", "DFM1",
                             "DFM1", "DFM1", "DFM1", "DFM1", "DFM1", "DFM1", "DFM1", "DFM1",
                             "DFM1", "DFM1", "DFM1", "DFM1", "DFM1", "DFM1", "DFM1", "DFM1",
                             "DFM1", "DFM1", "DFM1", "DFM6", "DFM6", "DFM6", "DFM6", "DFM6",
                             "DFM6", "DFM6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN6",
                             "DFN6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN6",
                             "DFN6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN6", "DFN9", "DFN9",
                             "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9",
                             "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9",
                             "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFN9", "DFP9", "DFP9",
                             "DFP9", "DFP9", "DFP9", "DFP9", "DFP9", "DFP9", "DFP9", "DFP9",
                             "DFP9", "DFP9", "DFP9", "DFP9", "DFP9", "DFP9", "DFP9", "DFX9",
                             "DFX9", "DFX9", "DFX9", "DFX9", "DFX9", "DFX9", "DFX9", "DFX9",
                             "DFX9", "DFX9", "DFX9", "DFX9", "DFX9", "DFX9", "DFX9", "DFX9",
                             "DG43", "DG43", "DG43", "DG43", "DG43", "DG43", "DG43", "DG43",
                             "DG43", "DG43", "DG43", "DG43", "DG43", "DG43", "DG43", "DG43",
                             "DG43", "DG43", "DGGV", "DGGV", "DGGV", "DGGV", "DGGV", "DGGV",
                             "DGGV", "DGGV", "DGGV", "DGGV", "DGGV", "DGGV", "DGGV", "DGGV",
                             "DGGV", "DGGV", "DGGV", "DGGV", "DGGV", "DGGV", "DF60", "DF60",
                             "DF60", "DF60", "DF60", "DF60", "DF60", "DF60", "DF60", "DF60",
                             "DF60", "DF60", "DF60", "DF60", "DF80", "DF80", "DF80", "DF80",
                             "DF80", "DF80", "DF80", "DFAM", "DFAM", "DFAM", "DFAM", "DFAM",
                             "DFAM", "DFAM", "DFAM", "DFAM", "DFAM", "DFAM", "DFAM", "DFAM",
                             "DFAM", "DFAM", "DFAM", "DFAM", "DFAM", "DFCP", "DFCP", "DFCP",
                             "DFCP", "DFCP", "DFCP", "DFCP", "DFCP", "DFCP", "DFCP", "DFDP",
                             "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP",
                             "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP",
                             "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFDP", "DFFN",
                             "DFFN", "DFFN", "DFFN", "DFFN", "DFFN", "DFFN", "DFFN", "DFFN",
                             "DFFN", "DFFN", "DFFN", "DFFN", "DFGP", "DFGP", "DFGP", "DFGP",
                             "DFGP", "DFGP", "DFGP", "DFGP", "DFIP", "DFIP", "DFIP", "DFIP",
                             "DFIP", "DFIP", "DFIP", "DFIP", "DFKV", "DFKV", "DFKV", "DFKV",
                             "DFKV", "DFKV", "DFKV", "DFKV", "DFKV", "DFKV", "DFKV", "DFKV",
                             "DFKV", "M57", "M57", "M57", "M57", "M57", "M57", "M57", "M57",
                             "M57", "M57", "M57", "M57", "M57", "M57", "M57", "M57", "M57",
                             "M57", "M57", "M57", "M57", "M57", "MMK", "MMK", "MMK", "MMK",
                             "MMK", "MMK", "MMK", "MMK", "MMK", "MMK", "MMK", "MMK", "MMK",
                             "MMK", "MMK", "MMK", "MMK", "MRF", "MRF", "MRF", "MRF", "MRF",
                             "MRF", "MRF", "MRF", "MRF", "MRF", "MRF", "MRF", "MRF", "MRF",
                             "MRF", "MRF"), ab_1010 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 176.874666666667, 84.7157713105304, 63.16, 60.9328166513353,
                                                        54.0693333333333, 45.7952014942986, 42.5813333333333, 33.6703651802256,
                                                        34.1146666666667, 24.7557266782829, 19.112, 18.201347092302,
                                                        15.6069291123939, 13.6373333333333, 12.2565407444802, 10.5094947968534,
                                                        5.704, 7.72697828862676, 6.62557566712235, 6.49066666666667,
                                                        4.87137396412591, 4.13866666666667, 3.581618487299, 3.07109498848885,
                                                        2.63334145213048, 2.25798525591058, 1.93613228994085, 1.69703118037814,
                                                        1.76266666666667, 1.33272418021303, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 110.565333333333, 64.8327283431107, 58.5466666666667,
                                                        48.0281364433086, 47.3706666666667, 37.0314156613325, 41.752,
                                                        27.9871272959603, 23.208, 21.1517512979709, 14.9173333333333,
                                                        15.9857986938084, 13.8972331564224, 9.08533333333333, 11.1525984222099,
                                                        9.69549683078581, 6.38666666666667, 7.32753796080553, 6.37018552869003,
                                                        5.24533333333333, 4.81437693121983, 3.92266666666667, 3.63854790907921,
                                                        3.16316685902419, 2.74989496580822, 2.39061758674033, 2.0782802678256,
                                                        1.84325189561268, 1.70133333333333, 1.3654824623971, 1.18708039018958,
                                                        1.77866666666667, 0.897156357580597, 0.68, 0.678041299139597,
                                                        0.589454315315267, 0.493333333333333, 0.454490417082317, 297.837333333333,
                                                        87.6440111687157, 78.608, 56.9723378778518, 34.2053333333333,
                                                        39.2235529902573, 26.0613333333333, 26.2397246562657, 14.3893333333333,
                                                        17.5538190078717, 14.44, 11.743132437464, 9.60483990024304, 8.55733333333333,
                                                        7.00348356238791, 5.72822785734048, 4.48266666666667, 3.83206288788524,
                                                        3.13428727145902, 2.82666666666667, 2.09677167739371, 2.14666666666667,
                                                        92.7706666666667, 63.4805764331769, 46.328, 37.0635650203169,
                                                        31.472, 23.2494495537988, 29.2693333333333, 14.0701452842328,
                                                        9.04, 8.51499678998081, 5.59466666666667, 5.1531216535933, 4.00879039023359,
                                                        3.41066666666667, 2.70170539268819, 2.10174945276102, 1.664,
                                                        1.27194065747558, 1.10191612236524, 35.6453333333333, 25.2753895026207,
                                                        19.5306666666667, 19.3539645908988, 11.4453333333333, 14.2235159313132,
                                                        12.2026666666667, 10.6699286878119, 7.68, 8.00416568960664, 6.93255536730842,
                                                        4.536, 5.30842189749124, 3.51466666666667, 3.98217173338668,
                                                        4.581, 2.92655713122911, 2.847, 2.19538938486209, 2.022, 1.64689576695226,
                                                        1.42640676509327, 1.263, 1.209, 1.00610972641938, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 45.864, 25.011273498022, 21.4266666666667,
                                                        16.5709050739645, 9.05066666666667, 10.3050609579135, 4.856,
                                                        6.61466638357424, 4.49066666666667, 4.24585662760079, 3.4016843733834,
                                                        2.584, 2.25374325091373, 1.08266666666667, 1.44664449632182,
                                                        1.446, 0.899634609726228, 0.999, 0.577462164926475, 0.579, 0.37066443233331,
                                                        0.337071304474315, 112.754666666667, 48.8079724091703, 52.096,
                                                        34.7374374027862, 23.0533333333333, 25.8700451371249, 16.5786666666667,
                                                        18.8343395815749, 12.0133333333333, 13.7120884634645, 10.904,
                                                        9.98290219922624, 8.51791734322212, 6.71466666666667, 6.63781630041499,
                                                        5.66372077558866, 3.864, 4.1233959901205, 3.51829012408259, 2.58933333333333,
                                                        2.56144396317188, 2.62133333333333, 1.86482494196822, 1.59116301034136,
                                                        1.35766080155836, 1.15842490059686, 0.988426747522294, 0.862714889383478,
                                                        1.06133333333333, 115.192, 67.5224199453182, 58.9146666666667,
                                                        42.3557184607401, 34.3493333333333, 28.2735863514763, 20.9946666666667,
                                                        18.2956292308792, 11.6826666666667, 11.8389667583265, 11.208,
                                                        7.66090808553322, 6.16259286927212, 5.888, 4.37761049465037,
                                                        3.52144039813335, 2.69333333333333, 2.27869811356717, 1.83303187939442,
                                                        1.368, 1.18614141187756, 1.04266666666667, 0.767543360694681,
                                                        0.617427749904269, 0.496671648630795, 0.399532943881904, 0.321392561236367,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 49.0506666666667, 25.7652994914264, 21.776, 5.1579757664245,
                                                        4.31466666666667, 0.806219432152068, 0.1, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 34.7066666666667, 20.7415678062884,
                                                        14.6746666666667, 15.0001549066231, 12.04, 10.3204040820631,
                                                        5.12533333333333, 7.27988070308962, 3.83466666666667, 5.13513449956142,
                                                        4.31286241053437, 2.86133333333333, 3.11903154347638, 2.632,
                                                        2.20012485607468, 2.235, 1.51372953726337, 1.224, 1.06776540534617,
                                                        0.942, 0.753188025197207, 0.632582519934986, 0.729, 0.446215598087893,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 23.1386666666667,
                                                        19.816202765992, 11.4426666666667, 14.0357081536252, 8.02133333333333,
                                                        9.42766635616912, 8.24533333333333, 6.50273548012195, 5.90933333333333,
                                                        4.48526359832057, 3.72506549756222, 2.66133333333333, 2.63844354008378,
                                                        1.85333333333333, 1.81986716247912, 2.058, 1.22238936807545,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 34, 21.2728434827817,
                                                        27.088, 14.4133368991981, 7.49866666666667, 9.19801568076967,
                                                        3.89866666666667, 6.04823206194949, 2.688, 3.9770655263911, 3.22500132001158,
                                                        2.76533333333333, 2.18508778872504, 1.32, 1.43682273227384, 1.068,
                                                        0.916922855156214, 0.885, 0.602930284461749, 0.882, 123.17, 51.46,
                                                        13.14, 0.82, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                                        118.36, 54.04, 19.31, 0.48, 0.1, 0.1, 0.1, 112.15, 65.71, 40.26,
                                                        29.05, 19.22, 10.25, 11.17, 8.04, 5.69, 3.78, 2.51, 1.71, 0.72,
                                                        0.91, 0.61, 0.45, 0.1, 0.143361912184393, 157.09, 75.66, 30.14,
                                                        9.41, 1.07, 0.1, 0.1, 0.1, 0.1, 0.1, 165.2, 121.26, 115.26, 66.49,
                                                        48.57, 43.68, 34.49, 18.74, 14.52, 13.05, 9.75, 6.48, 5.44, 4.42,
                                                        4.21, 3.72, 2.50284775349143, 1.95, 1.8, 1.09, 0.93, 0.78, 0.62,
                                                        0.46, 105.03, 70.47, 34.16, 26.9, 14.98, 7.47, 5.95, 3.19, 2.43,
                                                        1.05, 0.71, 0.65, 0.27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                        0, 0, 0, 0, 0, 0, 0, 0), ab_3bnc = c(0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 46.0156666666667, 28.0666144340487,
                                                                                             16.1266666666667, 16.4399025787169, 10.127, 10.3414290304201,
                                                                                             11.0153333333333, 6.2773424217952, 6.01333333333333, 3.81040451610282,
                                                                                             2.69233333333333, 2.31295054511055, 1.80203814077535, 1.312,
                                                                                             1.21735543823838, 0.948451291024466, 0.82, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.553, 7.29577495450582, 4.76966666666667,
                                                                                             5.14394612208937, 4.838, 3.43693556520261, 2.132, 2.35896431189561,
                                                                                             1.66733333333333, 1.61909134437587, 1.34136305758526, 1.24366666666667,
                                                                                             0.945739053274708, 0.779, 0.649114489555307, 0.84075, 0.43370685113918,
                                                                                             0.28, 0.297677673687092, 0.246616002203204, 0.204313114212999,
                                                                                             0.169266585568191, 0.140231708083332, 0.119342797970361, 0.107175670702966,
                                                                                             76.178, 36.0931730272363, 29.4653333333333, 7.32278990065103,
                                                                                             5.19333333333333, 1.83778435299725, 0.28, 0.414695471295901,
                                                                                             0.28, 0.0935759049383894, 0.0444509868617954, 0.021115374030189,
                                                                                             11.4526666666667, 5.28425971399422, 3.74466666666667, 1.13766393185544,
                                                                                             0.28, 0.193389562140827, 0.0845853691523291, 0.0369962297625233,
                                                                                             0.0161815338794171, 0.00707753304516357, 0.00309559491570195,
                                                                                             0.00135396158816497, 0.000666459802099168, 0.000291498402301028,
                                                                                             0.000127496539590858, 0.0000495513189855696, 0.0000216729205132958,
                                                                                             0.000010668050292164, 0.00000414612008041526, 0.00000181344377467708,
                                                                                             0.000000793170062644616, 0.000000556480719662418, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38.3896666666667, 30.1095399786788,
                                                                                             19.7483333333333, 17.587527926431, 10.8513333333333, 11.0366869067125,
                                                                                             10.5233333333333, 6.68199311126329, 6.58733333333333, 4.04551042503655,
                                                                                             3.53966666666667, 2.44929234834025, 1.90578663199964, 1.99533333333333,
                                                                                             1.28481749439918, 0.999712430018614, 0.779, 0.605260584729084,
                                                                                             0.470951347247148, 0.28, 0.285130282723854, 0.221858971454196,
                                                                                             0.172627764208349, 11.1246666666667, 4.47499290163423, 2.30966666666667,
                                                                                             0.972443332304158, 0.28, 0.167088794042389, 0.0734490305558706,
                                                                                             40.1936666666667, 17.6079321077649, 14.0766666666667, 9.70530450454661,
                                                                                             5.48033333333333, 5.79166191336373, 5.05666666666667, 3.32162519218938,
                                                                                             1.99533333333333, 1.90501346287652, 0.902, 1.09256044368714,
                                                                                             0.827406988121237, 0.464666666666667, 0.534572293591573, 0.404836962503397,
                                                                                             0.28, 0.232181483224685, 0.175833367245236, 0.28, 4.19566666666667,
                                                                                             3.71337028897751, 1.89966666666667, 1.88018449184731, 2.20033333333333,
                                                                                             0.857355041850858, 0.478333333333333, 0.411961783916307, 0.28,
                                                                                             0.197948927950702, 0.137214897654948, 0.0951150801036256, 0.0694757868309133,
                                                                                             0.0481594574832702, 0.033383333257205, 0.0219604798705237, 0.0152226386325159,
                                                                                             0.011119212594727, 0.00731452555924368, 0.00507030720699877,
                                                                                             0.0035146524494477, 0.00243629849949874, 0.00168880151424148,
                                                                                             0.00117064906253939, 33.5243333333333, 28.2996294168237, 20.3496666666667,
                                                                                             15.1925180532303, 9.07466666666667, 8.8613299013877, 8.84233333333333,
                                                                                             4.95858752878273, 4.31866666666667, 2.77470656822619, 2.63766666666667,
                                                                                             1.55265919882785, 1.16146380381363, 1.148, 0.736029050227363,
                                                                                             0.550585151615869, 0.28, 7.14766666666667, 6.00376699089126,
                                                                                             4.00433333333333, 2.27460600733218, 2.26866666666667, 0.742231986131856,
                                                                                             0.28, 0.260974205934907, 0.154748534261711, 0.0917604434137871,
                                                                                             0.0544107187552092, 0.0322636443909517, 0.02061421569685, 0.012223505585981,
                                                                                             0.00724810931484245, 0.00398868058878366, 0.00236514743877165,
                                                                                             22.9736666666667, 10.5859345966271, 9.17033333333333, 6.09718241198231,
                                                                                             4.12733333333333, 3.77986108330756, 1.886, 2.2586511305208, 1.09333333333333,
                                                                                             1.34965407906971, 0.751666666666667, 0.806484059682775, 0.623422462176537,
                                                                                             0.505666666666667, 0.415983175596057, 0.32156029922788, 0.341666666666667,
                                                                                             0.28, 8.13166666666667, 5.5197925670622, 2.71966666666667, 2.06968431108821,
                                                                                             1.98166666666667, 0.66733534425447, 0.28, 0.232035701775232,
                                                                                             0.136823329495751, 0.0806799270589715, 0.0475741282881373, 0.0280527978256826,
                                                                                             0.0178382476760464, 0.0105185901166656, 0.00620244432366449,
                                                                                             0.00339154474862207, 0.00199987519631704, 0.00127168310607597,
                                                                                             0.000695366203267293, 0.000410033104483262, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 72.6, 47.4,
                                                                                             19.6, 12.2, 7.08359222212496, 3.72030338068268, 2.29520753084227,
                                                                                             1.3, 78.1, 39, 29.2, 16.8, 5.47326978675105, 2.48414569655377,
                                                                                             1.37364566502148, 0.5, 65.5, 30.8, 23.2, 8.6, 3.47304158724727,
                                                                                             1.52998769865642, 0.827315965738102, 0.1, 0.2, 0.0707302627623429,
                                                                                             0.1, 0.0206810318982027, 0.00911066671817492, 63, 37.9, 23.4,
                                                                                             11, 10.1912187122521, 6.53042880587955, 4.67711859333549, 3.16850584098524,
                                                                                             2.1464987607248, 1.23062324328079, 0.7, 0.631245767602256, 0.404495837125289,
                                                                                             0.3, 0.2, 0.2, 0.0952225225408209, 0.0577157982314732, 0.0437010570820693,
                                                                                             0.0296052049693737, 0.0200559945182252, 0.0169731536882421, 88.7,
                                                                                             54.6, 38.8, 20.7, 13.1587525205418, 7.8114486685762, 5.28286157660654,
                                                                                             3.34731425777582, 2.12091734334378, 1.10514642613787, 0.7, 0.505469691246663,
                                                                                             0.300062679994212, 0.2, 0.1, 0.1, 0.1, 64.3, 35.1, 29.1, 14.7,
                                                                                             12.62655805643, 8.43722719132248, 6.23570866231284, 4.3821398624588,
                                                                                             3.07954569625899, 1.86050185517039, 1.3, 1.0162539647212, 0.679073865269921,
                                                                                             0.5, 0.4, 0.3), infected = c(0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
                                                                                                                          1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                                                                                                                          0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
                                                                                                                          0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                                                          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)), row.names = c(NA, -612L
                                                                                                                          ), class = c("data.table", "data.frame"))



## this fails:
library(lme4)
m1<-glmer(infected ~  ab_1010 + ab_3bnc + (1|monkey_id),
          data=ab_dat,
          family=binomial(link="probit"),
          nAGQ=100)
summary(m1)



library(gnlrim)
id <- as.numeric(as.factor(ab_dat$monkey_id))
o.value <- with(ab_dat, cbind(infected,1-infected))
attach(ab_dat)
## if you have NA's in the data it may complain about initial values.
## if your ids are not sequential you'll get a likelihood is NA INF or probs too small error
(rand.int.probit <-
    gnlrim::gnlrim(y=o.value,
                   mu = ~ pnorm(beta0 + rand1),
                   pmu = c(beta0=0),
                   pmix=c(var1=1),
                   p_uppb = c(  1,  4),
                   p_lowb = c( -2,  0.0005),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)




(rand.int.probit <-
    gnlrim::gnlrim(y=o.value,
                   mu = ~ pnorm(beta0 + ab_1010*beta_ab_1010 + rand1),
                   pmu = c(beta0=0, beta_ab_1010=0),
                   pmix=c(var1=1),
                   p_uppb = c(  1,  1, 4),
                   p_lowb = c( -2, -2, 0.0005),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)






(rand.int.probit <-
    gnlrim::gnlrim(y=o.value,
                   mu = ~ pnorm(beta0 + ab_1010*beta_ab_1010 + rand1),
                   pmu = c(beta0=0, beta_ab_1010=0),
                   pmix=c(var1=1),
                   p_uppb = c(  1,  1, 4),
                   p_lowb = c( -2, -2, 0.0005),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)






(rand.int.probit <-
    gnlrim::gnlrim(y=o.value,
                   mu = ~ stable_cdf2(beta0 + ab_1010*beta_ab_1010 +  ab_3bnc*beta_ab_3bnc + rand1, c(alp,0,1,0)),
                   pmu = c(beta0=0, beta_ab_1010=0, beta_ab_3bnc=0, alp=1.5),
                   pmix=c(alp=1.5, scl=1),
                   p_uppb = c(  1,  1, 1, 1.99, 4),
                   p_lowb = c( -2, -2,-2, 1.40, 0.0005),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-scl",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)




####################################################
## BEGIN Sameer XDR COP idea              ##########
####################################################
## none of these worked very well.  SO i copied and pasted
## to create chunk above this one and fiddled with how many outliers
## I had to see if that would give comp. advantage to stable

## code copied from:
## /Volumes/swihartbj$/projects/_XDR/_12_costs_cdf_log10_nofacet_2curves.R
## which supports this paper:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6710115/
## and improve this figure/analysis;
## https://pubmed.ncbi.nlm.nih.gov/30824387/#&gid=article-figures&pid=figure-4-uid-3
library(ggplot2)
library(data.table)
dtep.money <- readRDS("dtep.money.RDS")
mx=max(dtep.money$total_direct_cost,na.rm=TRUE)
mn=min(dtep.money$total_direct_cost,na.rm=TRUE)

ggplot(dtep.money, aes(total_direct_cost_adj, colour=factor(cc1_or_ncc0)))+
  stat_ecdf() +
  coord_cartesian(xlim=c(5e3, 1.2*mx))+
  scale_x_log10(
    breaks=c(5e3,1e4,73187,101190,1e6,1.2*mx),
    labels=c("$5,000","$10,000","$73,187","$101,190","$1,000,000","$2,469,965")

  )+
  xlab("Total Direct Cost (log10 USD)") + ylab("Percentile") +
  theme(legend.position="bottom")+
  ggtitle("Cumulative log10 hospitalization costs for the CC and NCC groups")+
  scale_color_discrete(name="Group",
                       breaks=c("0","1"),
                       labels=c("Non-Colistin case", "Colistin case")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



ggplot(dtep.money) +
  geom_jitter(aes(x=death_flag, y=total_direct_cost,color=group))

dtep.money[,mean(total_direct_cost,na.rm=TRUE),by=group]

dtep.money[,{j=list(min=min(total_direct_cost,na.rm=TRUE),
                    q10=quantile(total_direct_cost,0.10,na.rm=TRUE),
                    q25=quantile(total_direct_cost,0.25,na.rm=TRUE),
                    q50=quantile(total_direct_cost,0.50,na.rm=TRUE),
                    q75=quantile(total_direct_cost,0.75,na.rm=TRUE),
                    q90=quantile(total_direct_cost,0.90,na.rm=TRUE),
                    q99=quantile(total_direct_cost,0.99,na.rm=TRUE),
                    max=max(total_direct_cost,na.rm=TRUE))
},
by=group]


dtep.money[,{j=list(min=min(total_direct_cost,na.rm=TRUE),
                    q10=quantile(total_direct_cost,0.10,na.rm=TRUE),
                    q25=quantile(total_direct_cost,0.25,na.rm=TRUE),
                    q50=quantile(total_direct_cost,0.50,na.rm=TRUE),
                    q75=quantile(total_direct_cost,0.75,na.rm=TRUE),
                    q90=quantile(total_direct_cost,0.90,na.rm=TRUE),
                    q99=quantile(total_direct_cost,0.99,na.rm=TRUE),
                    max=max(total_direct_cost,na.rm=TRUE))
},
by=c("group","death_flag")]


dtep.money[,mean(death_flag),by=group]

table(dtep.money[,mean(death_flag),by=match_id]$V1)

## error about scaling
# probit.re <-
#   lme4::glmer(death_flag ~ -1 + total_direct_cost + group + (1|match_id), family=binomial(link="probit"), data=dtep.money,
#               nAGQ = 100)
# summary(probit.re)
#Error in pwrssUpdate(pp, resp, tol = tolPwrss, GQmat = GQmat, compDev = compDev,  :
# PIRLS loop resulted in NaN value
# In addition: Warning message:
#   Some predictor variables are on very different scales: consider rescaling

# probit.re <-
#   lme4::glmer(death_flag ~ I(total_direct_cost/1e6) + (1|match_id), family=binomial(link="probit"), data=dtep.money,
#               nAGQ = 100)
# summary(probit.re)
#
# probit.re <-
#   lme4::glmer(death_flag ~ scale(total_direct_cost) + group + (1|match_id), family=binomial(link="probit"), data=dtep.money,
#               nAGQ = 100)
# summary(probit.re)
#
#
# probit.re <-
#   lme4::glmer(death_flag ~ -1 + scale(total_direct_cost) * group + (1|match_id), family=binomial(link="probit"), data=dtep.money,
#               nAGQ = 100)
# summary(probit.re)

##remove match_ids with missing
bad.ids <- dtep.money[is.na(total_direct_cost)]$match_id
dtep.money2 <- copy(dtep.money[!(match_id %in% bad.ids)])
nrow(dtep.money)
nrow(dtep.money2)
rm(dtep.money)
dtep.money <- copy(dtep.money2)
rm(dtep.money2)
setkey(dtep.money, match_id, group)
dtep.money$match_id <- rep(1:(nrow(dtep.money)/3),each=3)  ## they have to be in order.

probit.re <-
  lme4::glmer(death_flag ~ -1 +
                I((group=="Colistin case")+0L) +
                I((group=="Non-Colistin case")+0L)+
                I((group=="Colistin case") * scale(total_direct_cost) ) +
                I((group=="Non-Colistin case") * scale(total_direct_cost) ) +
                (1|match_id),
              family=binomial(link="probit"),
              data=dtep.money,
              nAGQ = 100)
summary(probit.re)

mean(scale(dtep.money$total_direct_cost),na.rm=TRUE)
var(scale(dtep.money$total_direct_cost),na.rm=TRUE)

range(scale(dtep.money$total_direct_cost),na.rm=TRUE)

plot(seq(-1,19,0.1), pnorm(-0.59 + 0.156*seq(-1,19,0.1)), ylim=c(0,1))
points(seq(-1,19,0.1), pnorm(-0.91 + 0.142*seq(-1,19,0.1)), col="red")

## 18 times the central case!
## https://fooledbyrandomness.com/FT-MandelbrotTaleb.pdf
range(pnorm(-0.59 + 0.156*seq(-1,19,0.1))- pnorm(-0.91 + 0.142*seq(-1,19,0.1)))


## untransformed scale: error
# probit.re <-
#   lme4::glmer(death_flag ~ -1 +
#                 I((group=="Colistin case")+0L) +
#                 I((group=="Non-Colistin case")+0L)+
#                 I((group=="Colistin case") * (total_direct_cost) ) +
#                 I((group=="Non-Colistin case") * (total_direct_cost) ) +
#                 (1|match_id),
#               family=binomial(link="probit"),
#               data=dtep.money,
#               nAGQ = 100)
# summary(probit.re)
#
# Error in pwrssUpdate(pp, resp, tol = tolPwrss, GQmat = GQmat, compDev = compDev,  :
#                        PIRLS loop resulted in NaN value
#                      In addition: Warning message:
#                        Some predictor variables are on very different scales: consider rescaling


#can gnlrim help? A: had some hiccups that weren't untransformed scale so try untransformed scale
##A: okay htat bombed so try transfomred scale()

b0.cc <-   with(dtep.money, I((group=="Colistin case")+0L))
b0.nc <-   with(dtep.money, I((group=="Non-Colistin case")+0L))

b1.cc <-   as.numeric(with(dtep.money, I((group=="Colistin case")*scale(total_direct_cost))))
b1.nc <-   as.numeric(with(dtep.money, I((group=="Non-Colistin case")*scale(total_direct_cost))))

id <- rep(1:(nrow(dtep.money)/3),each=3) ##dtep.money$match_id
o.value <- with(dtep.money, cbind(death_flag,1-death_flag))

## if you have NA's in the data it may complain about initial values.
## if your ids are not sequential you'll get a likelihood is NA INF or probs too small error
(rand.int.probit <-
    gnlrim::gnlrim(y=o.value,
                   mu = ~ pnorm(beta0.nc*b0.nc + beta0.cc*b0.cc + slope1.nc*b1.nc + slope1.cc*b1.cc + rand1),
                   pmu = c(beta0.nc=0, beta0.cc=0, slope1.nc=0, slope1.cc=0),
                   pmix=c(var1=1),
                   p_uppb = c(  1,   1,  1, 1, 4),
                   p_lowb = c( -1,  -1, -1,-1, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
## errors out with conv code 9999 after 0: iteration if you dont use scale

## if you have NA's in the data it may complain about initial values.
## if your ids are not sequential you'll get a likelihood is NA INF or probs too small error
(rand.int.stable <-
    gnlrim::gnlrim(y=o.value,
                   mu = ~ stable_cdf2(beta0.nc*b0.nc + beta0.cc*b0.cc + beta1.nc*b1.nc + beta1.cc*b1.cc  + rand1, c(alp,0,1,0)),
                   pmu = c(beta0.nc=-1.31, beta0.cc=-0.83,beta1.nc=0.20, beta1.cc=0.20, alp=2),
                   pmix=c(alp=2,scl=0.25),
                   p_uppb = c(  1,1,1,1,    2, 3),
                   p_lowb = c( -2,-2,-1,-1, 0.800, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-scl",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
#errors with 9999 code. if you dont use scale
## use scale but the bound on alpha is wrong, causing you to re-run:
# beta0.nc beta0.cc beta1.nc beta1.cc      alp      scl
# -1.0      0.0      0.0      0.0      1.5      1.0
# [1] 5214.961
# fn is  fn
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -2   <? -1   <? 1     In Bounds
# par[ 2 ]:  -1   <? 0   <? 1     In Bounds
# par[ 3 ]:  -1   <? 0   <? 1     In Bounds
# par[ 4 ]:  -1   <? 0   <? 1     In Bounds
# par[ 5 ]:  0.1   <? 1.5   <? 1.9     In Bounds
# par[ 6 ]:  0.05   <? 1   <? 3     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.1760913   log bounds ratio= 0.2218487
# Method:  nlminb
# 0:     5214.9612: -1.00000  0.00000  0.00000  0.00000  1.50000  1.00000
# 1:     4736.6571: -1.48754 -0.604270 0.173822 0.187670  1.61612 0.435870
# 2:     4698.3961: -1.45533 -0.730714 0.194348 0.200567  1.65183 0.317871
# 3:     4684.4352: -1.36766 -0.854100 0.215510 0.228318  1.67765 0.228327
# 4:     4683.8708: -1.33907 -0.844905 0.206427 0.243941  1.68579 0.223985
# 5:     4683.7849: -1.34276 -0.841232 0.208917 0.238532  1.71119 0.199382
# 6:     4683.5561: -1.31762 -0.839319 0.200911 0.227781  1.72802 0.214046
# 7:     4683.3524: -1.32819 -0.835579 0.209620 0.248460  1.75226 0.223881
# 8:     4683.1410: -1.32729 -0.845423 0.196354 0.231659  1.77856 0.231979
# 9:     4683.0315: -1.31361 -0.839682 0.209551 0.235920  1.80343 0.248778
# 10:     4682.7915: -1.32730 -0.836349 0.203799 0.236994  1.83588 0.243509
# 11:     4682.6147: -1.31335 -0.843177 0.192263 0.230491  1.85977 0.261554
# 12:     4682.4350: -1.31542 -0.842499 0.208571 0.220794  1.89057 0.261364
# 13:     4682.3314: -1.31329 -0.837779 0.198387 0.226003  1.89518 0.264884
# 14:     4682.3041: -1.31345 -0.836506 0.199597 0.225166  1.89984 0.265632
# 15:     4682.3020: -1.31297 -0.836873 0.199407 0.224655  1.90000 0.266507
# 16:     4682.3015: -1.31341 -0.837402 0.199622 0.224267  1.90000 0.267526
# 17:     4682.3014: -1.31356 -0.837414 0.199339 0.224495  1.90000 0.267622
# 18:     4682.3014: -1.31360 -0.837451 0.199359 0.224474  1.90000 0.267779








####################################################
## END Sameer XDR COP idea                ##########
####################################################








####################################################
## the ol' cop == 1.89, 1.2 stable        ##########
####################################################
## is there better plogis approx than 1.89 1.2 ####
dens_diff <- function(THETA){

  a <- THETA[1]
  c <- THETA[2]

  s <- seq(0.0001,0.9999,0.0001)
  s <- seq(0.001 ,0.999 ,0.001 )
  s <- seq(0.01  ,0.99  ,0.01  )
  sqrt(mean(abs(qlogis(s) - libstable4u::stable_q(s, c(a,0,c,0),1))))

}

## depends on s in the function
optim(c(1.89,1.2),
      dens_diff,
      lower=c(1,0.1),
      upper=c(2,10),
      method="L-BFGS-B"
)


a <- 1.80
x <- seq(-3,3,1)
plot(x,dlogis(x))
points(x,libstable4u::stable_pdf(x, c(a,0,1.158,0),1), col="blue")

s <- seq(0.0001,0.9999,0.0001)
plot(qlogis(s),
     libstable4u::stable_q(s, c(a,0,1.2,0),1)
)
abline(a=0,b=1,col="blue")


####################################################
## the ol' cop == 1.89, 1.2 stable        ##########
####################################################

library(data.table)
library(mvtnorm)

## this file contains a subset of the models ran for the COP paper with coauthors
## Caffo and Crainiceanu.

## data from Kung-Yee Liang JRSS-B paper
## n.(left eye).(right eye).(race) is the number
## of people that had eye disease configuration
## and race b (black) or w (white)
n.1.1.b <- sum(10, 14, 28, 56)
n.1.1.w <- sum( 4,  9, 11, 79)
n.1.0.b <- sum(19, 24, 22, 29)
n.1.0.w <- sum(11, 15, 31, 60)
n.0.1.b <- sum(21, 23, 21, 37)
n.0.1.w <- sum(15, 16, 37, 67)
n.0.0.b <- sum(729, 551, 452, 307)
n.0.0.w <- sum(602, 541, 752, 606)

## Example 3.11
## note: can connect to past example of Restricted
T <- function(x, phi) qnorm( pnorm(x), sd = 1/phi )
## T <- function(x, phi) qnorm( pnorm(x), sd = 1 )

## FJ being evaluated over partition subsets with borders T
p.1.1 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(-Inf,-Inf), c(T(Beta0 + x*Beta1, phi), T(Beta0 + x*Beta1, phi)),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.1.0 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(-Inf,T(Beta0 + x*Beta1, phi)), c(T(Beta0 + x*Beta1, phi), Inf),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.0.1 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(T(Beta0 + x*Beta1, phi), -Inf), c(Inf, T(Beta0 + x*Beta1, phi)),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.0.0 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(T(Beta0 + x*Beta1, phi), T(Beta0 + x*Beta1, phi)), c(Inf, Inf),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.1.1.b <- function(Beta0, Beta1, rho, phi) p.1.1(x=1, Beta0, Beta1, rho, phi)
p.1.1.w <- function(Beta0, Beta1, rho, phi) p.1.1(x=0, Beta0, Beta1, rho, phi)
p.1.0.b <- function(Beta0, Beta1, rho, phi) p.1.0(x=1, Beta0, Beta1, rho, phi)
p.1.0.w <- function(Beta0, Beta1, rho, phi) p.1.0(x=0, Beta0, Beta1, rho, phi)
p.0.1.b <- function(Beta0, Beta1, rho, phi) p.0.1(x=1, Beta0, Beta1, rho, phi)
p.0.1.w <- function(Beta0, Beta1, rho, phi) p.0.1(x=0, Beta0, Beta1, rho, phi)
p.0.0.b <- function(Beta0, Beta1, rho, phi) p.0.0(x=1, Beta0, Beta1, rho, phi)
p.0.0.w <- function(Beta0, Beta1, rho, phi) p.0.0(x=0, Beta0, Beta1, rho, phi)

## function of log-likelihood depending on 'theta'
logl.theta <- function(theta){
  Beta0 <- theta[1]
  Beta1 <- theta[2]
  rho   <- theta[3] / (1+theta[3])
  phi   <- 1 / (sqrt(1+theta[3]))
  ##rho <- .4363368
  -(n.1.1.b*log(p.1.1.b(Beta0, Beta1, rho, phi)) +
      n.1.1.w*log(p.1.1.w(Beta0, Beta1, rho, phi)) +
      n.1.0.b*log(p.1.0.b(Beta0, Beta1, rho, phi)) +
      n.1.0.w*log(p.1.0.w(Beta0, Beta1, rho, phi)) +
      n.0.1.b*log(p.0.1.b(Beta0, Beta1, rho, phi)) +
      n.0.1.w*log(p.0.1.w(Beta0, Beta1, rho, phi)) +
      n.0.0.b*log(p.0.0.b(Beta0, Beta1, rho, phi)) +
      n.0.0.w*log(p.0.0.w(Beta0, Beta1, rho, phi)))
}

## maximize the likelihood, theta[3] is tau^2 (see logl.theta above)
probitprobitnormal <- optim(c( 0,0,1),
                            logl.theta,
                            method="L-BFGS-B",
                            control= list(maxit=100),
                            upper=c( 3, 1, 10),
                            lower=c(-3,-1,-10))
# > probitprobitnormal
# $par
# [1] -1.39879349  0.03942097  2.87840414
#
# $value
# [1] 2699.798
#
# $counts
# function gradient
# 21       21
#
# $convergence
# [1] 0
#
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

probitprobitnormal$val
probitprobitnormal$par
## rho and phi
(tau2 <- probitprobitnormal$par[3])
(rho   <- tau2 / (1+tau2))
(phi   <- 1 / (sqrt(1+tau2)))
(margBeta1 <- probitprobitnormal$par[2])
(condBeta1 <- phi^(-1)*probitprobitnormal$par[2])
(margBeta0 <- probitprobitnormal$par[1])
(condBeta0 <- phi^(-1)*probitprobitnormal$par[1])


data(bwVI)
lme4::glmer(value ~ black + (1|id), data=bwVI, family=binomial("probit"), nAGQ = 100)
1.697^2
setDT(bwVI)

attach(bwVI)
o.value <- value ## lesson learned!  got an error when I used "value"
(rand.int.stable <-
    gnlrim::gnlrim(y=cbind(o.value, 1-o.value), ## lesson learned!  got an error when I used "value"
                   mu = ~ stable_cdf2(Intercept + black*b_p + rand1, c(alpha, 0, 1.2, 0)),
                   ##pmu = c(Intercept = -2.7, b_p=0.08, alpha=1.89),
                   ## -4.75073 0.129755  1.89000  1.97264
                   pmu = c(Intercept = -4.88, b_p=0.134, alpha=1.89),
                   pmix=c(alpha=1.89,scl=2.03),
                   p_uppb = c(  50,   9, 1.89, 20.00),
                   p_lowb = c( -50,  -9, 1.89,  0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-scl",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
#
# Intercept       b_p     alpha       scl
# -2.70      0.08      1.89      2.88
# [1] 3461.13
# fn is  fn
# Looking for method =  nlminb
# Function has  4  arguments
# par[ 1 ]:  -50   <? -2.7   <? 50     In Bounds
# par[ 2 ]:  -9   <? 0.08   <? 9     In Bounds
# par[ 3 ]:  1.89   <? 1.89   <? 1.89     In Bounds
# par[ 4 ]:  0.05   <? 2.88   <? 20     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 1.556303   log bounds ratio= 0.7447275
# Method:  nlminb
# 0:     3461.1303: -2.70000 0.0800000  1.89000  2.88000
# 1:     2934.5180: -3.47493 -0.267695  1.89000  2.35219
# 2:     2715.0814: -3.98848 -0.466517  1.89000  1.78837
# 3:     2707.0199: -4.37736 0.201170  1.89000  1.63311
# 4:     2705.6610: -4.31865 0.218105  1.89000  1.79977
# 5:     2701.2813: -4.49575 0.209398  1.89000  1.80807
# 6:     2700.4145: -4.60898 0.0878555  1.89000  1.87063
# 7:     2700.1183: -4.59699 0.0963923  1.89000  1.90390
# 8:     2699.9953: -4.63191 0.104887  1.89000  1.90952
# 9:     2699.9209: -4.68627 0.127812  1.89000  1.95209
# 10:     2699.8230: -4.75394 0.148075  1.89000  1.96950
# 11:     2699.8154: -4.75496 0.130339  1.89000  1.97386
# 12:     2699.8150: -4.75091 0.129788  1.89000  1.97273
# 13:     2699.8150: -4.75071 0.129752  1.89000  1.97263
# 14:     2699.8150: -4.75073 0.129755  1.89000  1.97264
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p      alpha        scl
# -4.7507258  0.1297551  1.8900000  1.9726394
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2699.815
#
# $fevals
# function
# 23
#
# $gevals
# gradient
# 55
#
# $nitns
# [1] 14
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 15429.81
#
# Assemble the answers
#        Intercept       b_p alpha      scl    value fevals gevals niter convcode kkt1 kkt2
# nlminb -4.750726 0.1297551  1.89 1.972639 2699.815     23     55    14        0   NA   NA
# xtime
# nlminb 15429.81

lme4::glmer(value ~ black + (1|id), data=bwVI, family=binomial("logit"), nAGQ = 100)
# > lme4::glmer(value ~ black + (1|id), data=bwVI, family=binomial("logit"), nAGQ = 100)
# Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite
#                                                           Quadrature, nAGQ = 100) [glmerMod]
# Family: binomial  ( logit )
# Formula: value ~ black + (1 | id)
# Data: bwVI
# AIC       BIC    logLik  deviance  df.resid
# 5405.531  5427.279 -2699.765  5399.531     10395
# Random effects:
#   Groups Name        Std.Dev.
# id     (Intercept) 3.043
# Number of obs: 10398, groups:  id, 5199
# Fixed Effects:
#   (Intercept)        black
# -4.9409       0.1462
####################################################
## END                  ######            ##########
####################################################

####################################################
## the ol' cop == gnlmm trick for stable? ##########
####################################################



## Bruce Swihart (bswihart@jhsph.edu)
library(data.table)
library(mvtnorm)

## this file contains a subset of the models ran for the COP paper with coauthors
## Caffo and Crainiceanu.

## data from Kung-Yee Liang JRSS-B paper
## n.(left eye).(right eye).(race) is the number
## of people that had eye disease configuration
## and race b (black) or w (white)
n.1.1.b <- sum(10, 14, 28, 56)
n.1.1.w <- sum( 4,  9, 11, 79)
n.1.0.b <- sum(19, 24, 22, 29)
n.1.0.w <- sum(11, 15, 31, 60)
n.0.1.b <- sum(21, 23, 21, 37)
n.0.1.w <- sum(15, 16, 37, 67)
n.0.0.b <- sum(729, 551, 452, 307)
n.0.0.w <- sum(602, 541, 752, 606)

## Example 3.11
## note: can connect to past example of Restricted
T <- function(x, phi) qnorm( pnorm(x), sd = 1/phi )
## T <- function(x, phi) qnorm( pnorm(x), sd = 1 )

## FJ being evaluated over partition subsets with borders T
p.1.1 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(-Inf,-Inf), c(T(Beta0 + x*Beta1, phi), T(Beta0 + x*Beta1, phi)),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.1.0 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(-Inf,T(Beta0 + x*Beta1, phi)), c(T(Beta0 + x*Beta1, phi), Inf),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.0.1 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(T(Beta0 + x*Beta1, phi), -Inf), c(Inf, T(Beta0 + x*Beta1, phi)),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.0.0 <-function(x, Beta0, Beta1, rho, phi){
  pmvnorm(c(T(Beta0 + x*Beta1, phi), T(Beta0 + x*Beta1, phi)), c(Inf, Inf),
          sigma=(1/phi)^2*matrix(c(1,rho,rho,1), nrow=2))[1]
}
p.1.1.b <- function(Beta0, Beta1, rho, phi) p.1.1(x=1, Beta0, Beta1, rho, phi)
p.1.1.w <- function(Beta0, Beta1, rho, phi) p.1.1(x=0, Beta0, Beta1, rho, phi)
p.1.0.b <- function(Beta0, Beta1, rho, phi) p.1.0(x=1, Beta0, Beta1, rho, phi)
p.1.0.w <- function(Beta0, Beta1, rho, phi) p.1.0(x=0, Beta0, Beta1, rho, phi)
p.0.1.b <- function(Beta0, Beta1, rho, phi) p.0.1(x=1, Beta0, Beta1, rho, phi)
p.0.1.w <- function(Beta0, Beta1, rho, phi) p.0.1(x=0, Beta0, Beta1, rho, phi)
p.0.0.b <- function(Beta0, Beta1, rho, phi) p.0.0(x=1, Beta0, Beta1, rho, phi)
p.0.0.w <- function(Beta0, Beta1, rho, phi) p.0.0(x=0, Beta0, Beta1, rho, phi)

## function of log-likelihood depending on 'theta'
logl.theta <- function(theta){
  Beta0 <- theta[1]
  Beta1 <- theta[2]
  rho   <- theta[3] / (1+theta[3])
  phi   <- 1 / (sqrt(1+theta[3]))
  ##rho <- .4363368
  -(n.1.1.b*log(p.1.1.b(Beta0, Beta1, rho, phi)) +
      n.1.1.w*log(p.1.1.w(Beta0, Beta1, rho, phi)) +
      n.1.0.b*log(p.1.0.b(Beta0, Beta1, rho, phi)) +
      n.1.0.w*log(p.1.0.w(Beta0, Beta1, rho, phi)) +
      n.0.1.b*log(p.0.1.b(Beta0, Beta1, rho, phi)) +
      n.0.1.w*log(p.0.1.w(Beta0, Beta1, rho, phi)) +
      n.0.0.b*log(p.0.0.b(Beta0, Beta1, rho, phi)) +
      n.0.0.w*log(p.0.0.w(Beta0, Beta1, rho, phi)))
}

## maximize the likelihood, theta[3] is tau^2 (see logl.theta above)
probitprobitnormal <- optim(c( 0,0,1),
                            logl.theta,
                            method="L-BFGS-B",
                            control= list(maxit=100),
                            upper=c( 3, 1, 10),
                            lower=c(-3,-1,-10))
# > probitprobitnormal
# $par
# [1] -1.39879349  0.03942097  2.87840414
#
# $value
# [1] 2699.798
#
# $counts
# function gradient
# 21       21
#
# $convergence
# [1] 0
#
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

probitprobitnormal$val
probitprobitnormal$par
## rho and phi
(tau2 <- probitprobitnormal$par[3])
(rho   <- tau2 / (1+tau2))
(phi   <- 1 / (sqrt(1+tau2)))
(margBeta1 <- probitprobitnormal$par[2])
(condBeta1 <- phi^(-1)*probitprobitnormal$par[2])
(margBeta0 <- probitprobitnormal$par[1])
(condBeta0 <- phi^(-1)*probitprobitnormal$par[1])


data(bwVI)
lme4::glmer(value ~ black + (1|id), data=bwVI, family=binomial("probit"), nAGQ = 100)
1.697^2
setDT(bwVI)
#summed_binom <- bwVI[,{j=list(eyes_impaired=sum(value),eyes_tested=.N)}, by=c("id","black")]
#attach(summed_binom)
## PPN -- initially I didn't get the "expected likelihood" so I tested to make sure
##        gnlrim and gnlrem gave same answers.  See below where I don't sum things...
# (rand.int <-
#     gnlrim::gnlrem(y=cbind(eyes_impaired, eyes_tested - eyes_impaired),
#                    mu = ~ pnorm(Intercept + black*b_p + rand1),
#                    pmu = c(Intercept=-0.95, b_p=0.55),
#                    pmix=c(var1=1),
#                    p_uppb = c(  0,   2, 4.00),
#                    p_lowb = c( -4,  -2, 0.05),
#                    distribution="binomial",
#                    nest=id,
#                    random=c("rand1"),
#                    mixture="normal-var",
#                    ooo=TRUE,
#                    compute_hessian = FALSE,
#                    compute_kkt = FALSE,
#                    trace=1,
#                    method='nlminb',
#     )
# )
# Assemble the answers
# Intercept        b_p     var1    value fevals gevals niter convcode kkt1 kkt2
# nlminb -2.754239 0.07762926 2.876914 2389.273     31     72    20        0   NA   NA
# xtime
# nlminb 1658.171
# (rand.int.gnlrim <-
#     gnlrim::gnlrim(y=cbind(eyes_impaired, eyes_tested - eyes_impaired),
#                    mu = ~ pnorm(Intercept + black*b_p + rand1),
#                    pmu = c(Intercept=-0.95, b_p=0.55),
#                    pmix=c(var1=1),
#                    p_uppb = c(  0,   2, 4.00),
#                    p_lowb = c( -4,  -2, 0.05),
#                    distribution="binomial",
#                    nest=id,
#                    random=c("rand1"),
#                    mixture="normal-var",
#                    ooo=TRUE,
#                    compute_hessian = FALSE,
#                    compute_kkt = FALSE,
#                    trace=1,
#                    method='nlminb',
#     )
# )
# Assemble the answers
# Intercept        b_p     var1    value fevals gevals niter convcode kkt1 kkt2
# nlminb -2.754237 0.07762872 2.876913 2389.273     21     64    17        0   NA   NA
# xtime
# nlminb 1333.852

#detach(summed_binom)
attach(bwVI)
o.value <- value ## lesson learned!  got an error when I used "value"
(rand.int.gnlrim <-
    gnlrim::gnlrim(y=cbind(o.value, 1-o.value),
                   mu = ~ pnorm(Intercept + black*b_p + rand1),
                   pmu = c(Intercept=-0.95, b_p=0.55),
                   pmix=c(var1=1),
                   p_uppb = c(  0,   2, 4.00),
                   p_lowb = c( -4,  -2, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
# Assemble the answers
# Intercept        b_p     var1    value fevals gevals niter convcode kkt1 kkt2
# nlminb -2.754237 0.07762872 2.876913 2699.803     21     64    17        0   NA   NA

(rand.int.stable <-
    gnlrim::gnlrim(y=cbind(o.value, 1-o.value), ## lesson learned!  got an error when I used "value"
                   mu = ~ stable_cdf2(Intercept + black*b_p + rand1, c(alpha, 0, 1, 0)),
                   pmu = c(Intercept = -2.7, b_p=0.08, alpha=1.8),
                   pmix=c(alpha=1.8,scl=2.88),
                   p_uppb = c(  50,   9, 1.99, 20.00),
                   p_lowb = c( -50,  -9, 1.05,  0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-scl",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)

# Method:  nlminb
# 0:     2745.8764: -5.00000 -3.24000  1.25000  1.52000
# 1:     2737.7009: -5.20990 -2.56549  1.01000  1.14402
# 2:     2725.4479: -5.34348 -2.24165  1.13678 0.904028
# 3:     2713.3670: -5.52718 -1.84601  1.06002 0.918473
# 4:     2711.8555: -5.54324 -1.88506  1.11478  1.09971
# 5:     2711.5921: -5.61401 -1.71030  1.06926  1.10168
# 6:     2708.1815: -5.68973 -1.53588  1.10251  1.08248
# 7:     2703.3790: -5.96875 -0.813365  1.08192  1.04080
# 8:     2700.0922: -5.75068 -0.0843337  1.16032  1.17062
# 9:     2699.5341: -5.95969 0.0846123  1.16950  1.17735
# 10:     2699.5101: -5.82060 0.313741  1.19095  1.18469
# 11:     2699.4617: -5.97930 0.183156  1.17878  1.23352
# 12:     2699.4160: -5.94655 0.181168  1.17476  1.19521
# 13:     2699.4112: -5.94868 0.207219  1.17623  1.20046
# 14:     2699.4103: -5.96215 0.207706  1.17457  1.19995
# 15:     2699.4049: -6.11683 0.219983  1.15559  1.19040
# 16:     2699.4042: -6.12561 0.222487  1.15524  1.18992
# 17:     2699.4040: -6.15028 0.225216  1.15259  1.18817
# 18:     2699.4040: -6.14895 0.225096  1.15276  1.18826
# 19:     2699.4040: -6.14913 0.224992  1.15274  1.18826
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p      alpha        scl
# -6.1491302  0.2249925  1.1527411  1.1882609
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2699.404
#
# $fevals
# function
# 26
#
# $gevals
# gradient
# 102
#
# $nitns
# [1] 19
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 12956.43
#
# Assemble the answers
#        Intercept       b_p    alpha      scl    value fevals gevals niter convcode kkt1
# nlminb  -6.14913 0.2249925 1.152741 1.188261 2699.404     26    102    19        0   NA

## tried a re-run:
# > o.value <- value ## lesson learned!  got an error when I used "value"
# > (rand.int.stable <-
#      +     gnlrim::gnlrim(y=cbind(o.value, 1-o.value), ## lesson learned!  got an error when I used "value"
#                           +                    mu = ~ stable_cdf2(Intercept + black*b_p + rand1, c(alpha, 0, 1, 0)),
#                           +                    pmu = c(Intercept = -2.7, b_p=0.08, alpha=1.8),
#                           +                    pmix=c(alpha=1.8,scl=2.88),
#                           +                    p_uppb = c(  50,   9, 1.99, 20.00),
#                           +                    p_lowb = c( -50,  -9, 1.05,  0.05),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1"),
#                           +                    mixture="libstable4u-subgauss-scl",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +     )
#    + )
# [1] 4
# Intercept       b_p     alpha       scl
# -2.70      0.08      1.80      2.88
# [1] 3430.29
# fn is  fn
# Looking for method =  nlminb
# Function has  4  arguments
# par[ 1 ]:  -50   <? -2.7   <? 50     In Bounds
# par[ 2 ]:  -9   <? 0.08   <? 9     In Bounds
# par[ 3 ]:  1.05   <? 1.8   <? 1.99     In Bounds
# par[ 4 ]:  0.05   <? 2.88   <? 20     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 1.556303   log bounds ratio= 2.026872
# Method:  nlminb
# 0:     3430.2898: -2.70000 0.0800000  1.80000  2.88000
# 1:     2884.1373: -3.42798 -0.246897  1.99000  2.32575
# 2:     2727.5734: -3.71540 -0.355767  1.99000  1.94951
# 3:     2708.5006: -3.83659 -0.373315  1.99000  1.74881
# 4:     2703.0968: -3.89486 -0.161453  1.90695  1.73893
# 5:     2702.4809: -3.91506 -0.150820  1.93647  1.69825
# 6:     2701.4776: -3.89743 -0.102593  1.93890  1.71837
# 7:     2700.9524: -3.89844 -0.0231217  1.99000  1.69087
# 8:     2700.2269: -3.86544 0.170356  1.99000  1.67962
# 9:     2700.1937: -4.04618 0.179886  1.99000  1.75697
# 10:     2700.1016: -4.02595 0.154927  1.99000  1.73926
# 11:     2699.9558: -4.00552 0.138957  1.96788  1.72570
# 12:     2699.8339: -3.93927 0.137515  1.99000  1.70717
# 13:     2699.8253: -3.88489 0.0811044  1.99000  1.69049
# 14:     2699.7956: -3.90269 0.110038  1.99000  1.69314
# 15:     2699.7956: -3.90272 0.109646  1.99000  1.69323
# 16:     2699.7956: -3.90269 0.109644  1.99000  1.69321
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p      alpha        scl
# -3.9026857  0.1096438  1.9900000  1.6932125
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2699.796
#
# $fevals
# function
# 23
#
# $gevals
# gradient
# 78
#
# $nitns
# [1] 16
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 21196.23
#
# Assemble the answers
# Intercept       b_p alpha      scl    value fevals gevals niter convcode kkt1 kkt2
# nlminb -3.902686 0.1096438  1.99 1.693212 2699.796     23     78    16        0   NA   NA

library(stable, lib.loc="~/RLIBS")
#
# T <- function(eta_m, g_g=gam_g, g_d=gam_d, g_h=gam_h, aa=a){
#   phi <- g_d / (g_h^aa + g_g^aa)^(1/aa)
#   eta_m / phi
# }
#
# ## Here we do for Ji=2.  We have graphical evidence suggesting a spectral measure at the usual basis (axis)
# ## as well as (1/sqrt(2),1/sqrt(2)) and (-1/sqrt(2),-1/sqrt(2))
#
# ## assign the lambda weights...
# lambda.vec <- c(rep(gam_h^a,Ji),
#                 Ji^(a/2)*gam_g^a  )
# ## assign the positions of the weights on the unit sphere
# s.vec <- matrix(cbind(diag(Ji),
#                       rep(1,Ji)/sqrt(Ji)),
#                 nrow=Ji)
# ## define the corresponding discrete spectral measure stable dist.
#
#
# dist.specJ <- mvstable.discrete.spec.meas  (alpha=a,
#                                             s=s.vec,
#                                             lambda=lambda.vec,
#                                             beta=rep(0,length(lambda.vec)) )#read the manual - length(beta) == length(lambda)
# mvstable.info(dist.specJ)


## data from Kung-Yee Liang JRSS-B paper
## n.(left eye).(right eye).(race) is the number
## of people that had eye disease configuration
## and race b (black) or w (white)
n.1.1.b <- sum(10, 14, 28, 56)
n.1.1.w <- sum( 4,  9, 11, 79)
n.1.0.b <- sum(19, 24, 22, 29)
n.1.0.w <- sum(11, 15, 31, 60)
n.0.1.b <- sum(21, 23, 21, 37)
n.0.1.w <- sum(15, 16, 37, 67)
n.0.0.b <- sum(729, 551, 452, 307)
n.0.0.w <- sum(602, 541, 752, 606)

## Example 3.11
## note: can connect to past example of Restricted
## T <- function(x, phi) qnorm( pnorm(x), sd = 1/phi )
## T <- function(x, phi) qnorm( pnorm(x), sd = 1 )
T <- function(x, a, gam_g) (1+gam_g^a)^(1/a) * x

## FJ being evaluated over partition subsets with borders T
p.1.1 <-function(x, Beta0, Beta1, a, gam_g){
  big.val <- 100000
  gam_h <- 1
  Ji <- 2

  ## assign the lambda weights...
  lambda.vec <- c(rep(gam_h^a,Ji),
                  Ji^(a/2)*gam_g^a  )
  ## assign the positions of the weights on the unit sphere
  s.vec <- matrix(cbind(diag(Ji),
                        rep(1,Ji)/sqrt(Ji)),
                  nrow=Ji)
  ## define the corresponding discrete spectral measure stable dist.


  dist.specJ <- mvstable.discrete.spec.meas  (alpha=a,
                                              s=s.vec,
                                              lambda=lambda.vec,
                                              beta=rep(0,length(lambda.vec)) )#read the manual - length(beta) == length(lambda)
  mvstable.info(dist.specJ)

  ## pmvstable(dist.specJ, c(-big.val,-big.val), c(T(Beta0 + x*Beta1, a, gam_g), T(Beta0 + x*Beta1, a, gam_g)), epsabs=1e-4)
  pmvstable.MC(dist.specJ, c(-big.val,-big.val), c(T(Beta0 + x*Beta1, a, gam_g), T(Beta0 + x*Beta1, a, gam_g)), n=4e6)

}
#
# T(x=-6.1, a=1.15, gam_g=1.18)
# p.1.1(x=0,-6.14,0.225, 1.15, 1.18)
# p.1.1(x=1,-6.14,0.225, 1.15, 1.18)
#
p.1.0 <-function(x, Beta0, Beta1, a, gam_g){
  big.val <- 100000
  gam_h <- 1
  Ji <- 2

  ## assign the lambda weights...
  lambda.vec <- c(rep(gam_h^a,Ji),
                  Ji^(a/2)*gam_g^a  )
  ## assign the positions of the weights on the unit sphere
  s.vec <- matrix(cbind(diag(Ji),
                        rep(1,Ji)/sqrt(Ji)),
                  nrow=Ji)
  ## define the corresponding discrete spectral measure stable dist.


  dist.specJ <- mvstable.discrete.spec.meas  (alpha=a,
                                              s=s.vec,
                                              lambda=lambda.vec,
                                              beta=rep(0,length(lambda.vec)) )#read the manual - length(beta) == length(lambda)
  mvstable.info(dist.specJ)

  ## pmvstable(dist.specJ, c(-big.val,-big.val), c(T(Beta0 + x*Beta1, a, gam_g), T(Beta0 + x*Beta1, a, gam_g)), epsabs=1e-4)

  pmvstable.MC(dist.specJ,c( -big.val,T(Beta0 + x*Beta1, a, gam_g)), c(T(Beta0 + x*Beta1, a, gam_g),big.val),  n=4e6)
}
p.0.1 <-function(x, Beta0, Beta1, a, gam_g){
  big.val <-100000
  gam_h <- 1
  Ji <- 2

  ## assign the lambda weights...
  lambda.vec <- c(rep(gam_h^a,Ji),
                  Ji^(a/2)*gam_g^a  )
  ## assign the positions of the weights on the unit sphere
  s.vec <- matrix(cbind(diag(Ji),
                        rep(1,Ji)/sqrt(Ji)),
                  nrow=Ji)
  ## define the corresponding discrete spectral measure stable dist.


  dist.specJ <- mvstable.discrete.spec.meas  (alpha=a,
                                              s=s.vec,
                                              lambda=lambda.vec,
                                              beta=rep(0,length(lambda.vec)) )#read the manual - length(beta) == length(lambda)
  mvstable.info(dist.specJ)

  ## pmvstable(dist.specJ, c(-big.val,-big.val), c(T(Beta0 + x*Beta1, a, gam_g), T(Beta0 + x*Beta1, a, gam_g)), epsabs=1e-4)

  pmvstable.MC(dist.specJ,c(T(Beta0 + x*Beta1, a, gam_g), -big.val), c(big.val,T(Beta0 + x*Beta1, a, gam_g)),  n=4e6)
}
p.0.0 <-function(x, Beta0, Beta1, a, gam_g){
  big.val <- 100000
  gam_h <- 1
  Ji <- 2

  ## assign the lambda weights...
  lambda.vec <- c(rep(gam_h^a,Ji),
                  Ji^(a/2)*gam_g^a  )
  ## assign the positions of the weights on the unit sphere
  s.vec <- matrix(cbind(diag(Ji),
                        rep(1,Ji)/sqrt(Ji)),
                  nrow=Ji)
  ## define the corresponding discrete spectral measure stable dist.


  dist.specJ <- mvstable.discrete.spec.meas  (alpha=a,
                                              s=s.vec,
                                              lambda=lambda.vec,
                                              beta=rep(0,length(lambda.vec)) )#read the manual - length(beta) == length(lambda)
  mvstable.info(dist.specJ)

  ## pmvstable(dist.specJ, c(-big.val,-big.val), c(T(Beta0 + x*Beta1, a, gam_g), T(Beta0 + x*Beta1, a, gam_g)), epsabs=1e-4)

  pmvstable.MC(dist.specJ,c(T(Beta0 + x*Beta1, a, gam_g), T(Beta0 + x*Beta1, a, gam_g)), c(big.val,big.val),  n=4e6)
}
T(x=-6.1, a=1.15, gam_g=1.18)

p.1.1(x=0,-6.14,0.225, 1.15, 1.18)
p.1.1(x=1,-6.14,0.225, 1.15, 1.18)

p.0.0(x=0,-6.14,0.225, 1.15, 1.18)
p.0.0(x=1,-6.14,0.225, 1.15, 1.18)

p.1.0(x=0,-6.14,0.225, 1.15, 1.18)
p.1.0(x=1,-6.14,0.225, 1.15, 1.18)

p.0.1(x=0,-6.14,0.225, 1.15, 1.18)
p.0.1(x=1,-6.14,0.225, 1.15, 1.18)



p.1.1.b <- function(Beta0, Beta1, a, gam_g) p.1.1(x=1, Beta0, Beta1, a, gam_g)
p.1.1.w <- function(Beta0, Beta1, a, gam_g) p.1.1(x=0, Beta0, Beta1, a, gam_g)
p.1.0.b <- function(Beta0, Beta1, a, gam_g) p.1.0(x=1, Beta0, Beta1, a, gam_g)
p.1.0.w <- function(Beta0, Beta1, a, gam_g) p.1.0(x=0, Beta0, Beta1, a, gam_g)
p.0.1.b <- function(Beta0, Beta1, a, gam_g) p.0.1(x=1, Beta0, Beta1, a, gam_g)
p.0.1.w <- function(Beta0, Beta1, a, gam_g) p.0.1(x=0, Beta0, Beta1, a, gam_g)
p.0.0.b <- function(Beta0, Beta1, a, gam_g) p.0.0(x=1, Beta0, Beta1, a, gam_g)
p.0.0.w <- function(Beta0, Beta1, a, gam_g) p.0.0(x=0, Beta0, Beta1, a, gam_g)

## function of log-likelihood depending on 'theta'
logl.theta <- function(theta){
  Beta0 <- theta[1]
  Beta1 <- theta[2]
    a   <- theta[3]
  gam_g <- theta[4]
  ##rho <- .4363368
  # -(n.1.1.b*log(p.1.1.b(Beta0, Beta1, a, gam_g)) +
  #     n.1.1.w*log(p.1.1.w(Beta0, Beta1, a, gam_g)) +
  #     n.1.0.b*log(p.1.0.b(Beta0, Beta1, a, gam_g)) +
  #     n.1.0.w*log(p.1.0.w(Beta0, Beta1, a, gam_g)) +
  #     n.0.1.b*log(p.0.1.b(Beta0, Beta1, a, gam_g)) +
  #     n.0.1.w*log(p.0.1.w(Beta0, Beta1, a, gam_g)) +
  #     n.0.0.b*log(p.0.0.b(Beta0, Beta1, a, gam_g)) +
  #     n.0.0.w*log(p.0.0.w(Beta0, Beta1, a, gam_g)))



  ## smarten this up a bit; it will quicken and be more accurate
  ## if you realize 1.0 nad 0.1 are the same and should make 0.0+1.1
  ## equal 1.

  Kp.1.1.b <- p.1.1.b(Beta0, Beta1, a, gam_g)
  Kp.1.1.w <- p.1.1.w(Beta0, Beta1, a, gam_g)

  Kp.0.0.b <- p.0.0.b(Beta0, Beta1, a, gam_g)
  Kp.0.0.w <- p.0.0.w(Beta0, Beta1, a, gam_g)

  remainder.b <- 1 - Kp.1.1.b - Kp.0.0.b
  remainder.w <- 1 - Kp.1.1.w - Kp.0.0.w

  Kp.0.1.b <- Kp.1.0.b <- remainder.b / 2
  Kp.0.1.w <- Kp.1.0.w <- remainder.w / 2


  -(n.1.1.b*log(Kp.1.1.b) +
      n.1.1.w*log(Kp.1.1.w) +
      n.1.0.b*log(Kp.1.0.b) +
      n.1.0.w*log(Kp.1.0.w) +
      n.0.1.b*log(Kp.0.1.b) +
      n.0.1.w*log(Kp.0.1.w) +
      n.0.0.b*log(Kp.0.0.b) +
      n.0.0.w*log(Kp.0.0.w))

}

# Hoping for 2699.404?
#        Intercept       b_p    alpha      scl    value fevals gevals niter convcode kkt1
# nlminb  -6.14913 0.2249925 1.152741 1.188261 2699.404     26    102    19        0   NA
# T <- function(x, a, gam_g) (1+gam_g^a)^(1/a) * x
# logl.theta(c(-6.14913,0.2249925,1.152741,1.188261))
#
# logl.theta(c(-5,0.22,1.15,1.18))


T <- function(x, a, gam_g)  x
logl.theta(c(-6.14913,0.2249925,1.152741,1.188261))
## 2700.102


# Assemble the answers
# Intercept       b_p alpha      scl    value fevals gevals niter convcode kkt1 kkt2
# nlminb -3.902686 0.1096438  1.99 1.693212 2699.796     23     78    16        0   NA   NA
T <- function(x, a, gam_g) (1+gam_g^a)^(1/a) * x
T <- function(x, a, gam_g)  x
logl.theta(c(-3.902686,0.1096438,1.99,1.693212 ))


T <- function(x, a, gam_g) sqrt(2) * (1/sqrt(2)^2+(gam_g)^2)^(1/2) * x
logl.theta(c(-1.398793,0.03942097,1.99999,sqrt(2.87)/sqrt(2))) ## if approaching PPN


## maximize the likelihood
T <- function(x, a, gam_g)  x
sss <- optim(c( -5.1,0.21,1.15,1.1),
                            logl.theta,
                            method="L-BFGS-B",
                            control= list(maxit=100),
                            upper=c(-5, 0.23,  1.2, 3),
                            lower=c(-7, 0.20,  1.1, 1))

sss
sss$val
sss$par
#$par
# [1] -6.4589171  0.2041944  1.1000000  1.1043536
#
# $value
# [1] 2700.217
#
# $counts
# function gradient
# 33       33
#
# $convergence
# [1] 0
#
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"


##https://www.r-bloggers.com/2013/07/optimising-a-noisy-objective-function/#google_vignette

library(DEoptim)
set.seed(1)
R0 = DEoptim(logl.theta,
             lower = c(-4, 0, 1.01, 0.05),
             upper = c(0, 1, 1.99, 5),
             control = DEoptim.control(trace = 10)
             )
# Error in x[, i] : invalid subscript type 'closure'
# Error during wrapup: incorrect number of dimensions
# Error: no more error handlers available (recursive errors?); invoking 'abort' restart



########################################$$$$#########
## 2022-02-17E                                     ##
## START two random parameters: CAUCHY MARGINAL COEFF
#####################################################

## 2022-02-17D was ambition and it bit me.  Let's
## start by doing a cauchy-cauchy-cauchy using
## dmvt(df=1) and then verify it with an alpha=1 locked

library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.0, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvpd::rmvss(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
set.seed(1709)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat
## glmer -- only random intercept
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)




attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## SSS
(rand.int <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1, c(alpha, 0, 1, 0)),
                   pmu = c(Intercept=-1.2, b_p=1, alpha=1),
                   pmix=c(alpha=1, scl=0.25),
                   p_uppb = c(  0,   4, 1, 4.00),
                   p_lowb = c( -4,  -2, 1, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-scl",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
# Assemble the answers
#        Intercept      b_p alpha      scl    value fevals gevals niter convcode kkt1 kkt2  xtime
# nlminb -2.262357 2.035221     1 1.140446 3727.372     15     51    12        0   NA   NA 20.332

## excuse the 1^alpha; had to do it to get alpha to appear before phi in `mu` argument.
(rand.int.marg.phi <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2((Intercept + period_numeric*b_p)*1^alpha/phi + rand1, c(alpha, 0, 1, 0)),
                   pmu = c(Intercept=-0.95, b_p=0.55,  alpha=1.0, phi=0.3),
                   pmix=c(alpha=1.0, phi=0.3),
                   p_uppb = c(  0,   4, 1, 0.99),
                   p_lowb = c( -4,  -2, 1, 0.01),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-phi",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
# Assemble the answers
# Intercept       b_p alpha       phi    value fevals gevals niter convcode kkt1 kkt2  xtime
# nlminb -1.056956 0.9508394     1 0.4671921 3727.372     21     60    16        0   NA   NA 26.287
1/(1+1.140446 ^(1))^(1/1)
#> This is phi from the marginal phi model:  [1]




## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)



## run it with bivariate CAUCHY and lock alpha=1 ... note pmix doesn't have alpha in it
## we'll verify this with bivariate stable in next chunk
## this finally ran successfully, after many dead-ends of programming a t(df=1) distribution
## irony of all ironies -- making a product distribution based function in mvpd was what
## came through.  Need further sims to make sure it is estimating what it should.
## about 44 minutes compared to ## 3hr20min for bivariate-subgauss-corr (alpha lock at 1) approach below
(rand.int.rand.slopes.nonzero.corr.CUBA.cauchy <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),

                   pmu = c(Intercept=-1.2, b_p=1, alpha=1.0),
                   pmix=c(var1=1, var2=1.3, corr12= 0.30),

                   p_uppb = c(  0,   4, 1.0, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-cauchy-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
# > (rand.int.rand.slopes.nonzero.corr.CUBA.cauchy <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, alpha=1.0),
#                           +                    pmix=c(var1=1, var2=1.3, corr12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   4, 1.0, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-cauchy-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.1,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 4     In Bounds
# par[ 3 ]:  1   <? 1   <? 1     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-09-30 08:58:01 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 08:59:50 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.5228787
# Method:  nlminb
# [1] "2023-09-30 08:59:51 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:01:38 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 09:01:38 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:03:25 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69603"
# [1] "2023-09-30 09:03:25 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:05:13 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 09:05:13 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:07:03 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 09:07:03 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:08:54 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 09:08:54 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:10:43 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# 0:     3718.6960: -1.20000  1.00000  1.00000  1.00000  1.30000 0.300000
# [1] "2023-09-30 09:10:43 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:12:33 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22428"
# [1] "2023-09-30 09:12:33 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:14:24 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22655"
# [1] "2023-09-30 09:14:24 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:16:15 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22426"
# [1] "2023-09-30 09:16:15 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:18:04 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22385"
# [1] "2023-09-30 09:18:04 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:19:54 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22464"
# [1] "2023-09-30 09:19:54 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:21:43 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22485"
# 1:     3681.2243: -1.69075 0.948139  1.00000  1.02489  1.33680 0.367106
# [1] "2023-09-30 09:21:43 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:23:35 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# [1] "2023-09-30 09:23:35 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:25:27 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27727"
# [1] "2023-09-30 09:25:27 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:27:18 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27739"
# [1] "2023-09-30 09:27:18 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:29:11 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# [1] "2023-09-30 09:29:11 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:30:57 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27749"
# [1] "2023-09-30 09:30:57 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:32:44 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27757"
# 2:     3672.2773: -2.06343  1.03255  1.00000 0.780039  1.45334 0.541597
# [1] "2023-09-30 09:32:44 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:41:17 ... ending   pcubature -- tol=0.1 -- ret.val is: 3678.6125"
# 3:     3672.2773: -2.06343  1.03255  1.00000 0.780039  1.45334 0.541597
# [1] "2023-09-30 09:41:17 ... starting pcubature for bivariate normal or cauchy"
# [1] "2023-09-30 09:43:11 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p      alpha       var1       var2     corr12
# -2.0634324  1.0325505  1.0000000  0.7800393  1.4533409  0.5415967
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 3672.277
#
# $fevals
# function
# 4
#
# $gevals
# gradient
# 15
#
# $nitns
# [1] 3
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 2539.906
#
# Assemble the answers
# Intercept      b_p alpha      var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb -2.063432 1.032551     1 0.7800393 1.453341 0.5415967 3672.277      4     15     3        0   NA   NA 2539.906



## test the newly added bivariate-t-corr locking df=1
## about 50 minutes to run.
(rand.int.rand.slopes.nonzero.corr.CUBA.tcorr <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=1),

                   pmu = c(Intercept=-1.2, b_p=1),
                   pmix=c(df=1.0, var1=1, var2=1.3, corr12= 0.30),

                   p_uppb = c(  0,   2, 1.0, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-t-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
# > (rand.int.rand.slopes.nonzero.corr.CUBA.tcorr <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ pt(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, df=1),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1),
#                           +                    pmix=c(df=1.0, var1=1, var2=1.3, corr12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2, 1.0, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-t-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    tol.pcubature = 0.1,
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2
#                           +     )
#    + )
# fn is  fn1
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  1   <? 1   <? 1     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-09-30 11:23:55 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:25:48 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.3467875
# Method:  nlminb
# [1] "2023-09-30 11:25:48 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:27:38 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 11:27:38 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:29:29 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69603"
# [1] "2023-09-30 11:29:29 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:31:25 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 11:31:25 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:33:16 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 11:33:16 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:35:04 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-30 11:35:04 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:36:54 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# 0:     3718.6960: -1.20000  1.00000  1.00000  1.00000  1.30000 0.300000
# [1] "2023-09-30 11:36:54 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:38:42 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22428"
# [1] "2023-09-30 11:38:42 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:40:30 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22655"
# [1] "2023-09-30 11:40:30 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:42:19 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22426"
# [1] "2023-09-30 11:42:19 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:44:08 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22385"
# [1] "2023-09-30 11:44:08 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:45:56 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22464"
# [1] "2023-09-30 11:45:56 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:47:45 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22485"
# 1:     3681.2243: -1.69075 0.948139  1.00000  1.02489  1.33680 0.367106
# [1] "2023-09-30 11:47:45 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:49:35 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# [1] "2023-09-30 11:49:35 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:51:52 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27727"
# [1] "2023-09-30 11:51:52 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 11:59:09 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27739"
# [1] "2023-09-30 11:59:09 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:01:05 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# [1] "2023-09-30 12:01:05 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:03:02 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27749"
# [1] "2023-09-30 12:03:02 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:04:51 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27757"
# 2:     3672.2773: -2.06343  1.03255  1.00000 0.780039  1.45334 0.541597
# [1] "2023-09-30 12:04:51 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:13:26 ... ending   pcubature -- tol=0.1 -- ret.val is: 3678.6125"
# 3:     3672.2773: -2.06343  1.03255  1.00000 0.780039  1.45334 0.541597
# [1] "2023-09-30 12:13:26 ... starting pcubature for bivariate-t-corr"
# [1] "2023-09-30 12:15:17 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p         df       var1       var2     corr12
# -2.0634324  1.0325505  1.0000000  0.7800393  1.4533409  0.5415967
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 3672.277
#
# $fevals
# function
# 4
#
# $gevals
# gradient
# 15
#
# $nitns
# [1] 3
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 2600.974
#
# Assemble the answers
# Intercept      b_p df      var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb -2.063432 1.032551  1 0.7800393 1.453341 0.5415967 3672.277      4     15     3        0   NA   NA 2600.974



## 3hr20min
(rand.int.rand.slopes.nonzero.corr.CUBA <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),

                   pmu = c(Intercept=-1.2, b_p=1, alpha=1.0),
                   pmix=c(alpha=1.0, var1=1, var2=1.3, corr12= 0.30),

                   p_uppb = c(  0,   2, 1.0, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-subgauss-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   tol.pcubature = 0.1,
                   abs.tol.nlminb = 1e-2,
                   xf.tol.nlminb =  1e-2,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2
    )
)
##
# > (rand.int.rand.slopes.nonzero.corr.CUBA <-
#      +   gnlrim::gnlrem(y=ybind,
#                         +                  mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),
#                         +
#                           +                  pmu = c(Intercept=-1.2, b_p=1, alpha=1.0),
#                         +                  pmix=c(alpha=1.0, var1=1, var2=1.3, corr12= 0.30),
#                         +
#                           +                  p_uppb = c(  0,   2, 1.0, 4.00, 4.00, 0.90),
#                         +                  p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
#                         +                  distribution="binomial",
#                         +                  nest=id,
#                         +                  random=c("rand1", "rand2"),
#                         +                  mixture="bivariate-subgauss-corr",
#                         +                  ooo=TRUE,
#                         +                  compute_hessian = FALSE,
#                         +                  compute_kkt = FALSE,
#                         +                  trace=1,
#                         +                  method='nlminb',
#                         +                  int2dmethod="cuba",
#                         +                  tol.pcubature = 0.1,
#                         +                  abs.tol.nlminb = 1e-2,
#                         +                  xf.tol.nlminb =  1e-2,
#                         +                  x.tol.nlminb =   1e-2,
#                         +                  rel.tol.nlminb = 1e-2
#                         +   )
#    + )
# fn is  fn
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  1   <? 1   <? 1     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2023-09-23 11:17:12 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:21:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.3467875
# Method:  nlminb
# [1] "2023-09-23 11:21:37 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:26:02 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-23 11:26:02 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:30:31 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69603"
# [1] "2023-09-23 11:30:31 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:34:58 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-23 11:34:58 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:39:24 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-23 11:39:24 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:43:49 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# [1] "2023-09-23 11:43:49 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:48:14 ... ending   pcubature -- tol=0.1 -- ret.val is: 3718.69602"
# 0:     3718.6960: -1.20000  1.00000  1.00000  1.00000  1.30000 0.300000
# [1] "2023-09-23 11:48:14 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:52:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22428"
# [1] "2023-09-23 11:52:37 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 11:57:01 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22655"
# [1] "2023-09-23 11:57:01 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:01:26 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22426"
# [1] "2023-09-23 12:01:26 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:05:42 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22385"
# [1] "2023-09-23 12:05:42 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:09:58 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22464"
# [1] "2023-09-23 12:09:58 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:14:14 ... ending   pcubature -- tol=0.1 -- ret.val is: 3681.22485"
# 1:     3681.2243: -1.69075 0.948139  1.00000  1.02489  1.33680 0.367106
# [1] "2023-09-23 12:14:14 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:18:31 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# [1] "2023-09-23 12:18:31 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:22:49 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27727"
# [1] "2023-09-23 12:22:49 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:27:07 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27739"
# [1] "2023-09-23 12:27:07 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:31:25 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# [1] "2023-09-23 12:31:25 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:35:42 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27749"
# [1] "2023-09-23 12:35:42 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 12:40:02 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27757"
# 2:     3672.2773: -2.06343  1.03255  1.00000 0.780039  1.45334 0.541597
# [1] "2023-09-23 12:40:02 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 14:40:10 ... ending   pcubature -- tol=0.1 -- ret.val is: 3678.6125"
# 3:     3672.2773: -2.06343  1.03255  1.00000 0.780039  1.45334 0.541597
# [1] "2023-09-23 14:40:10 ... starting pcubature for bivariate-subgauss-corr"
# [1] "2023-09-23 14:44:54 ... ending   pcubature -- tol=0.1 -- ret.val is: 3672.27732"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p      alpha       var1       var2     corr12
# -2.0634324  1.0325505  1.0000000  0.7800393  1.4533409  0.5415967
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 3672.277
#
# $fevals
# function
# 4
#
# $gevals
# gradient
# 15
#
# $nitns
# [1] 3
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 4469.231
#
# Assemble the answers
# Intercept      b_p alpha      var1     var2    corr12    value fevals
# nlminb -2.063432 1.032551     1 0.7800393 1.453341 0.5415967 3672.277      4
# gevals niter convcode kkt1 kkt2    xtime
# nlminb     15     3        0   NA   NA 4469.231










## takes 3.5 hours to run.
(rand.int.rand.slopes.nonzero.corr.CUBA <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),

                   pmu = c(Intercept=-1.2, b_p=1, alpha=1.0),
                   pmix=c(alpha=1.7, var1=1, var2=1.3, corr12= 0.30),

                   p_uppb = c(  0,   2, 1.0, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2, 1.0, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-subgauss-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba"
    )
)
# fn is  fn
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -4   <? -0.95   <? 0     In Bounds
# par[ 2 ]:  -2   <? 0.55   <? 2     In Bounds
# par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 4 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.3467875
# Method:  nlminb
#  0:     2712.2660: -0.950000 0.550000  1.00000  1.30000 0.300000
#  1:     2551.8588: -1.91097 0.564613  1.24340  1.36070 0.415727
#  2:     2540.4907: -2.08561 0.766148  1.12182  1.43309 0.534328
#  3:     2539.6569: -1.91468 0.976685  1.00080  1.50592 0.642532
#  4:     2537.8782: -2.02461 0.944054 0.972437  1.50899 0.618941
#  5:     2537.8313: -1.98169 0.912124 0.923186  1.56179 0.538668
#  6:     2537.5632: -2.01604 0.896569 0.949136  1.58691 0.568717
#  7:     2537.5229: -1.99913 0.892145 0.958570  1.59318 0.574194
#  8:     2537.4951: -2.01024 0.882361 0.956124  1.60724 0.580601
#  9:     2537.4433: -1.99579 0.863435 0.952010  1.64238 0.586637
# 10:     2537.4051: -2.00440 0.859836 0.952438  1.68363 0.578498
# 11:     2537.3828: -1.99446 0.848237 0.931713  1.71592 0.590764
# 12:     2537.3725: -1.99310 0.839836 0.944576  1.72247 0.586854
# 13:     2537.3615: -1.99591 0.841539 0.941970  1.73897 0.584327
# 14:     2537.3502: -1.99413 0.831356 0.945324  1.77149 0.587329
# 15:     2537.3444: -1.99242 0.834011 0.940837  1.78887 0.581282
# 16:     2537.3435: -1.99459 0.832976 0.940084  1.78959 0.581511
# 17:     2537.3430: -1.99327 0.831268 0.940081  1.79106 0.581831
# 18:     2537.3399: -1.99139 0.821393 0.932831  1.81239 0.584840
# 19:     2537.3387: -1.99184 0.819531 0.937044  1.83662 0.581810
# 20:     2537.3384: -1.99164 0.819911 0.935593  1.83029 0.582834
# 21:     2537.3384: -1.99164 0.819914 0.935584  1.83022 0.582816
#
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -1.9916400  0.8199140  0.9355840  1.8302160  0.5828159
#
# $message
# [1] "relative convergence (4)"

# Assemble the answers
#        Intercept      b_p     var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb  -1.99164 0.819914 0.935584 1.830216 0.5828159 2537.338     32    143    21        0   NA   NA



## takes about 3hr45min to run:
(rand.int.rand.slopes.nonzero.corr.CUBA.MARGINAL <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm(
                     (Intercept + period_numeric*b_p)*
                       sqrt(1 + var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric ) +
                       rand1 + period_numeric*rand2
                   ),
                   pmu = c(Intercept=-0.95, b_p=0.55, var1=1, var2=1, corr12= 0.20),
                   pmix=c(var1=1, var2=1, corr12= 0.20),
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
                   method='nlminb',
                   int2dmethod="cuba"
    )
)
## [1] "2022-02-16 21:42:37 EST"
# fn is  fn
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -4   <? -0.95   <? 0     In Bounds
# par[ 2 ]:  -2   <? 0.55   <? 2     In Bounds
# par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? 0.2   <? 0.9     In Bounds
# [1] "2022-02-16 21:43:24 EST"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.69897   log bounds ratio= 0.3467875
# Method:  nlminb
#  0:     2624.4797: -0.950000 0.550000  1.00000  1.00000 0.200000
#  1:     2562.9606: -1.37963 0.615494  1.15864  1.03380 0.262517
#  2:     2552.5119: -1.23506  1.05587  1.16743  1.08807 0.297671
#  3:     2542.3105: -1.44351  1.00734  1.20491  1.10583 0.347031
#  4:     2539.4509: -1.39190  1.05166  1.12312  1.19576 0.521837
#  5:     2539.0653: -1.47390  1.06191 0.950419  1.30500 0.558917
#  6:     2538.9093: -1.39239  1.08395 0.930294  1.36602 0.555530
#  7:     2537.8924: -1.41833  1.06487 0.987916  1.44436 0.583334
#  8:     2537.7425: -1.41637  1.09756 0.968858  1.54184 0.565155
#  9:     2537.5635: -1.42303  1.06617 0.949298  1.64112 0.565285
# 10:     2537.5148: -1.41318  1.08047 0.955722  1.66005 0.582748
# 11:     2537.5036: -1.42424  1.07913 0.956294  1.66253 0.584436
# 12:     2537.4869: -1.42027  1.07978 0.951819  1.67225 0.582759
# 13:     2537.4773: -1.42564  1.08940 0.945028  1.69029 0.589160
# 14:     2537.4556: -1.42156  1.08470 0.943423  1.71247 0.588605
# 15:     2537.4479: -1.42067  1.07906 0.939057  1.73378 0.583282
# 16:     2537.4390: -1.42280  1.08308 0.945738  1.75541 0.582606
# 17:     2537.4311: -1.42020  1.08397 0.944559  1.77821 0.580317
# 18:     2537.4289: -1.42185  1.08430 0.937058  1.78136 0.582546
# 19:     2537.4275: -1.42060  1.08443 0.937152  1.78925 0.585739
# 20:     2537.4258: -1.42065  1.08632 0.936922  1.80843 0.581730
# 21:     2537.4255: -1.42129  1.08492 0.936574  1.81418 0.582024
# 22:     2537.4254: -1.42155  1.08627 0.934905  1.81566 0.583622
# 23:     2537.4253: -1.42110  1.08577 0.935267  1.81651 0.583067
# 24:     2537.4253: -1.42124  1.08580 0.935359  1.81593 0.583069
# 25:     2537.4253: -1.42124  1.08581 0.935325  1.81599 0.583094
# [1] "2022-02-17 01:23:51 EST"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -1.4212356  1.0858135  0.9353249  1.8159932  0.5830943
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2537.425
#
# $fevals
# function
# 34
#
# $gevals
# gradient
# 164
#
# $nitns
# [1] 25
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 12181.54
#
# Assemble the answers
#        Intercept      b_p      var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb -1.421236 1.085814 0.9353249 1.815993 0.5830943 2537.425     34    164    25        0   NA   NA
# xtime
# nlminb 12181.54
glm.fit <-
  glm(cbind(r, n_r) ~ x1, summed_binom_dat, family=binomial(link="probit"))
summary(glm.fit)


phi_x <- function(x, var1=0.935584, var2=1.830216, corr12=0.5828159){

  1/sqrt(1 + var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x )

}

xvalueset <- c(0, sort(unique(period_numeric)), 1, 10)
cbind(
  xvalue = xvalueset,

  "phi(xvalue)" = phi_x(x=xvalueset),

  beta0_c = -1.99164,
  beta1_c = 0.819914,

  beta0_m =  -1.99164*  phi_x(x=xvalueset),
  beta1_m =  0.819914 *  phi_x(x=xvalueset)
)

xdomain <- seq(-2,4,0.05)
plot(xdomain, phi_x(xdomain), ylim=c(0,1))

library(data.table)
dat.interim <- data.table::copy(summed_binom_dat[, r/(r+n_r), by=c("id","period_numeric")])

## do it either way:
dat.interim2 <- data.table::copy(dat.interim[order(period_numeric),mean(V1), by=period_numeric])
summed_binom_dat[, mean(r/(r+n_r)), by=c("period_numeric")]

dat.interim2[,marginal.fit:=pnorm( -1.4212364+1.085814*period_numeric)]
dat.interim2[,conditional.fit:=pnorm( -1.99164+0.819914*period_numeric)]
setnames(dat.interim2, "V1", "empirical.avg")
setcolorder(dat.interim2, c("period_numeric", "conditional.fit","empirical.avg","marginal.fit" ))

dat.interim2[, glm.fit:=pnorm(-1.42435 + 1.10648 *period_numeric)]

print(dat.interim2, digits=2)
plot(x1, r/(r+n_r),col="grey", ylim=c(0,1),xlim=range(xdomain))
lines(xdomain, pnorm( -1.99164+0.819914*xdomain), ylim=c(0,1))
lines(xdomain, pnorm((-1.99164+0.819914*xdomain)*phi_x(xdomain)), col="red")
lines(xdomain, pnorm( -1.4212364+1.085814*xdomain), ylim=c(0,1),col="orange")
lines(xdomain, phi_x(xdomain), ylim=c(0,1), lty=2, col='cyan')
text(3.6,0.52, "(b0c+b1c*x)*phi(x)", col="red")
text(3.6,0.68, "(b0c+b1c*x)       ", col="black")
text(3.6,0.95, "(b0m+b1m*x)       ", col="orange")
text(3.6,0.19, "            phi(x)", col="cyan")





########################################$$$$#########
## 2022-02-17E                                     ##
## END two random parameters: CAUCHY MARGINAL COEFF##
#####################################################

########################################$$$$#########
## 2022-02-17D                                     ##
## START two random parameters: MARGINAL COEFF     ##
#####################################################
## next step: do 2022-02-16B with mvsubgauss freely
## estimating alpha.  Then you are ready to
## write gapr paper


##library(mvsubgaussPD)
library(mvpd)
library(libstable4u)
library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Stabit-BIVARIATE_STABLE", alpha=1.7, gamma1=1, gamma2=2, gamma12=-4, delta=1){

  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
  }
  if(mrim=="Stabit-BIVARIATE_STABLE"){
    print("STABLE for the win")
    print(paste0("alpha set to ", alpha))
    G <- function(n, v1, v2, rho){mvsubgaussPD::rmvsubgaussPD(n=n, alpha=alpha, Q=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) stable_cdf(x, c(alpha,0,1,0))
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
set.seed(39)##set.seed(1709)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1)
#  sim_mrim_data(10,10, J=100, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat
## glmer -- only random intercept
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)




attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## SSS
(rand.int <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1, c(alpha, 0, 1, 0)),
                   pmu = c(Intercept=-1.2, b_p=1, alpha=1.7),
                   pmix=c(alpha=1.7, scl=0.25),
                   p_uppb = c(  0,   2, 1.9, 4.00),
                   p_lowb = c( -4,  -2, 1.4, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-scl",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)

## excuse the 1^alpha; had to do it to get alpha to appear before phi in `mu` argument.
(rand.int.marg.phi <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2((Intercept + period_numeric*b_p)*1^alpha/phi + rand1, c(alpha, 0, 1, 0)),
                   pmu = c(Intercept=-0.95, b_p=0.55,  alpha=1.7, phi=0.3),
                   pmix=c(alpha=1.7, phi=0.3),
                   p_uppb = c(  0,   2, 1.9, 0.99),
                   p_lowb = c( -4,  -2, 1.4, 0.01),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="libstable4u-subgauss-phi",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)

#> This is phi from the marginal phi model:
rand.int.marg.phi[4]
#> This is phi calculated from conditional model:
1/(1 + rand.int[4]^rand.int[3])^(1/rand.int[3])


## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)



## takes X.XX hours to run.
(rand.int.rand.slopes.nonzero.corr.CUBA <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),

                   #pmu = c(Intercept=-1.2, b_p=1, alpha=1.7),
                   #pmix=c(alpha=1.7, var1=1, var2=1.3, corr12= 0.30),

                   p_uppb = c(  0,   2, 1.90, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2, 1.40, 0.05, 0.05,-0.90),

                   ##pmu = c(Intercept= -1.362697, b_p=0.9523764, alpha=1.747751),
                   ## pmu = c(Intercept= -1.96, b_p=0.90, alpha=1.747751),
                   #pmu = c(Intercept= -2, b_p=1, alpha=1.747751),
                   #pmix=c( alpha=1.747751, var1=0.8093904, var2=1.305351, corr12= 0.6295481),

                   pmu = c(Intercept= 0, b_p=-2, alpha=1.7),
                   pmix=c( alpha=1.7, var1=1, var2=2, corr12= -0.2),

                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-subgauss-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                  int2dmethod="cuba",
                   abs.tol.nlminb = 0,
                   xf.tol.nlminb =  1e-1,
                   x.tol.nlminb =   1e-2,
                   rel.tol.nlminb = 1e-2,
                  tol.pcubature = 0.1
    )
)
# > (rand.int.rand.slopes.nonzero.corr.CUBA <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),
#                           +
#                             +
#                             +                    p_uppb = c(  0,   2, 1.90, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2, 1.40, 0.05, 0.05,-0.90),
#                           +
#                             +
#                             +                    pmu = c(Intercept= 0, b_p=-2, alpha=1.7),
#                           +                    pmix=c( alpha=1.7, var1=1, var2=2, corr12= -0.2),
#                           +
#                             +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-subgauss-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                   int2dmethod="cuba",
#                           +                    abs.tol.nlminb = 0,
#                           +                    xf.tol.nlminb =  1e-1,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2,
#                           +     )
#    + )
# fn is  fn
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? 0   <? 0     In Bounds
# par[ 2 ]:  -2   <? -2   <? 2     In Bounds
# par[ 3 ]:  1.4   <? 1.7   <? 1.9     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 2   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? -0.2   <? 0.9     In Bounds
# [1] "2022-08-10 19:21:01 ... starting pcubature"
# [1] "2022-08-10 21:22:53 ... ending   pcubature -- tol=0.1 -- ret.val is: 3700.26629"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 1   log bounds ratio= 0.90309
# Method:  nlminb
# [1] "2022-08-10 21:22:53 ... starting pcubature"
# [1] "2022-08-10 23:23:43 ... ending   pcubature -- tol=0.1 -- ret.val is: 3700.26629"
# [1] "2022-08-10 23:23:43 ... starting pcubature"
# [1] "2022-08-11 01:23:54 ... ending   pcubature -- tol=0.1 -- ret.val is: 3700.26629"
# [1] "2022-08-11 01:23:54 ... starting pcubature"
# [1] "2022-08-11 03:24:37 ... ending   pcubature -- tol=0.1 -- ret.val is: 3700.26629"
# [1] "2022-08-11 03:24:37 ... starting pcubature"
# [1] "2022-08-11 05:25:18 ... ending   pcubature -- tol=0.1 -- ret.val is: 3700.26629"
# [1] "2022-08-11 05:25:18 ... starting pcubature"
# Timing stopped at: 5.914e+04 2.147e+04 3.35e+04
# >
# > (rand.int.rand.slopes.nonzero.corr.CUBA <-
#      +     gnlrim::gnlrem(y=ybind,
#                           +                    mu = ~ stable_cdf2(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2, c(alpha, 0, 1, 0)),
#                           +
#                             +                    pmu = c(Intercept=-1.2, b_p=1, alpha=1.7),
#                           +                    pmix=c(alpha=1.7, var1=1, var2=1.3, corr12= 0.30),
#                           +
#                             +                    p_uppb = c(  0,   2, 1.90, 4.00, 4.00, 0.90),
#                           +                    p_lowb = c( -4,  -2, 1.40, 0.05, 0.05,-0.90),
#                           +                    distribution="binomial",
#                           +                    nest=id,
#                           +                    random=c("rand1", "rand2"),
#                           +                    mixture="bivariate-subgauss-corr",
#                           +                    ooo=TRUE,
#                           +                    compute_hessian = FALSE,
#                           +                    compute_kkt = FALSE,
#                           +                    trace=1,
#                           +                    method='nlminb',
#                           +                    int2dmethod="cuba",
#                           +                    abs.tol.nlminb = 1e-2,
#                           +                    xf.tol.nlminb =  1e-2,
#                           +                    x.tol.nlminb =   1e-2,
#                           +                    rel.tol.nlminb = 1e-2,
#                           +     )
#    + )
# [1] "2022-02-18 14:48:46 ... starting pcubature"
# [1] "2022-02-18 14:56:13 ... ending   pcubature tol=0.01"
# [1] 6
# Intercept       b_p     alpha      var1      var2    corr12
# -1.2       1.0       1.7       1.0       1.3       0.3
# [1] 3432.649
# fn is  fn
# Looking for method =  nlminb
# Function has  6  arguments
# par[ 1 ]:  -4   <? -1.2   <? 0     In Bounds
# par[ 2 ]:  -2   <? 1   <? 2     In Bounds
# par[ 3 ]:  1.4   <? 1.7   <? 1.9     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 6 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# [1] "2022-02-18 14:56:13 ... starting pcubature"
# [1] "2022-02-18 15:03:46 ... ending   pcubature tol=0.01"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.7533277   log bounds ratio= 0.90309
# Method:  nlminb
# [1] "2022-02-18 15:03:47 ... starting pcubature"
# [1] "2022-02-18 15:11:48 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 15:11:48 ... starting pcubature"
# [1] "2022-02-18 15:19:43 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 15:19:43 ... starting pcubature"
# [1] "2022-02-18 15:27:51 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 15:27:51 ... starting pcubature"
# [1] "2022-02-18 15:35:57 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 15:35:57 ... starting pcubature"
# [1] "2022-02-18 15:43:37 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 15:43:37 ... starting pcubature"
# [1] "2022-02-18 15:51:14 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 15:51:14 ... starting pcubature"
# [1] "2022-02-18 15:59:29 ... ending   pcubature tol=0.01"
# 0:     3432.6490: -1.20000  1.00000  1.70000  1.00000  1.30000 0.300000
# [1] "2022-02-18 15:59:29 ... starting pcubature"
# [1] "2022-02-18 16:08:36 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 16:08:36 ... starting pcubature"
# [1] "2022-02-18 16:17:37 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 16:17:37 ... starting pcubature"
# [1] "2022-02-18 16:27:10 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 16:27:10 ... starting pcubature"
# [1] "2022-02-18 16:35:14 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 16:35:14 ... starting pcubature"
# [1] "2022-02-18 16:43:16 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 16:43:16 ... starting pcubature"
# [1] "2022-02-18 16:51:24 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 16:51:24 ... starting pcubature"
# [1] "2022-02-18 16:59:27 ... ending   pcubature tol=0.01"
# 1:     3400.0784: -2.09980 0.908341  1.90000  1.16001  1.33176 0.517222
# [1] "2022-02-18 16:59:27 ... starting pcubature"
# [1] "2022-02-18 17:05:05 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:05:05 ... starting pcubature"
# [1] "2022-02-18 17:11:53 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:11:53 ... starting pcubature"
# [1] "2022-02-18 17:18:41 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:18:41 ... starting pcubature"
# [1] "2022-02-18 17:25:41 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:25:41 ... starting pcubature"
# [1] "2022-02-18 17:32:31 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:32:31 ... starting pcubature"
# [1] "2022-02-18 17:39:19 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:39:19 ... starting pcubature"
# [1] "2022-02-18 17:46:06 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 17:46:06 ... starting pcubature"
# [1] "2022-02-18 17:52:51 ... ending   pcubature tol=0.01"
# 2:     3381.6564: -2.05560 0.907948  1.66244  1.11724  1.32523 0.499059
# [1] "2022-02-18 17:52:51 ... starting pcubature"
# [1] "2022-02-18 18:00:18 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:00:18 ... starting pcubature"
# [1] "2022-02-18 18:07:44 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:07:44 ... starting pcubature"
# [1] "2022-02-18 18:15:11 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:15:11 ... starting pcubature"
# [1] "2022-02-18 18:22:35 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:22:35 ... starting pcubature"
# [1] "2022-02-18 18:30:01 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:30:01 ... starting pcubature"
# [1] "2022-02-18 18:37:44 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:37:44 ... starting pcubature"
# [1] "2022-02-18 18:45:23 ... ending   pcubature tol=0.01"
# 3:     3379.3622: -1.93967 0.910037  1.73290 0.918896  1.30546 0.548499
# [1] "2022-02-18 18:45:23 ... starting pcubature"
# [1] "2022-02-18 18:53:08 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 18:53:08 ... starting pcubature"
# [1] "2022-02-18 19:00:13 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 19:00:13 ... starting pcubature"
# [1] "2022-02-18 19:07:18 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 19:07:18 ... starting pcubature"
# [1] "2022-02-18 19:14:22 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 19:14:22 ... starting pcubature"
# [1] "2022-02-18 19:21:26 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 19:21:26 ... starting pcubature"
# [1] "2022-02-18 19:28:20 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 19:28:20 ... starting pcubature"
# [1] "2022-02-18 19:35:14 ... ending   pcubature tol=0.01"
# [1] "2022-02-18 19:35:14 ... starting pcubature"
# [1] "2022-02-18 19:42:09 ... ending   pcubature tol=0.01"
# 4:     3379.1129: -1.96450 0.893471  1.74994 0.895177  1.30239 0.556624
# [1] "2022-02-18 19:42:09 ... starting pcubature"
# [1] "2022-02-18 19:50:00 ... ending   pcubature tol=0.01"
# 5:     3379.1129: -1.96450 0.893471  1.74994 0.895177  1.30239 0.556624
# [1] "2022-02-18 19:50:00 ... starting pcubature"
# nlminb function evaluation failure
# Post processing for method  nlminb
# Save results from method  nlminb
# $fevals
# [1] NA
#
# $convcode
# [1] 9999
#
# $value
# [1] 8.988466e+307
#
# $par
# [1] NA NA NA NA NA NA
#
# $nitns
# [1] NA
#
# $gevals
# [1] NA
#
# $message
# [1] "nlminb failure"
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 31125.53
#
# Assemble the answers
# Intercept b_p alpha var1 var2 corr12         value fevals gevals niter convcode kkt1 kkt2    xtime
# nlminb        NA  NA    NA   NA   NA     NA 8.988466e+307     NA     NA    NA     9999   NA   NA 31125.53

## calculate bm_1, bm_0:
 0.893471* 1/(1+  (0.895177 + 1.30239*0.33^2 + 2*0.5566*sqrt(0.895177*1.30239)*0.33)^(1.74994/2)     )^(1/1.74994)
-1.96450 * 1/(1+  (0.895177 + 1.30239*0.00^2 + 2*0.5566*sqrt(0.895177*1.30239)*0.00)^(1.74994/2)     )^(1/1.74994)

 ## takes about XhrXXmin to run:
(rand.int.rand.slopes.nonzero.corr.CUBA.MARGINAL <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ stable_cdf2(
                     (Intercept + period_numeric*b_p)* 1^alpha *
                       (1^alpha + (var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric)^(alpha/2) )^(1/alpha) +
                       rand1 + period_numeric*rand2
                     , c(alpha, 0, 1, 0)),
#                   pmu = c(Intercept= -1.268806, b_p=0.7400231, alpha=1.677604, var1=0.8507984, var2=1.301004, corr12= 0.5395131),
#                   pmix=c(alpha=1.677604, var1=0.8507984, var2=1.301004 , corr12= 0.5395131),
                   pmu = c(Intercept= -1.362697, b_p=0.9523764, alpha=1.747751, var1=0.8093904, var2=1.305351, corr12= 0.6295481),
                   pmix=c( alpha=1.747751, var1=0.8093904, var2=1.305351, corr12= 0.6295481),
                   p_uppb = c(  0,   2, 1.90, 4.00, 4.00, 0.90),
                   p_lowb = c( -4,  -2, 1.40, 0.05, 0.05,-0.90),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1", "rand2"),
                   mixture="bivariate-subgauss-corr",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
                   int2dmethod="cuba",
                   abs.tol.nlminb = 0,
                   xf.tol.nlminb =  1e-3,
                   x.tol.nlminb =   1.5e-8,
                   rel.tol.nlminb = 1e-1,
    )
)
 # >  ## takes about XhrXXmin to run:
 #   > (rand.int.rand.slopes.nonzero.corr.CUBA.MARGINAL <-
 #        +     gnlrim::gnlrem(y=ybind,
 #                             +                    mu = ~ stable_cdf2(
 #                               +                      (Intercept + period_numeric*b_p)* 1^alpha *
 #                                 +                        (1^alpha + (var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric)^(alpha/2) )^(1/alpha) +
 #                                 +                        rand1 + period_numeric*rand2
 #                               +                      , c(alpha, 0, 1, 0)),
 #                             +                    pmu = c(Intercept= -1.35819, b_p=0.5456027, alpha=1.74994, var1=0.895177, var2=1.30239, corr12= 0.556624),
 #                             +                    pmix=c(alpha=1.74994, var1=0.895177, var2=1.30239, corr12= 0.556624),
 #                             +                    p_uppb = c(  0,   2, 1.90, 4.00, 4.00, 0.90),
 #                             +                    p_lowb = c( -4,  -2, 1.40, 0.05, 0.05,-0.90),
 #                             +                    distribution="binomial",
 #                             +                    nest=id,
 #                             +                    random=c("rand1", "rand2"),
 #                             +                    mixture="bivariate-subgauss-corr",
 #                             +                    ooo=TRUE,
 #                             +                    compute_hessian = FALSE,
 #                             +                    compute_kkt = FALSE,
 #                             +                    trace=1,
 #                             +                    method='nlminb',
 #                             +                    int2dmethod="cuba",
 #                             +                    abs.tol.nlminb = 1e-1,
 #                             +                    xf.tol.nlminb =  1e-1,
 #                             +                    x.tol.nlminb =   1e-1,
 #                             +                    rel.tol.nlminb = 1e-1,
 #                             +     )
 #      + )
 # [1] "2022-02-18 22:21:05 ... starting pcubature"
 # [1] "2022-02-18 22:28:36 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] 6
 # Intercept        b_p      alpha       var1       var2     corr12
 # -1.3581900  0.5456027  1.7499400  0.8951770  1.3023900  0.5566240
 # [1] 3389.83
 # fn is  fn
 # Looking for method =  nlminb
 # Function has  6  arguments
 # par[ 1 ]:  -4   <? -1.35819   <? 0     In Bounds
 # par[ 2 ]:  -2   <? 0.5456027   <? 2     In Bounds
 # par[ 3 ]:  1.4   <? 1.74994   <? 1.9     In Bounds
 # par[ 4 ]:  0.05   <? 0.895177   <? 4     In Bounds
 # par[ 5 ]:  0.05   <? 1.30239   <? 4     In Bounds
 # par[ 6 ]:  -0.9   <? 0.556624   <? 0.9     In Bounds
 # [1] "2022-02-18 22:28:36 ... starting pcubature"
 # [1] "2022-02-18 22:36:18 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # Analytic gradient not made available.
 # Analytic Hessian not made available.
 # Scale check -- log parameter ratio= 0.5061466   log bounds ratio= 0.90309
 # Method:  nlminb
 # [1] "2022-02-18 22:36:18 ... starting pcubature"
 # [1] "2022-02-18 22:43:41 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] "2022-02-18 22:43:41 ... starting pcubature"
 # [1] "2022-02-18 22:51:08 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] "2022-02-18 22:51:08 ... starting pcubature"
 # [1] "2022-02-18 22:58:29 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] "2022-02-18 22:58:29 ... starting pcubature"
 # [1] "2022-02-18 23:05:54 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] "2022-02-18 23:05:54 ... starting pcubature"
 # [1] "2022-02-18 23:13:43 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] "2022-02-18 23:13:43 ... starting pcubature"
 # [1] "2022-02-18 23:20:44 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # [1] "2022-02-18 23:20:44 ... starting pcubature"
 # [1] "2022-02-18 23:27:47 ... ending   pcubature -- tol=0.01 -- ret.val is: 3389.83037"
 # 0:     3389.8304: -1.35819 0.545603  1.74994 0.895177  1.30239 0.556624
 # [1] "2022-02-18 23:27:47 ... starting pcubature"
 # [1] "2022-02-18 23:30:52 ... ending   pcubature -- tol=0.01 -- ret.val is: 3463.72179"
 # [1] "2022-02-18 23:30:52 ... starting pcubature"
 # [1] "2022-02-18 23:39:11 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81883"
 # [1] "2022-02-18 23:39:11 ... starting pcubature"
 # [1] "2022-02-18 23:47:23 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81746"
 # [1] "2022-02-18 23:47:23 ... starting pcubature"
 # [1] "2022-02-18 23:55:57 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81996"
 # [1] "2022-02-18 23:55:57 ... starting pcubature"
 # [1] "2022-02-19 00:03:32 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.83836"
 # [1] "2022-02-19 00:03:32 ... starting pcubature"
 # [1] "2022-02-19 00:10:49 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81799"
 # [1] "2022-02-19 00:10:49 ... starting pcubature"
 # [1] "2022-02-19 00:17:59 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81871"
 # [1] "2022-02-19 00:17:59 ... starting pcubature"
 # [1] "2022-02-19 00:25:05 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81831"
 # 1:     3383.8188: -1.26881 0.740023  1.67760 0.850798  1.30100 0.539513
 # [1] "2022-02-19 00:25:05 ... starting pcubature"
 # [1] "2022-02-19 00:30:42 ... ending   pcubature -- tol=0.01 -- ret.val is: 3414.27633"
 # 2:     3383.8188: -1.26881 0.740023  1.67760 0.850798  1.30100 0.539513
 # [1] "2022-02-19 00:30:42 ... starting pcubature"
 # [1] "2022-02-19 00:37:46 ... ending   pcubature -- tol=0.01 -- ret.val is: 3383.81883"
 # Post processing for method  nlminb
 # Save results from method  nlminb
 # $par
 # Intercept        b_p      alpha       var1       var2     corr12
 # -1.2688059  0.7400231  1.6776039  0.8507984  1.3010037  0.5395131
 #
 # $message
 # [1] "false convergence (8)"
 #
 # $convcode
 # [1] 1
 #
 # $value
 # [1] 3383.819
 #
 # $fevals
 # function
 # 4
 #
 # $gevals
 # gradient
 # 12
 #
 # $nitns
 # [1] 2
 #
 # $kkt1
 # [1] NA
 #
 # $kkt2
 # [1] NA
 #
 # $xtimes
 # user.self
 # 12701.6
 #
 # Assemble the answers
 # Intercept       b_p    alpha      var1     var2    corr12    value fevals
 # nlminb -1.268806 0.7400231 1.677604 0.8507984 1.301004 0.5395131 3383.819      4
 # gevals niter convcode kkt1 kkt2   xtime
 # nlminb     12     2        1   NA   NA 12701.6




glm.fit <-
  glm(cbind(r, n_r) ~ x1, summed_binom_dat, family=binomial(link="probit"))
summary(glm.fit)


phi_x <- function(x, var1=0.935584, var2=1.830216, corr12=0.5828159){

  1/sqrt(1 + var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x )

}

xvalueset <- c(0, sort(unique(period_numeric)), 1, 10)
cbind(
  xvalue = xvalueset,

  "phi(xvalue)" = phi_x(x=xvalueset),

  beta0_c = -1.99164,
  beta1_c = 0.819914,

  beta0_m =  -1.99164*  phi_x(x=xvalueset),
  beta1_m =  0.819914 *  phi_x(x=xvalueset)
)

xdomain <- seq(-2,4,0.05)
plot(xdomain, phi_x(xdomain), ylim=c(0,1))

library(data.table)
dat.interim <- data.table::copy(summed_binom_dat[, r/(r+n_r), by=c("id","period_numeric")])

## do it either way:
dat.interim2 <- data.table::copy(dat.interim[order(period_numeric),mean(V1), by=period_numeric])
summed_binom_dat[, mean(r/(r+n_r)), by=c("period_numeric")]

dat.interim2[,marginal.fit:=pnorm( -1.4212364+1.085814*period_numeric)]
dat.interim2[,conditional.fit:=pnorm( -1.99164+0.819914*period_numeric)]
setnames(dat.interim2, "V1", "empirical.avg")
setcolorder(dat.interim2, c("period_numeric", "conditional.fit","empirical.avg","marginal.fit" ))

dat.interim2[, glm.fit:=pnorm(-1.42435 + 1.10648 *period_numeric)]

print(dat.interim2, digits=2)
plot(x1, r/(r+n_r),col="grey", ylim=c(0,1),xlim=range(xdomain))
lines(xdomain, pnorm( -1.99164+0.819914*xdomain), ylim=c(0,1))
lines(xdomain, pnorm((-1.99164+0.819914*xdomain)*phi_x(xdomain)), col="red")
lines(xdomain, pnorm( -1.4212364+1.085814*xdomain), ylim=c(0,1),col="orange")
lines(xdomain, phi_x(xdomain), ylim=c(0,1), lty=2, col='cyan')
text(3.6,0.52, "(b0c+b1c*x)*phi(x)", col="red")
text(3.6,0.68, "(b0c+b1c*x)       ", col="black")
text(3.6,0.95, "(b0m+b1m*x)       ", col="orange")
text(3.6,0.19, "            phi(x)", col="cyan")



########################################$$$$#########
## 2022-02-17D                                     ##
## END two random parameters     MARGINAL COEFF    ##
#####################################################






########################################$$$$#########
## 2022-02-17C                                     ##
## START two random parameters: univariate dnorm   ##
#####################################################
# TRIED SOME stuff cannot get it to jive.  Not crucial
# to the overall mission.  Set aside.
########################################$$$$#########
## 2022-02-17C                                     ##
## END two random parameters: univariate dnorm    ##
#####################################################

########################################$$$$#########
## 2022-02-16B                                     ##
## START two random parameters: MARGINAL COEFF     ##
#####################################################
## see the hard coded picture at the end of this chunk!
## explains everything!

library(data.table)
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Probit-BIVARIATE_NORM", alpha=1.89, gamma=1.2, delta=1){

  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
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
set.seed(1709)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(200,200, J=100, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat
## glmer -- only random intercept
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)




attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## PPN
(rand.int <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm(Intercept + period_numeric*b_p + rand1),
                   pmu = c(Intercept=-0.95, b_p=0.55),
                   pmix=c(var1=1),
                   p_uppb = c(  0,   2, 4.00),
                   p_lowb = c( -4,  -2, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
# Assemble the answers
#        Intercept      b_p  var1   value fevals gevals niter convcode kkt1 kkt2  xtime
# nlminb -2.209654 1.618575 1.472 2571.28     16     43    11        0   NA   NA 35.234
(rand.int.marg.var1 <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm( (Intercept + period_numeric*b_p)*sqrt(1+var1) + rand1),
                   pmu = c(Intercept=-0.95, b_p=0.55, var1=1),
                   pmix=c(var1=1),
                   p_uppb = c(  0,   2, 4.00),
                   p_lowb = c( -4,  -2, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
# Assemble the answers
# Intercept      b_p     var1   value fevals gevals niter convcode kkt1 kkt2  xtime
# nlminb -1.405397 1.029446 1.471998 2571.28     20     70    14        0   NA   NA 59.092
(rand.int.marg.phi <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm( (Intercept + period_numeric*b_p)/phi + rand1),
                   pmu = c(Intercept=-0.95, b_p=0.55, phi=0.3),
                   pmix=c(phi=0.3),
                   p_uppb = c(  0,   2, 0.99),
                   p_lowb = c( -4,  -2, 0.01),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-phi",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)
# Assemble the answers
#        Intercept      b_p       phi   value fevals gevals niter convcode kkt1 kkt2  xtime
# nlminb -1.405398 1.029448 0.6360282 2571.28     34     97    24        0   NA   NA 80.871
1/sqrt(1+1.471998)
#> [1] 0.6360276


## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
## Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite Quadrature, nAGQ =  0)
# [glmerMod]
# Family: binomial  ( probit )
# Formula: cbind(r, n_r) ~ x1 + (x1 | id)
# Data: summed_binom_dat
# AIC       BIC    logLik  deviance  df.resid
# 5114.722  5141.611 -2552.361  5104.722      1595
# Random effects:
#   Groups Name        Std.Dev. Corr
# id     (Intercept) 0.9379
# x1          1.3202   0.60
# Number of obs: 1600, groups:  id, 400
# Fixed Effects:
#   (Intercept)           x1
# -1.888        0.907
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)
##  Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( probit )
# Formula: cbind(r, n_r) ~ x1 + (x1 | id)
# Data: summed_binom_dat
# AIC       BIC    logLik  deviance  df.resid
# 5109.936  5136.825 -2549.968  5099.936      1595
# Random effects:
#   Groups Name        Std.Dev. Corr
# id     (Intercept) 0.9423
# x1          1.3336   0.60
# Number of obs: 1600, groups:  id, 400
# Fixed Effects:
#   (Intercept)           x1
# -1.9973       0.7912

## takes 3.5 hours to run.
(rand.int.rand.slopes.nonzero.corr.CUBA <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
                   pmu = c(Intercept=-0.95, b_p=0.55),
                   pmix=c(var1=1, var2=1.3, corr12= 0.30),
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
                   method='nlminb',
                   int2dmethod="cuba"
    )
)
# fn is  fn
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -4   <? -0.95   <? 0     In Bounds
# par[ 2 ]:  -2   <? 0.55   <? 2     In Bounds
# par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 4 ]:  0.05   <? 1.3   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? 0.3   <? 0.9     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.6368221   log bounds ratio= 0.3467875
# Method:  nlminb
#  0:     2712.2660: -0.950000 0.550000  1.00000  1.30000 0.300000
#  1:     2551.8588: -1.91097 0.564613  1.24340  1.36070 0.415727
#  2:     2540.4907: -2.08561 0.766148  1.12182  1.43309 0.534328
#  3:     2539.6569: -1.91468 0.976685  1.00080  1.50592 0.642532
#  4:     2537.8782: -2.02461 0.944054 0.972437  1.50899 0.618941
#  5:     2537.8313: -1.98169 0.912124 0.923186  1.56179 0.538668
#  6:     2537.5632: -2.01604 0.896569 0.949136  1.58691 0.568717
#  7:     2537.5229: -1.99913 0.892145 0.958570  1.59318 0.574194
#  8:     2537.4951: -2.01024 0.882361 0.956124  1.60724 0.580601
#  9:     2537.4433: -1.99579 0.863435 0.952010  1.64238 0.586637
# 10:     2537.4051: -2.00440 0.859836 0.952438  1.68363 0.578498
# 11:     2537.3828: -1.99446 0.848237 0.931713  1.71592 0.590764
# 12:     2537.3725: -1.99310 0.839836 0.944576  1.72247 0.586854
# 13:     2537.3615: -1.99591 0.841539 0.941970  1.73897 0.584327
# 14:     2537.3502: -1.99413 0.831356 0.945324  1.77149 0.587329
# 15:     2537.3444: -1.99242 0.834011 0.940837  1.78887 0.581282
# 16:     2537.3435: -1.99459 0.832976 0.940084  1.78959 0.581511
# 17:     2537.3430: -1.99327 0.831268 0.940081  1.79106 0.581831
# 18:     2537.3399: -1.99139 0.821393 0.932831  1.81239 0.584840
# 19:     2537.3387: -1.99184 0.819531 0.937044  1.83662 0.581810
# 20:     2537.3384: -1.99164 0.819911 0.935593  1.83029 0.582834
# 21:     2537.3384: -1.99164 0.819914 0.935584  1.83022 0.582816
#
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -1.9916400  0.8199140  0.9355840  1.8302160  0.5828159
#
# $message
# [1] "relative convergence (4)"

# Assemble the answers
#        Intercept      b_p     var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb  -1.99164 0.819914 0.935584 1.830216 0.5828159 2537.338     32    143    21        0   NA   NA



## takes about 3hr45min to run:
(rand.int.rand.slopes.nonzero.corr.CUBA.MARGINAL <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm(
                     (Intercept + period_numeric*b_p)*
                       sqrt(1 + var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric ) +
                       rand1 + period_numeric*rand2
                   ),
                   pmu = c(Intercept=-0.95, b_p=0.55, var1=1, var2=1, corr12= 0.20),
                   pmix=c(var1=1, var2=1, corr12= 0.20),
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
                   method='nlminb',
                   int2dmethod="cuba"
    )
)
## [1] "2022-02-16 21:42:37 EST"
# fn is  fn
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -4   <? -0.95   <? 0     In Bounds
# par[ 2 ]:  -2   <? 0.55   <? 2     In Bounds
# par[ 3 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 4 ]:  0.05   <? 1   <? 4     In Bounds
# par[ 5 ]:  -0.9   <? 0.2   <? 0.9     In Bounds
# [1] "2022-02-16 21:43:24 EST"
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.69897   log bounds ratio= 0.3467875
# Method:  nlminb
#  0:     2624.4797: -0.950000 0.550000  1.00000  1.00000 0.200000
#  1:     2562.9606: -1.37963 0.615494  1.15864  1.03380 0.262517
#  2:     2552.5119: -1.23506  1.05587  1.16743  1.08807 0.297671
#  3:     2542.3105: -1.44351  1.00734  1.20491  1.10583 0.347031
#  4:     2539.4509: -1.39190  1.05166  1.12312  1.19576 0.521837
#  5:     2539.0653: -1.47390  1.06191 0.950419  1.30500 0.558917
#  6:     2538.9093: -1.39239  1.08395 0.930294  1.36602 0.555530
#  7:     2537.8924: -1.41833  1.06487 0.987916  1.44436 0.583334
#  8:     2537.7425: -1.41637  1.09756 0.968858  1.54184 0.565155
#  9:     2537.5635: -1.42303  1.06617 0.949298  1.64112 0.565285
# 10:     2537.5148: -1.41318  1.08047 0.955722  1.66005 0.582748
# 11:     2537.5036: -1.42424  1.07913 0.956294  1.66253 0.584436
# 12:     2537.4869: -1.42027  1.07978 0.951819  1.67225 0.582759
# 13:     2537.4773: -1.42564  1.08940 0.945028  1.69029 0.589160
# 14:     2537.4556: -1.42156  1.08470 0.943423  1.71247 0.588605
# 15:     2537.4479: -1.42067  1.07906 0.939057  1.73378 0.583282
# 16:     2537.4390: -1.42280  1.08308 0.945738  1.75541 0.582606
# 17:     2537.4311: -1.42020  1.08397 0.944559  1.77821 0.580317
# 18:     2537.4289: -1.42185  1.08430 0.937058  1.78136 0.582546
# 19:     2537.4275: -1.42060  1.08443 0.937152  1.78925 0.585739
# 20:     2537.4258: -1.42065  1.08632 0.936922  1.80843 0.581730
# 21:     2537.4255: -1.42129  1.08492 0.936574  1.81418 0.582024
# 22:     2537.4254: -1.42155  1.08627 0.934905  1.81566 0.583622
# 23:     2537.4253: -1.42110  1.08577 0.935267  1.81651 0.583067
# 24:     2537.4253: -1.42124  1.08580 0.935359  1.81593 0.583069
# 25:     2537.4253: -1.42124  1.08581 0.935325  1.81599 0.583094
# [1] "2022-02-17 01:23:51 EST"
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -1.4212356  1.0858135  0.9353249  1.8159932  0.5830943
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 2537.425
#
# $fevals
# function
# 34
#
# $gevals
# gradient
# 164
#
# $nitns
# [1] 25
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 12181.54
#
# Assemble the answers
#        Intercept      b_p      var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb -1.421236 1.085814 0.9353249 1.815993 0.5830943 2537.425     34    164    25        0   NA   NA
# xtime
# nlminb 12181.54
glm.fit <-
glm(cbind(r, n_r) ~ x1, summed_binom_dat, family=binomial(link="probit"))
summary(glm.fit)


phi_x <- function(x, var1=0.935584, var2=1.830216, corr12=0.5828159){

  1/sqrt(1 + var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x )

}

xvalueset <- c(0, sort(unique(period_numeric)), 1, 10)
cbind(
  xvalue = xvalueset,

  "phi(xvalue)" = phi_x(x=xvalueset),

  beta0_c = -1.99164,
  beta1_c = 0.819914,

  beta0_m =  -1.99164*  phi_x(x=xvalueset),
  beta1_m =  0.819914 *  phi_x(x=xvalueset)
)

xdomain <- seq(-2,4,0.05)
plot(xdomain, phi_x(xdomain), ylim=c(0,1))

library(data.table)
dat.interim <- data.table::copy(summed_binom_dat[, r/(r+n_r), by=c("id","period_numeric")])

## do it either way:
dat.interim2 <- data.table::copy(dat.interim[order(period_numeric),mean(V1), by=period_numeric])
summed_binom_dat[, mean(r/(r+n_r)), by=c("period_numeric")]

dat.interim2[,marginal.fit:=pnorm( -1.4212364+1.085814*period_numeric)]
dat.interim2[,conditional.fit:=pnorm( -1.99164+0.819914*period_numeric)]
setnames(dat.interim2, "V1", "empirical.avg")
setcolorder(dat.interim2, c("period_numeric", "conditional.fit","empirical.avg","marginal.fit" ))

dat.interim2[, glm.fit:=pnorm(-1.42435 + 1.10648 *period_numeric)]

print(dat.interim2, digits=2)
plot(x1, r/(r+n_r),col="grey", ylim=c(0,1),xlim=range(xdomain))
lines(xdomain, pnorm( -1.99164+0.819914*xdomain), ylim=c(0,1))
lines(xdomain, pnorm((-1.99164+0.819914*xdomain)*phi_x(xdomain)), col="red")
lines(xdomain, pnorm( -1.4212364+1.085814*xdomain), ylim=c(0,1),col="orange")
lines(xdomain, phi_x(xdomain), ylim=c(0,1), lty=2, col='cyan')
text(3.6,0.52, "(b0c+b1c*x)*phi(x)", col="red")
text(3.6,0.68, "(b0c+b1c*x)       ", col="black")
text(3.6,0.95, "(b0m+b1m*x)       ", col="orange")
text(3.6,0.19, "            phi(x)", col="cyan")



########################################$$$$#########
## 2022-02-16B                                     ##
## END two random parameters     MARGINAL COEFF    ##
#####################################################





########################################$$$$#########
## 2022-02-16                                      ##
## START two random parameters: MARGINAL COEFF     ##
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
sim_mrim_data <- function(n1, n2, J, a0, a1, v1=1,v2=2,rho=0.5, mrim="Probit-BIVARIATE_NORM", alpha=1.89, gamma=1.2, delta=1){

  if(mrim=="Logistic-BIVARIATE_NORM"){
    print("HERE")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) plogis(x)
  }
  if(mrim=="Probit-BIVARIATE_NORM"){
    print("Probit for the win")
    G <- function(n, v1, v2, rho){mvtnorm::rmvnorm(n=n, sigma=matrix(c(v1,rho*sqrt(v1*v2),rho*sqrt(v1*v2),v2),nrow=2))}
    H <- function(x) pnorm(x)
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
set.seed(1126)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(12,12, J=100, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat
## glmer -- only random intercept
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)




attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1

## PPN
(rand.int <-
    gnlrim::gnlrem(y=ybind,
                   mu = ~ pnorm(Intercept + period_numeric*b_p + rand1),
                   pmu = c(Intercept=-0.95, b_p=0.55),
                   pmix=c(var1=1),
                   p_uppb = c(  0,   2, 4.00),
                   p_lowb = c( -4,  -2, 0.05),
                   distribution="binomial",
                   nest=id,
                   random=c("rand1"),
                   mixture="normal-var",
                   ooo=TRUE,
                   compute_hessian = FALSE,
                   compute_kkt = FALSE,
                   trace=1,
                   method='nlminb',
    )
)




## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial(link="probit"), nAGQ = 1)





# (rand.int.nonzero.corr <-
#     gnlrem(y=ybind,
#            mu = ~ plogis(Intercept + period_numeric*b_p + rand1),
#            pmu = c(Intercept=-0.95, b_p=0.55),
#            pmix=c(var1=1),
#            p_uppb = c(  0,   2, 4.00),
#            p_lowb = c( -4,  -2, 0.05),
#            distribution="binomial",
#            nest=id,
#            random=c("rand1"),
#            mixture="normal-var",
#            ooo=TRUE,
#            compute_hessian = FALSE,
#            compute_kkt = FALSE,
#            trace=1,
#            method='nlminb'
#     )
# )

# (rand.int.rand.slopes.nonzero.corr <-
#     gnlrem(y=ybind,
#            mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
#            pmu = c(Intercept=-0.95, b_p=0.55),
#            pmix=c(var1=1, var2=1, corr12= 0.20),
#            p_uppb = c(  0,   2, 4.00, 4.00, 0.90),
#            p_lowb = c( -4,  -2, 0.05, 0.05,-0.90),
#            distribution="binomial",
#            nest=id,
#            random=c("rand1", "rand2"),
#            mixture="bivariate-normal-corr",
#            ooo=TRUE,
#            compute_hessian = FALSE,
#            compute_kkt = FALSE,
#            trace=1,
#            method='nlminb'
#     )
# )

(rand.int.rand.slopes.nonzero.corr.CUBA <-
    gnlrim::gnlrem(y=ybind,
           mu = ~ pnorm(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-1.95, b_p=1.55),
           pmix=c(var1=1, var2=1, corr12= 0.20),
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
           method='nlminb',
           int2dmethod="cuba"
    )
)
# Assemble the answers
#         Intercept     b_p    var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb  -1.76073 1.83619 1.05234 2.521599 0.1297366 195.0142     18     90    13        0   NA   NA

(rand.int.rand.slopes.nonzero.corr.CUBA.MARGINAL <-
    gnlrim::gnlrem(y=ybind,
           mu = ~ pnorm(
             (Intercept + period_numeric*b_p)*
               sqrt(1 + var1 + var2*period_numeric^2 + 2*corr12*sqrt(var1*var2)*period_numeric ) +
               rand1 + period_numeric*rand2
           ),
           pmu = c(Intercept=-1.95, b_p=1.55, var1=1, var2=1, corr12= 0.20),
           pmix=c(var1=1, var2=1, corr12= 0.20),
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
           method='nlminb',
           int2dmethod="cuba"
    )
)
#
# Assemble the answers
#        Intercept     b_p     var1     var2    corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb -1.229578 1.50625 1.029324 2.531012 0.1703418 194.9702     16     97    14        0   NA   NA

phi_x <- function(x, var1=1.05234, var2=2.521599, corr12=0.1297366){

  1/sqrt(1 + var1 + var2*x^2 + 2*corr12*sqrt(var1*var2)*x )

}
xvalueset <- c(0, sort(unique(period_numeric)), 1, 10)
cbind(
  xvalue = xvalueset,

  "phi(xvalue)" = phi_x(x=xvalueset),

  beta0_c = -1.76073,
  beta1_c = 1.83619,

  beta0_m =  -1.76073*  phi_x(x=xvalueset),
  beta1_m =  1.83619 *  phi_x(x=xvalueset)
)

xdomain <- seq(-2,4,0.05)
plot(xdomain, phi_x(xdomain), ylim=c(0,1))

plot(x1, r/(r+n_r),col="grey", ylim=c(0,1), xlim=range(xdomain))
lines(xdomain, pnorm( -1.76073+1.83619*xdomain), ylim=c(0,1))
lines(xdomain, pnorm(-1.229578+1.50625*xdomain), col="red")
lines(xdomain, phi_x(xdomain), ylim=c(0,1), lty=2, col='cyan')


########################################$$$$#########
## 2022-02-16                                      ##
## END two random parameters     MARGINAL COEFF    ##
#####################################################







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
set.seed(6)
binom_dat <-
  #  sim_mrim_data(800,400, J=100, a0 = -2, a1 = 1)
  #  sim_mrim_data(4000,2000, J=100, a0 = -2, a1 = 1)
  sim_mrim_data(8,8, J=100, a0 = -2, a1 = 1)
data.table::setDT(binom_dat)

summed_binom_dat <-
  binom_dat[, {j=list(r=sum(y), n_r=sum(y==0))}, by=c("id","x1")]
data.table::setkey(summed_binom_dat,id, x1)
summed_binom_dat
## glmer -- only random intercept
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial, nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id), summed_binom_dat, binomial, nAGQ = 1)




attach(summed_binom_dat)

ybind <- cbind(r,n_r)
period_numeric <- x1


(rand.int <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=1),
           p_uppb = c(  0,   2, 4.00),
           p_lowb = c( -4,  -2, 0.05),
           distribution="binomial",
           nest=id,
           random=c("rand1"),
           mixture="normal-var",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb',
    )
)




## glmer with correlation between random intercept and random slope
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial, nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial, nAGQ = 1)





# (rand.int.nonzero.corr <-
#     gnlrem(y=ybind,
#            mu = ~ plogis(Intercept + period_numeric*b_p + rand1),
#            pmu = c(Intercept=-0.95, b_p=0.55),
#            pmix=c(var1=1),
#            p_uppb = c(  0,   2, 4.00),
#            p_lowb = c( -4,  -2, 0.05),
#            distribution="binomial",
#            nest=id,
#            random=c("rand1"),
#            mixture="normal-var",
#            ooo=TRUE,
#            compute_hessian = FALSE,
#            compute_kkt = FALSE,
#            trace=1,
#            method='nlminb'
#     )
# )

## this is very slow compared to `CUBA` below
(rand.int.rand.slopes.nonzero.corr <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=1, var2=1, corr12= 0.20),
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

(rand.int.rand.slopes.nonzero.corr.CUBA <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=1, var2=1, corr12= 0.20),
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
           method='nlminb',
           int2dmethod="cuba"
    )
)

## wout corr
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id)+(0+x1 | id), summed_binom_dat, binomial, nAGQ = 0)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id)+(0+x1 | id), summed_binom_dat, binomial, nAGQ = 1)

(rand.int.rand.slopes.nonzero.corr.CUBA.wout.corr <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=.3, var2=.3, corr12= 0),
           p_uppb = c(  0,   2, 4.00, 4.00, 0),
           p_lowb = c( -4,  -2, 0.05, 0.05,-0),
           distribution="binomial",
           nest=id,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-corr",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb',
           int2dmethod="cuba"
    )
)

(take2 <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=.3, var2=.3, corr=0),
           p_uppb = c(  0,   2, 4.00, 4.00, 0.99),
           p_lowb = c( -4,  -2, 0.05, 0.05,-0.99),
           distribution="binomial",
           nest=id,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-corr",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb',
           int2dmethod="cuba"
    )
)

## lock down PQL estimated s1,s2,rho
(take3 <-
    gnlrem(y=ybind,
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-0.95, b_p=0.55),
           pmix=c(var1=0.952^2, var2=1.4878^2, corr=0.553),
           p_uppb = c(  0,   2, 0.952^2, 1.4878^2, 0.553),
           p_lowb = c( -4,  -2, 0.952^2, 1.4878^2, 0.553),
           distribution="binomial",
           nest=id,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-corr",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb',
           int2dmethod="cuba"
    )
)

## https://stats.stackexchange.com/questions/49039/reml-use-in-glmm
##cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial, nAGQ = 1
## ?MASS::glmmPQL
library(nlme) # will be loaded automatically if omitted
summary(MASS::glmmPQL(cbind(r, n_r) ~ x1, random = ~ x1 | id,
                      family = binomial, data = summed_binom_dat))
lme4::glmer(cbind(r, n_r) ~ x1 + (x1 | id), summed_binom_dat, binomial, nAGQ = 0)
##?GMMAT::glmmkin
##GMMAT::glmmkin(cbind(r, n_r) ~ x1, id="id", random.slope="x1")



g1 <- hglm::hglm(
  fixed = r/(r+n_r) ~ x1, weights=(r+n_r),
  random = ~ 1 | id,
  data = summed_binom_dat,
  family = binomial(link = logit),
  rand.family = gaussian(link = identity)
)
summary(g1)
yy <- r/(r+n_r)
g2 <- hglm::hglm2(
  meanmodel= yy ~ x1 + (x1 | id),
  weights=(r+n_r),
  data = summed_binom_dat,
  family = binomial(link = logit),
  rand.family = gaussian(link = identity),
  fix.disp=NULL
)
summary(g2)
lme4::glmer(cbind(r, n_r) ~ x1 + (1 | id)+(0+x1 | id), summed_binom_dat, binomial, nAGQ = 1)

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
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 1))

## using  nAGQ = 9  provides a better evaluation of the deviance
## Currently the internal calculations use the sum of deviance residuals,
## which is not directly comparable with the nAGQ=0 or nAGQ=1 result.
## 'verbose = 1' monitors iteratin a bit; (verbose = 2 does more):
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, verbose = 2, nAGQ =9))

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
         pmu = c(Intercept=2,b2=2, b3=2, b4=2),#pmu = c(Intercept=-1,b2=-1, b3=-1, b4=-1),
         pmix=c(var=1.00),
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
         abs.tol.nlminb = 1e-10,
         xf.tol.nlminb =  1e-10,
         x.tol.nlminb =   1e-10,
         rel.tol.nlminb = 1e-10,
         steps=10,
         points=5,
         eps=1e-04,

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
         abs.tol.nlminb = 1e-2,
         xf.tol.nlminb =  1e-2,
         x.tol.nlminb =   1e-2,
         rel.tol.nlminb = 1e-2,
         points = 2, #5 default
         steps = 2   # 10 is default
  )

 gm1.redux$value
gm1a.redux$value

## tips on covariance...
## https://stats.stackexchange.com/questions/414864/lme4glmer-get-the-covariance-matrix-of-the-fixed-and-random-effect-estimates


#
# (gm2.redux <-
#   gnlrem(y=cbind(incidence, size - incidence),
#          mu = ~ plogis(Intercept + period2*b2 + period3*b3 + period4*b4 + rand1 + rand2),
#          pmu = c(Intercept=-1,b2=-1, b3=-1, b4=-1),
#          pmix=c(var1=0.78, var2=0.40),
#          p_uppb = c(  2,   2, 2, 2,  1.0, 1.0),
#          p_lowb = c( -2,  -2,-2,-2,  0.1, 0.1),
#          distribution="binomial",
#          nest=herd,
#          random=c("rand1", "rand2"),
#          mixture="bivariate-normal-indep",
#          ooo=TRUE,
#          compute_hessian = FALSE,
#          compute_kkt = FALSE,
#          trace=1,
#          method='nlminb',
#   )
# )


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
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
           pmu = c(Intercept=-1, b_p=-1),
           pmix=c(var1=0.78, var2=0.40, corr12=0),
           p_uppb = c(  2,   2, 4.00, 4.00, 0),
           p_lowb = c( -2,  -2, 0.001, 0.00004, 0),
           distribution="binomial",
           nest=herd,
           random=c("rand1", "rand2"),
           mixture="bivariate-normal-corr",
           ooo=TRUE,
           compute_hessian = FALSE,
           compute_kkt = FALSE,
           trace=1,
           method='nlminb',
           int2dmethod="cuba"

    )
)
# [1] 5
# Intercept       b_p      var1      var2    corr12
# -1.00     -1.00      0.78      0.40      0.00
# [1] 100.0948
# fn is  fn
# Looking for method =  nlminb
# Function has  5  arguments
# par[ 1 ]:  -2   <? -1   <? 2     In Bounds
# par[ 2 ]:  -2   <? -1   <? 2     In Bounds
# par[ 3 ]:  0.001   <? 0.78   <? 4     In Bounds
# par[ 4 ]:  4e-05   <? 0.4   <? 4     In Bounds
# par[ 5 ]:  0   <? 0   <? 0     In Bounds
# Analytic gradient not made available.
# Analytic Hessian not made available.
# Scale check -- log parameter ratio= 0.39794   log bounds ratio= 0.0001085872
# Method:  nlminb
# 0:     100.09476: -1.00000 -1.00000 0.780000 0.400000  0.00000
# 1:     96.840530: -0.889769 -0.753339 0.781888 0.240449  0.00000
# 2:     93.411291: -0.921116 -0.611218 0.414012 4.00000e-05  0.00000
# 3:     93.361374: -0.915466 -0.592099 0.416240 0.0188371  0.00000
# 4:     93.269964: -0.909596 -0.586202 0.431892 4.00000e-05  0.00000
# 5:     93.268432: -0.902234 -0.569471 0.450896 0.00781833  0.00000
# 6:     93.245441: -0.904989 -0.573779 0.453331 4.00000e-05  0.00000
# 7:     93.244993: -0.906060 -0.560295 0.450888 4.00000e-05  0.00000
# 8:     93.240816: -0.909700 -0.566084 0.450205 4.00000e-05  0.00000
# 9:     93.238758: -0.922579 -0.565109 0.445501 4.00000e-05  0.00000
# 10:     93.238174: -0.922553 -0.562195 0.445604 4.00000e-05  0.00000
# 11:     93.237931: -0.925453 -0.562172 0.445910 4.00000e-05  0.00000
# 12:     93.237854: -0.927833 -0.560123 0.450825 4.00000e-05  0.00000
# 13:     93.237692: -0.931224 -0.558374 0.446415 4.00000e-05  0.00000
# 14:     93.237642: -0.932079 -0.560479 0.446628 4.00000e-05  0.00000
# 15:     93.237548: -0.933500 -0.559038 0.447681 4.00000e-05  0.00000
# 16:     93.237548: -0.933253 -0.559058 0.447741 4.00000e-05  0.00000
# 17:     93.237548: -0.933305 -0.559046 0.447725 4.00000e-05  0.00000
# 18:     93.237548: -0.933304 -0.559047 0.447725 4.00000e-05  0.00000
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -0.9333043 -0.5590465  0.4477254  0.0000400  0.0000000
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 93.23755
#
# $fevals
# function
# 25
#
# $gevals
# gradient
# 90
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
# 2779.862
#
# Assemble the answers
# Intercept        b_p      var1  var2 corr12    value fevals gevals niter convcode kkt1 kkt2
# nlminb -0.9333043 -0.5590465 0.4477254 4e-05      0 93.23755     25     90    18        0   NA   NA
# xtime
# nlminb 2779.862

(rand.int.rand.slopes.nonzero.corr <-
    gnlrem(y=cbind(incidence, size - incidence),
           mu = ~ plogis(Intercept + period_numeric*b_p + rand1 + period_numeric*rand2),
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
           method='nlminb',
           int2dmethod="cuba"
    )
)
# Intercept        b_p       var1       var2     corr12
# -0.9333689 -0.5589960  0.3686549  0.2531917 -0.1000000
# [1] 96.66891
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
# 0:     96.668908: -0.933369 -0.558996 0.368655 0.253192 -0.100000
# 1:     93.764320: -0.925143 -0.600387 0.401532 0.0705280 -0.139629
# 2:     93.663329: -0.794393 -0.508991 0.587035 0.0400000 -0.403583
# 3:     93.368641: -0.929231 -0.476850 0.561841 0.0400000 -0.516573
# 4:     93.135918: -1.05063 -0.488291 0.606646 0.0400000 -0.391118
# 5:     93.061465: -0.987166 -0.544726 0.447698 0.0400000 -0.379357
# 6:     93.003421: -1.00495 -0.518751 0.510988 0.0400000 -0.413143
# 7:     92.980272: -1.01102 -0.512586 0.519505 0.0400000 -0.443200
# 8:     92.843109: -1.02957 -0.491499 0.572051 0.0400000 -0.613687
# 9:     92.690070: -1.03509 -0.488232 0.643033 0.0400000 -0.769081
# 10:     92.417551: -0.995270 -0.521391 0.751831 0.0400000 -0.878620
# 11:     92.301027: -0.921752 -0.570864 0.813047 0.0400000 -0.816547
# 12:     92.116208: -0.843831 -0.621220 0.952859 0.106774 -0.818105
# 13:     92.108905: -0.829141 -0.637993 0.988450 0.0980816 -0.819566
# 14:     92.098351: -0.851070 -0.621007 0.955995 0.0958164 -0.823257
# 15:     92.080683: -0.868195 -0.607672 0.990312 0.0975466 -0.837036
# 16:     92.078781: -0.876284 -0.610849 0.987859 0.0932515 -0.837657
# 17:     92.075977: -0.872079 -0.611615 0.995128 0.0986009 -0.836853
# 18:     92.075055: -0.867671 -0.619043  1.01318 0.0987198 -0.835824
# 19:     92.072255: -0.897652 -0.600761 0.995213 0.0929377 -0.839901
# 20:     92.071222: -0.898548 -0.599909 0.998049 0.0981050 -0.839632
# 21:     92.069745: -0.895347 -0.602814  1.00165 0.0959737 -0.839253
# 22:     92.068704: -0.888963 -0.603302  1.01083 0.100449 -0.839293
# 23:     92.067546: -0.912729 -0.599588  1.00949 0.100118 -0.840250
# 24:     92.065667: -0.911289 -0.597023  1.01329 0.0960900 -0.839790
# 25:     92.065181: -0.911061 -0.596072  1.01740 0.100429 -0.841468
# 26:     92.063686: -0.908216 -0.600430  1.02047 0.0988163 -0.842168
# 27:     92.062576: -0.904611 -0.599756  1.02507 0.100845 -0.843054
# 28:     92.061936: -0.908856 -0.597960  1.02725 0.0975579 -0.842813
# 29:     92.060866: -0.906395 -0.599946  1.03095 0.100915 -0.844201
# 30:     92.060197: -0.905387 -0.601117  1.03622 0.0986433 -0.844509
# 31:     92.059329: -0.908543 -0.596741  1.03836 0.0999210 -0.844653
# 32:     92.058422: -0.906826 -0.600893  1.04156 0.101522 -0.846202
# 33:     92.057351: -0.907661 -0.597961  1.04666 0.101453 -0.846465
# 34:     92.056886: -0.915429 -0.596897  1.05498 0.0984724 -0.847833
# 35:     92.053833: -0.912204 -0.599004  1.06449 0.103531 -0.851084
# 36:     92.053320: -0.910584 -0.595329  1.07546 0.105772 -0.851029
# 37:     92.050592: -0.909178 -0.598999  1.08645 0.106157 -0.853318
# 38:     92.047510: -0.902631 -0.602922  1.13649 0.114179 -0.862802
# 39:     92.046690: -0.917413 -0.599439  1.15437 0.114772 -0.861209
# 40:     92.046615: -0.918075 -0.597644  1.15287 0.113208 -0.864805
# 41:     92.046506: -0.914415 -0.599314  1.15109 0.114132 -0.862900
# 42:     92.046501: -0.915416 -0.599059  1.15193 0.114079 -0.862989
# 43:     92.046501: -0.915315 -0.599069  1.15184 0.114072 -0.863020
# 44:     92.046501: -0.915313 -0.599070  1.15185 0.114073 -0.863017
#
# Post processing for method  nlminb
# Successful convergence!
#   Save results from method  nlminb
# $par
# Intercept        b_p       var1       var2     corr12
# -0.9153133 -0.5990704  1.1518453  0.1140729 -0.8630170
#
# $message
# [1] "relative convergence (4)"
#
# $convcode
# [1] 0
#
# $value
# [1] 92.0465
#
# $fevals
# function
# 53
#
# $gevals
# gradient
# 264
#
# $nitns
# [1] 44
#
# $kkt1
# [1] NA
#
# $kkt2
# [1] NA
#
# $xtimes
# user.self
# 328.956
#
# Assemble the answers
# Intercept        b_p     var1      var2    corr12   value fevals gevals niter convcode kkt1 kkt2
# nlminb -0.9153133 -0.5990704 1.151845 0.1140729 -0.863017 92.0465     53    264    44        0   NA   NA
# xtime
# nlminb 328.956



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
#set.seed(10) # converges, LRT_free < LRT_2, but fails confidence interval test
#set.seed(11) # converges, LRT_free < LRT_2, but coef too close to 0
#set.seed(12) #converges, LRT_free < LRT_2, but fails confidence interval test
#set.seed(13) #converges, LRT_free < LRT_2, and has reverse confidence interval test (ppn sig, free nonsig)
set.seed(14) # passes everything but CIs are both significant -- but PPN closer to 0.  Gives hope.
#set.seed(15) # FAIL
#set.seed(16) # FAIL
#set.seed(17) #FAIL
#set.seed(18) # FAIL
#set.seed(19) # FAIL
#set.seed(20) # FAIL
#set.seed(21) # FAIL
#set.seed(22) # this was close enough to make me want to tinker.
#set.seed(23) # FAIL
#set.seed(24) ## NICE Likelihood example >400 vs <400
dputted_dat<-
  dput(data.table::data.table(sim_mrim_data2(n1=25,n2=100,J=JJ,a0=a0.true,a1=a1.true, v = gam.g.true / sqrt(alp.true), mrim="SSS", alpha=1.69, gam=1/sqrt(1.69)))[,sum(y), by=c("id","group")])

robby <- data.table::data.table(id=dputted_dat$id, group=dputted_dat$group, y1=dputted_dat$V1, y0=JJ-dputted_dat$V1); data.table::setorder(robby, group, y1, y0); robby

robby[, mean(y1)  , by=group]
robby[, median(y1), by=group]
robby[, median(y1/(y1+y0)), by=group]
#
# tail(robby,25)
# robby[id==23 & group==1, y1:=45]
# robby[id==23 & group==1, y0:=15]
# robby[id== 4 & group==1, y1:=55]
# robby[id== 4 & group==1, y0:= 5]
# robby[id==10 & group==1, y1:=57]
# robby[id==10 & group==1, y0:= 3]

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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl-over-sqrt",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-scl",
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
         mixture="libstable4u-subgauss-phi",
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
         mixture="libstable4u-subgauss-phi",
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
         mixture="libstable4u-subgauss-phi",
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
         mixture="libstable4u-subgauss-phi",
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
         mixture="libstable4u-subgauss-phi",
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
## when using libstable4u-subgauss-phi
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

cond_PPN_phi_libstable4u <- gnlrim(y=y_cbind,
                                  mu=~stable_cdf2( a_cond+b_cond*dose+rand, c(alpha,0,1/sqrt(2),0)),
                                  pmu=c(a_cond=0,b_cond=0,alpha=2),
                                  pmix=c(alpha=2, phi=0.5),
                                  p_uppb = c(Inf ,  Inf, 2, 1-1e-5),
                                  p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                  distribution="binomial",
                                  nest=id,
                                  random="rand",
                                  mixture="libstable4u-subgauss-phi")



## next two lines the same? No, because phi assumes scale=1 (this is hardcoded).
## and to get same beta coeffs you must assume the 1/sqrt(2)
## so not surprising that betas are the same but phi is different.
cond_PPN_phi$coefficients
cond_PPN_phi_libstable4u$coefficients

cond_PPN_phi$maxlike
cond_PPN_phi_libstable4u$maxlike



marg_PPN_phi_libstable4u1 <- gnlrim(y=y_cbind,
                                   mu=~stable_cdf2( (a_marg+b_marg*dose)/phi+rand, c(alpha,0,1,0)),
                                   pmu=c(a_marg=0,b_marg=0,phi=0.5,alpha=2),
                                   pmix=c(alpha=2, phi=0.5),
                                   p_uppb = c(Inf ,  Inf, 2, 1-1e-5),
                                   p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                   distribution="binomial",
                                   nest=id,
                                   random="rand",
                                   mixture="libstable4u-subgauss-phi")



## next two lines the same? No, because phi assumes scale=1
## and to get same beta coeffs you must assume the 1/sqrt(2)
cond_PPN_phi_libstable4u1$coefficients
marg_PPN_phi_libstable4u1$coefficients
## no agreement.

##If alpha is 1?
cond_CCC_phi_libstable4u1 <- gnlrim(y=y_cbind,
                                   mu=~stable_cdf2( a_cond+b_cond*dose+rand, c(alpha,0,1,0)),
                                   pmu=c(a_cond=0,b_cond=0,alpha=1),
                                   pmix=c(alpha=1, phi=0.5),
                                   p_uppb = c(Inf ,  Inf, 1, 1-1e-5),
                                   p_lowb = c(-Inf, -Inf, 1, 0+1e-5),
                                   distribution="binomial",
                                   nest=id,
                                   random="rand",
                                   mixture="libstable4u-subgauss-phi")

marg_CCC_phi_libstable4u1 <- gnlrim(y=y_cbind,
                                   mu=~stable_cdf2( (a_marg+b_marg*dose)*((1^alpha)/phi)+rand, c(alpha,0,1,0)),
                                   pmu=c(a_marg=0,b_marg=0,alpha=1,phi=0.5),
                                   pmix=c(alpha=1, phi=0.5),
                                   p_uppb = c(Inf ,  Inf, 1, 1-1e-5),
                                   p_lowb = c(-Inf, -Inf, 1, 0+1e-5),
                                   distribution="binomial",
                                   nest=id,
                                   random="rand",
                                   mixture="libstable4u-subgauss-phi",
                                   compute_hessian=TRUE,
                                   compute_kkt=FALSE,
                                   trace=1)

cond_CCC_phi_libstable4u1$coefficients
marg_CCC_phi_libstable4u1$coefficients
marg_CCC_phi_libstable4u1$se

cond_CCC_phi_libstable4u1$maxlike
marg_CCC_phi_libstable4u1$maxlike



## no agreement.  I'm coming to the realization that `mixture=blah-blah-phi` should only be used when
## marginal coefficients and phi (explicitly) appear in the `mu` statement.  Furthermore, the
## parameters have to be listed correctly in `pmu` and `pmix`.  That is, the order of pmu must
## match the order the parameters appear in `mu`.  This means you might have to do a `1^alpha/phi` trick
## as opposed to `1/phi` because phi has to be last in pmix and pmu.  Note:  even the "trick" doesn't work for
## conditional models -- I put in alpha and phi explicitly:
# cond_CCC_phi_libstable4u1_trick <- gnlrim(y=y_cbind,
#                                          mu=~stable_cdf2( (a_cond+b_cond*dose)*1^alpha*1^phi +rand, c(alpha,0,1,0)),
#                                          pmu=c(a_cond=0,b_cond=0,alpha=1,phi=0.5),
#                                          pmix=c(alpha=1, phi=0.5),
#                                          p_uppb = c(Inf ,  Inf, 1, 1-1e-5),
#                                          p_lowb = c(-Inf, -Inf, 1, 0+1e-5),
#                                          distribution="binomial",
#                                          nest=id,
#                                          random="rand",
#                                          mixture="libstable4u-subgauss-phi")
# > cond_CCC_phi_libstable4u1_trick$maxlike
# [1] 9.989502 ## should be 10.03177
# > cond_CCC_phi_libstable4u1$maxlike
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

cond_PPN_scl_libstable4u <- gnlrim(y=y_cbind,
                                  mu=~stable_cdf2( a_cond+b_cond*dose+rand, c(alpha,0,1/sqrt(2),0)),
                                  pmu=c(a_cond=0,b_cond=0,alpha=2),
                                  pmix=c(alpha=2, scl=0.5),
                                  p_uppb = c(Inf ,  Inf, 2, Inf),
                                  p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                  distribution="binomial",
                                  nest=id,
                                  random="rand",
                                  mixture="libstable4u-subgauss-scl")



## next two lines the same? Yes!
c(cond_PPN_var$coefficients[1:2], 2, cond_PPN_var$coefficients[3],
  sqrt(cond_PPN_var$coefficients[3]/2))

c(cond_PPN_scl_libstable4u$coefficients[1:3],
  2*cond_PPN_scl_libstable4u$coefficients[4]^2 ,
  cond_PPN_scl_libstable4u$coefficients[4])

## can we do marginal with libstable4u?
## YES!, if we are sure to make sure to have alpha appear before scl in the `mu` statement
marg_PPN_scl_libstable4u <- gnlrim(y=y_cbind,
                                  mu=~stable_cdf2( (a_marg+b_marg*dose)*sqrt((1/sqrt(2))^alpha+scl^2)/(1/sqrt(2)) +rand, c(alpha,0,1/sqrt(2),0)),
                                  pmu=c(a_marg=0,b_marg=0,alpha=2, scl=0.5),
                                  pmix=c(alpha=2, scl=0.5),
                                  p_uppb = c(Inf ,  Inf, 2, Inf),
                                  p_lowb = c(-Inf, -Inf, 2, 0+1e-5),
                                  distribution="binomial",
                                  nest=id,
                                  random="rand",
                                  mixture="libstable4u-subgauss-scl")

## next two lines the same? Yes!
marg_PPN_scl_libstable4u$coefficients[1:2] * sqrt((1/sqrt(2))^2+marg_PPN_scl_libstable4u$coefficients[4]^2)/(1/sqrt(2))

c(cond_PPN_var$coefficients[1:2], 2, cond_PPN_var$coefficients[3],
  sqrt(cond_PPN_var$coefficients[3]/2))


cond_PPN_var$maxlike
marg_PPN_var$maxlike
marg_PPN_scl_libstable4u$maxlike
cond_PPN_scl_libstable4u$maxlike



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


fit_PPN_libstable4u_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf2(a+b*dose+rand, c(alpha, 0, 1/sqrt(2), 0)),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstable4u-subgauss-scl",
                                 ooo=TRUE,
                                 compute_hessian = FALSE,
                                 compute_kkt = FALSE,
                                 trace=1
)

fit_PPN_libstable4u_scl

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
## this section tests the libstable4u implementation of
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

#Now try reproducing PPN with **libstable4u** subgauss:

fit_PPN_libstable4u_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf2(a+b*dose+rand, c(alpha, 0, 1/sqrt(2), 0)),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstable4u-subgauss-scl"#,
                                 # ooo=TRUE,
                                 # compute_hessian = FALSE,
                                 # compute_kkt = FALSE,
                                 # trace=1
)

## should be 0s, not NAs.
attr(fit_PPN_libstable4u_scl, "details")[1,]$nhatend
bonk <- attr(fit_PPN_libstable4u_scl, "details")[1,]$nhatend
bonk_no_na <- matrix(bonk[!is.na(bonk)], nrow=3)
bonk_no_na_cov <- solve(bonk_no_na)
bonk_se <- sqrt(diag(bonk_no_na_cov))
## compared these two:
bonk_se
fit_PPN_scl$se



fit_PPN_libstable4u_scl$coefficients
fit_PPN_libstable4u_scl$se


# next two function are now in the package gnlrim
# ## what if we made a wrapper for libstable4u
# ## that would take out of bound pars and force them in bounds
# ## ***internally***
# stable_cdf2 <- function(x, pars, parameterization=0L, tol=1e-12){
#   libstable4u::stable_cdf(x,
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
#   libstable4u::stable_pdf(x,
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


fit2_PPN_libstable4u_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf2(a+b*dose+rand, c(alpha, 0, 1/sqrt(2), 0)),
                                 pmu=c(a=0,b=0, alpha=2),
                                 pmix=c(alpha=2, scl=1),
                                 p_uppb = c(Inf ,  Inf, 2, Inf),
                                 p_lowb = c(-Inf, -Inf, 2,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstable4u-subgauss-scl"#,
                                 #ooo=TRUE,
                                 #trace=1
)

## should be 0s, not NAs.
attr(fit2_PPN_libstable4u_scl, "details")[1,]$nhatend
bonk <- attr(fit2_PPN_libstable4u_scl, "details")[1,]$nhatend
bonk_no_na <- matrix(bonk[!is.na(bonk)], nrow=3)
bonk_no_na_cov <- solve(bonk_no_na)
bonk_se <- sqrt(diag(bonk_no_na_cov))
## compared these two:
bonk_se
fit_PPN_scl$se

fit2_PPN_libstable4u_scl$coefficients
fit2_PPN_libstable4u_scl$se

fit_PPN_scl$coefficients
fit_PPN_scl$se




## all the same?
fit_PPN_scl$coefficients
fit_PPN_stabledist_scl$coefficients
fit2_PPN_libstable4u_scl$coefficients

fit_PPN_scl$se
fit_PPN_stabledist_scl$se
fit2_PPN_libstable4u_scl$se


##############################################################
##############################################################
##############################################################
## this section tests the libstable4u implementation of
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

#Now try reproducing CCC with **libstable4u** subgauss:

fit_CCC_libstable4u_scl <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf(a+b*dose+rand, c(alpha, 0, 1, 0)),
                                 pmu=c(a=0,b=0, alpha=1),
                                 pmix=c(alpha=1, scl=1),
                                 p_uppb = c(Inf ,  Inf, 1, Inf),
                                 p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstable4u-subgauss-scl")
fit_CCC_libstable4u_scl$coefficients
fit_CCC_libstable4u_scl$se


## all the same?
fit_CCC_scl$coefficients
fit_CCC_stabledist_scl$coefficients
fit_CCC_libstable4u_scl$coefficients

fit_CCC_scl$se
fit_CCC_stabledist_scl$se
fit_CCC_libstable4u_scl$se



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

#Now try reproducing CCC with **libstable4u** subgauss:

fit_CCC_libstable4u_phi <- gnlrim(y=y_cbind,
                                 mu=~ stable_cdf(a+b*dose+rand, c(alpha, 0, 1, 0)),
                                 pmu=c(a=0,b=0, alpha=1),
                                 pmix=c(alpha=1, phi=0.5),
                                 p_uppb = c(Inf ,  Inf, 1,   1-1e-5),
                                 p_lowb = c(-Inf, -Inf, 1,   0+1e-5),
                                 distribution="binomial",
                                 nest=id,
                                 random="rand",
                                 mixture="libstable4u-subgauss-phi")
fit_CCC_libstable4u_phi$coefficients
fit_CCC_libstable4u_phi$se


## all the same?
fit_CCC_phi$coefficients
fit_CCC_stabledist_phi$coefficients
fit_CCC_libstable4u_phi$coefficients

fit_CCC_phi$se
fit_CCC_stabledist_phi$se
fit_CCC_libstable4u_phi$se



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
