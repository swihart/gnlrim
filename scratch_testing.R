dose <- c(c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8), 3+c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8))
y_01 <- rep(as.numeric(c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
       1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
       20.972267, 17.178012) > 4),2)
id <- rep(1:8, each=5)

y_cbind = cbind(y_01, 1-y_01)

fit_PPN_var <- gnlrim(y=y_cbind,
                      mu=~pnorm(a+b*dose+rand),
                      pmu=c(0,0),
                      pmix=1,
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="normal-var")

fit_PPN_phi <- gnlrim(y=y_cbind,
                      mu=~pnorm(a+b*dose+rand),
                      pmu=c(0,0),
                      pmix=0.5,
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
                      pmu=c(0,0),
                      pmix=1,
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="logit-bridge-var")

fit_LLB_phi <- gnlrim(y=y_cbind,
                      mu=~plogis(a+b*dose+rand),
                      pmu=fit_LLB_var$coefficients[1:2],
                      pmix=0.25,
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
                      pmu=c(0,0),
                      pmix=1,
                      p_uppb = c(Inf ,  Inf, Inf),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="Cauchy-scl")

fit_CCC_phi <- gnlrim(y=y_cbind,
                      mu=~pcauchy(a+b*dose+rand),
                      pmu=fit_CCC_scl$coefficients[1:2],
                      pmix=0.25,
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


fit_nlminb <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(8.7,0.25),
                            pshape=3.44, pmix=2.3, trace=1)
fit_nlminb$coefficients

fit_nlminb_bounds <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(8.7,0.25),
                     pshape=3.44, pmix=8, trace=1, p_lowb=c(8,.2,3,8), p_uppb=c(9,.3,4,10))
fit_nlminb_bounds$coefficients



fit_nlm <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(8.7,0.25),
                  pshape=3.44, pmix=2.3, method="nlm")

fit_both <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(8.7,0.25),
                  pshape=3.44, pmix=2.3, method=c("nlminb","nlm"), ooo=TRUE)
fit_both
## with `ooo=TRUE` hessian extraction (rows or by name):
attr(fit_both, "details")[1,]$nhatend
attr(fit_both, "details")[2,]$nhatend
attr(fit_both, "details")["nlminb",]$nhatend
attr(fit_both, "details")["nlm",]$nhatend


fit_both[1, 1:attr(fit_both, "npar")]



fit_both <- gnlrim(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(8.7,0.25),
                   pshape=3.44, pmix=2.3, method=c("nlminb","nlm"))
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
