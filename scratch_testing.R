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
