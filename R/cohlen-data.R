#' Cohlen et al data
#'
#' A dataset containing the fertility crossover trial.
#' Acquired from Supplementary materials of
#' Reference
#' Makubate, B., and Senn, S. (2010),Planning and Analysis of Cross-over Trials in Infertility, Statistics in Medicine, 29, 3203-3210.
#'
#' @format A data frame with 320 rows and 10 variables for 74 unique patients:
#' \describe{
#'   \item{group}{group id}
#'   \item{patient}{patient id}
#'   \item{period2}{binary variable for period 2}
#'   \item{period3}{binary variable for period 3}
#'   \item{period4}{binary variable for period 4}
#'   \item{period5}{binary variable for period 5}
#'   \item{period6}{binary variable for period 6}
#'   \item{treatment}{binary}
#'   \item{response}{outcome -- }
#'   \item{period}{integer variable for period}
#' }
#' @source \url{http://senns.uk/InfCros/InfCROSInt.html}
#' @examples
#' \dontrun{
#' #Proceed to fit various models
#' fit1 <- glmer(response~(1|patient),
#'               family=binomial,data=cohlen)#null model
#' fit2 <- glmer(response~treatment+(1|patient),family=binomial,data=cohlen)#treatment only
#' fit3 <- glmer(response~period2 + period3 + period4 + period5 +
#'                       period6 +(1|patient),family=binomial,data=cohlen)#period only, as a factor
#' #full model,period as a factor
#' fit4 <- glmer(response~treatment + period2 + period3 + period4 + period5 + period6+(1|patient),
#'                family=binomial,data=cohlen)
#' #period only, Period having a linear effect:
#' fit5 <- glmer(response~period+(1|patient),family=binomial,data=cohlen)
#' #full model,period having a linear effect:
#' fit6 <- glmer(response~treatment+period+(1|patient),family=binomial,data=cohlen)
#' summary(fit2)#Summary of model with treatment only
#' summary(fit4)#Summary of model with treatment and term
#carry out analysis of deviance to check effect of adding treatment to model with term
#' anova(fit3,fit4)
#' }
"cohlen"
