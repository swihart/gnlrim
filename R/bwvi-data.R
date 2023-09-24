#' Black White Visual Impairment data set
#'
#' Visual impairment dataset acquired from
#' supplementary materials of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4203426/
#'
#' References
#'
#'  - Tielsch JM, Sommer A, Katz J, Quigley H, Ezrine S. Socioeconomic status and visual impairment among urban Americans. Archives of ophthalmology. 1991;109:637–641.
#'
#'  - Liang KY, Zeger SL. Longitudinal data analysis using generalized linear models. Biometrika. 1986;73:13–22.
#'
#'  - Swihart, B. J., Caffo, B. S. and Crainiceanu, C. M. (2014) A unifying framework for marginalised random-476
#' intercept models of correlated binary outcomes. International Statistical Review, 82, 275–295
#'
#' @format A data frame with 10398 rows and 4 variables for 5199 unique patients:
#' \describe{
#'   \item{id}{subject id}
#'   \item{black}{1 if subject was black; 0 if white}
#'   \item{variable}{eye1 for left eye, eye2 for right eye}
#'   \item{value}{1 if the eye is impaired, 0 if healthy}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4203426/}
#' @examples
#' \donttest{
#' ## Example 3.5 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4203426/
#' attach(bwVI)
#' head(bwVI)
#' o.value <- value ## lesson learned!  got an error when I used "value"
#' (LLB <-
#'     gnlrim::gnlrim(y=cbind(o.value, 1-o.value),
#'                    mu = ~ plogis(Intercept + black*b_p + rand1),
#'                    pmu = c(Intercept = -4.88, b_p=0.134),
#'                    pmix=c(var=2.03),
#'                    p_uppb = c(  50,   9, 20.00),
#'                    p_lowb = c( -50,  -9,  0.05),
#'                    distribution="binomial",
#'                    nest=id,
#'                    random=c("rand1"),
#'                    mixture="logit-bridge-var",
#'                    ooo=TRUE,
#'                    compute_hessian = FALSE,
#'                    compute_kkt = FALSE,
#'                    trace=1,
#'                    method='nlminb',
#'     )
#' )
#' # > head(bwVI)
#' # id black variable value
#' # 1  1     1     eye1     1
#' # 2  1     1     eye2     1
#' # 3  2     1     eye1     1
#' # 4  2     1     eye2     1
#' # 5  3     1     eye1     1
#' # 6  3     1     eye2     1
#' # > o.value <- value ## lesson learned!  got an error when I used "value"
#' # > (LLB <-
#' #      +     gnlrim::gnlrim(y=cbind(o.value, 1-o.value),
#' #                           +                    mu = ~ plogis(Intercept + black*b_p + rand1),
#' #                           +                    pmu = c(Intercept = -4.88, b_p=0.134),
#' #                           +                    pmix=c(var=2.03),
#' #                           +                    p_uppb = c(  50,   9, 20.00),
#' #                           +                    p_lowb = c( -50,  -9,  0.05),
#' #                           +                    distribution="binomial",
#' #                           +                    nest=id,
#' #                           +                    random=c("rand1"),
#' #                           +                    mixture="logit-bridge-var",
#' #                           +                    ooo=TRUE,
#' #                           +                    compute_hessian = FALSE,
#' #                           +                    compute_kkt = FALSE,
#' #                           +                    trace=1,
#' #                           +                    method='nlminb',
#' #                           +     )
#' #    + )
#' # [1] 3
#' # Intercept       b_p       var
#' # -4.880     0.134     2.030
#' # [1] 3119.63
#' # fn is  fn
#' # Looking for method =  nlminb
#' # Function has  3  arguments
#' # par[ 1 ]:  -50   <? -4.88   <? 50     In Bounds
#' # par[ 2 ]:  -9   <? 0.134   <? 9     In Bounds
#' # par[ 3 ]:  0.05   <? 2.03   <? 20     In Bounds
#' # Analytic gradient not made available.
#' # Analytic Hessian not made available.
#' # Scale check -- log parameter ratio= 1.561315   log bounds ratio= 0.7447275
#' # Method:  nlminb
#' # 0:     3119.6299: -4.88000 0.134000  2.03000
#' # 1:     2791.2694: -4.06317 0.507976  2.46924
#' # 2:     2743.1744: -3.16525 0.142932  2.71519
#' # 3:     2729.2153: -3.48196 -0.400397  3.49268
#' # 4:     2710.6152: -3.83306 0.175929  4.23063
#' # 5:     2709.7961: -3.77051 0.184101  4.28020
#' # 6:     2709.2022: -3.80107 0.130731  4.33171
#' # 7:     2708.0638: -3.76409 0.123961  4.48770
#' # 8:     2702.5219: -3.97565 0.0273767  5.52687
#' # 9:     2700.4295: -4.24106 0.170311  6.49712
#' # 10:     2700.0148: -4.28328 0.109940  6.91920
#' # 11:     2699.9724: -4.30180 0.104920  7.08190
#' # 12:     2699.9702: -4.30653 0.106829  7.11938
#' # 13:     2699.9701: -4.30680 0.107369  7.12115
#' # 14:     2699.9701: -4.30677 0.107430  7.12092
#' # Post processing for method  nlminb
#' # Successful convergence!
#' #   Save results from method  nlminb
#' # $par
#' # Intercept        b_p        var
#' # -4.3067653  0.1074302  7.1209235
#' #
#' # $message
#' # [1] "relative convergence (4)"
#' #
#' # $convcode
#' # [1] 0
#' #
#' # $value
#' # [1] 2699.97
#' #
#' # $fevals
#' # function
#' # 18
#' #
#' # $gevals
#' # gradient
#' # 51
#' #
#' # $nitns
#' # [1] 14
#' #
#' # $kkt1
#' # [1] NA
#' #
#' # $kkt2
#' # [1] NA
#' #
#' # $xtimes
#' # user.self
#' # 752.63
#' #
#' # Assemble the answers
#' # Intercept       b_p      var   value fevals gevals niter convcode kkt1
#' # nlminb -4.306765 0.1074302 7.120923 2699.97     18     51    14        0   NA
#' # kkt2  xtime
#' # nlminb   NA 752.63
#'
#' }
"bwVI"
