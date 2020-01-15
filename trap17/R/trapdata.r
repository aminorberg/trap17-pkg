#' TRAP17 DATA SET
#'
#' The processed trap17 data set of five viruses on 320 sentinel plants in Åland
#' @format A list of data frames, character vectors and formulas.
#' \describe{
#'   \item{Y}{Virus occurrences.}
#'   \item{Y_pooled}{Virus occurrences, pooled over the two sampling time points.}
#'   \item{X}{Explanatory variables.}
#'   \item{X_pooled}{Explanatory variables, pooled over the two sampling time points.}
#'   \item{PI}{Sampling design.}
#'   \item{PI_pooled}{Sampling design after pooling over the time points, i.e. running number of sampling units.}
#'   \item{yvars}{Names of viruses.}
#'   \item{xvars}{Names of explanatory variables.}
#'   \item{pivars}{sampling design variables.}
#'   \item{pool_}{Formulas for pooling.}
#'   \item{gropuping_pooled}{The unit over which variables were grouped.}
#'   \item{pivars_pooled}{The unit over which variables were pooled.}
#' }
"trapdata"

