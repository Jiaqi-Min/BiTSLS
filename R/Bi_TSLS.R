#' Bidirectional Two-Stage Least Squares Estimation
#'
#' @description
#' Performs bidirectional two-stage least squares (Bi-TSLS) estimation to identify
#' bidirectional causal effects between X and Y using proxy variables Z and W.
#' Requires at least one covariate in addition to X, Y, Z, and W.
#'
#' @param data A data frame containing at least these variables:
#'   \itemize{
#'     \item X: Treatment/exposure variable
#'     \item Y: Outcome variable
#'     \item Z: Negative control exposure (proxy variable)
#'     \item W: Negative control outcome (proxy variable)
#'     \item At least one additional covariate (required)
#'   }
#' @param R_w Sensitivity parameter for W (default = 0)
#' @param R_z Sensitivity parameter for Z (default = 0)
#'
#' @return A list containing:
#'   \itemize{
#'     \item beta_XY: Estimated effect of X on Y
#'     \item beta_YX: Estimated effect of Y on X
#'   }
#'
#' @examples
#' # Create example data with covariate
#' n <- 1000
#' data <- data.frame(
#'   X = rnorm(n),
#'   Y = rnorm(n),
#'   Z = rnorm(n),
#'   W = rnorm(n),
#'   V = rnorm(n)  # At least one covariate is required
#' )
#'
#' # Run BiTSLS
#' result <- Bi_TSLS(data)
#' print(result)
#' @export
#' @importFrom stats lm coef as.formula
Bi_TSLS <- function(data, R_w = 0, R_z = 0) {

  # Check for required variables
  required_vars <- c("X", "Y", "Z", "W")
  if (!all(required_vars %in% names(data))) {
    stop("Data must contain variables: X, Y, Z, W")
  }

  # Check if variables are numeric
  for (var in required_vars) {
    if (!is.numeric(data[[var]])) {
      stop(paste("Variable", var, "must be numeric"))
    }
  }

  # Identify covariates
  covariate_names <- setdiff(names(data), required_vars)

  # Check if at least one covariate is present
  if (length(covariate_names) == 0) {
    stop("At least one covariate is required for Bi-TSLS estimation")
  }

  # Rename covariates to V1, V2, etc.
  V_names <- paste0("V", seq_along(covariate_names))
  names(data)[match(covariate_names, names(data))] <- V_names
  V_terms <- paste(V_names, collapse = " + ")

  # Create interaction terms
  ZV_interactions <- paste("Z:", V_names, collapse = " + ")
  WV_interactions <- paste("W:", V_names, collapse = " + ")

  # First stage regressions with interactions
  formula_W <- as.formula(paste("W ~ Z +", V_terms, "+", ZV_interactions))
  formula_Z <- as.formula(paste("Z ~ W +", V_terms, "+", WV_interactions))

  fit_W <- lm(formula_W, data = data)
  fit_Z <- lm(formula_Z, data = data)

  # Second stage regressions
  formula_X_Z <- as.formula(paste("X ~ Z + fit_W$fitted.values +", V_terms))
  formula_Y_Z <- as.formula(paste("Y ~ Z + fit_W$fitted.values +", V_terms))
  formula_X_W <- as.formula(paste("X ~ fit_Z$fitted.values + W +", V_terms))
  formula_Y_W <- as.formula(paste("Y ~ fit_Z$fitted.values + W +", V_terms))

  lm_X_Z <- lm(formula_X_Z, data = data)
  lm_Y_Z <- lm(formula_Y_Z, data = data)
  lm_X_W <- lm(formula_X_W, data = data)
  lm_Y_W <- lm(formula_Y_W, data = data)

  k1 <- coef(lm_X_W)["W"] / coef(lm_Y_W)["W"]
  k2 <- coef(lm_Y_Z)["Z"] / coef(lm_X_Z)["Z"]

  beta_xy_Bi_TSLS <- ((k2 * (1 + k1 * R_z - R_w * R_z)) - R_z) / (1 - k1 * k2 * R_w * R_z)
  beta_yx_Bi_TSLS <- ((k1 * (1 + k2 * R_w - R_w * R_z)) - R_w) / (1 - k1 * k2 * R_w * R_z)
  return(c(beta_xy_Bi_TSLS, beta_yx_Bi_TSLS))
}
