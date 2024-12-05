#' Calculate LPAI Using FPC Scores and Cox Coefficients
#'
#' This function applies Cox regression coefficients to compute the LPAI (risk score) for each individual.
#'
#' @param test_fpc_scores A data frame containing FPC scores for the test data. The first column must be `id`.
#' @param cox_coefs A named numeric vector of coefficients from the trained Cox model.
#' @return A data frame with `id` and the computed risk score (`LPAI`).

#' @examples
#' \dontrun{
#' #Load example data
#'data(ARIC_trained_models)
#'data(example.test.data)
#'data(aric_coef)

#' #Predict FPC scores
#'test_fpc_scores <- predict_protein_FPCscores(
#' test_data = example.test.data,
#' trained_models = ARIC_trained_models,
#' protein_indices = c(3, 4),
#' nFPCs = 2,Seqid = T
#' )

#' lpai<-get_LPAI(test_fpc_scores, aric_coef)
#' print(lpai)
#'
#'}
#'
#' @export
get_LPAI <- function(test_fpc_scores, cox_coefs) {
  # Validate inputs
  fpc_scores=test_fpc_scores[,-1]
  if (!is.data.frame(fpc_scores)) {
    stop("The `fpc_scores` parameter must be a data frame.")
  }
  if (!is.numeric(cox_coefs)) {
    stop("The `cox_coefs` parameter must be a numeric vector.")
  }


  # Filter non-zero coefficients
 non_zero_coefs <- cox_coefs[cox_coefs != 0]
 common_names <- intersect(names(non_zero_coefs), colnames(fpc_scores))

  # Check if all non-zero coefficients are in the FPC scores
  missing_coefs <- setdiff(names(non_zero_coefs), colnames(fpc_scores))
  if (length(missing_coefs) > 0) {
    message(sprintf(
      "The following coefficients are not present in the FPC scores and will be ignored: %s.",
      paste(missing_coefs, collapse = ", ")
    ))
  }


  # Subset FPC scores and coefficients to only common names
  fpc_scores_subset <- fpc_scores[, common_names, drop = FALSE]
  coefs_subset <- non_zero_coefs[common_names]

  # Calculate the linear predictor using the common names
  linear_predictor <- as.matrix(fpc_scores_subset) %*% as.matrix(as.vector(coefs_subset))

  # Return as a data frame with a risk score column
  risk_scores <- data.frame(
    id = test_fpc_scores$id,
    LPAI = as.vector(linear_predictor)
  )

  return(risk_scores)
}

