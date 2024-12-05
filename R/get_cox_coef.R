#' Perform Cox Regression Using FPCA Train Scores
#'
#' This function concatenates the first `nFPCs` principal component scores (e.g., FPC1, FPC2)
#' from FPCA results across all models and fits a Cox regression model using elastic net regularization.
#' Lambda is selected using 10-fold cross-validation.
#'
#' @param trained_models A list of trained FPCA models where each element contains `train_scores`.
#' @param id A string specifying the column name in `surv_data` containing IDs (default: `"id"`).
#' @param surv_data A data frame with survival data, including columns for time-to-event and event status.
#' @param time_to_event A string specifying the column name in `surv_data` for time-to-event data.
#' @param event A string specifying the column name in `surv_data` for event status (binary: 1 = event, 0 = censored).
#' @param alpha A numeric value (default: 0.5) for elastic net mixing (0 = Ridge, 1 = LASSO).
#' @param nFPCs A numeric value specifying the number of principal components to use (default: 2).
#' @return A named numeric vector of coefficients for the selected principal components.
#'
#' @importFrom dplyr rename inner_join arrange select mutate
#' @import glmnet
#' @importFrom survival Surv

#' @export
get_cox_coef <- function(trained_models, id = "id", surv_data, time_to_event, event, alpha = 0.5,nFPCs=2) {
  # Validate inputs
  if (!is.data.frame(surv_data)) {
    stop("The `surv_data` parameter must be a data frame.")
  }
  if (!all(c(time_to_event, event) %in% colnames(surv_data))) {
    stop("The `time_to_event` and `event` columns must exist in `surv_data`.")
  }

  # Extract IDs from `trained_models`
  trained_ids <- trained_models[[1]]$train_scores[, 1]

  # Concatenate the first two `train_scores` (FPC1 and FPC2) across all models
  concatenated_scores <- data.frame(
    id = trained_ids,
    do.call(cbind, lapply(trained_models, function(x) x$train_scores[, 2:(nFPCs+2-1)]))
  )

  # Rename ID column in `surv_data` to match `id`
  surv_data <- surv_data %>% dplyr::rename(id = !!id)%>%mutate(id=as.character(id))

  # Filter `surv_data` using inner_join to match `train_scores` IDs
  filtered_surv_data <- dplyr::inner_join(
    surv_data,
    data.frame(id = trained_ids),
    by = "id"
  )

  # Check if any IDs were dropped
  if (nrow(filtered_surv_data) < nrow(surv_data)) {
    dropped_ids <- setdiff(surv_data$id, trained_ids)
    warning(sprintf(
      "%d IDs in `surv_data` were not found in `trained_models` and have been removed: %s.",
      length(dropped_ids), paste(dropped_ids, collapse = ", ")
    ))
  }

  # Filter and align `concatenated_scores` with `surv_data` IDs
  concatenated_scores <- dplyr::inner_join(
    filtered_surv_data[, c("id"), drop = FALSE],
    concatenated_scores,
    by = "id"
  ) %>% dplyr::arrange(id) %>% dplyr::select(-id)

  concatenated_scores <- as.matrix(concatenated_scores)

  # Prepare Cox regression inputs
  time_vector <- filtered_surv_data[[time_to_event]]
  event_vector <- filtered_surv_data[[event]]
  y <- survival::Surv(time = time_vector, event = event_vector)

  # Train ELastic Net Cox regression
  cv_el <- glmnet::cv.glmnet(concatenated_scores, y, family = "cox", alpha = alpha)
  el_model <- glmnet::glmnet(concatenated_scores, y, family = "cox", alpha = alpha,
                                      lambda = cv_el$lambda.min)

  # Extract coefficients
  coefs <- as.vector(coef(el_model))
  protein_names <- names(trained_models)
  coef_names <- unlist(lapply(protein_names, function(protein) {
    paste0(protein, "_FPC", 1:nFPCs)
  }))
  names(coefs) <- coef_names

  # Return the coefficients
  return(coefs)
}
