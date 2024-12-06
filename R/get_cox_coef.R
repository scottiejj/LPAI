#' Perform Cox Regression Using FPCA Train Scores
#'
#' This function performs Cox regression using the concatenated functional principal component (FPC) scores
#' from FPCA results across all models. It uses elastic net regularization to fit the model, with an option
#' to specify a custom lambda or use cross-validation to select the optimal lambda.
#'
#' @param trained_models A named list of trained FPCA models, where each element contains `fpca_model` and associated train scores.
#' @param id A string specifying the column name in `surv_data` containing IDs (default: `"id"`).
#' @param surv_data A data frame with survival data, including columns for time-to-event and event status.
#' @param time_to_event A string specifying the column name in `surv_data` for time-to-event data.
#' @param event A string specifying the column name in `surv_data` for event status (binary: 1 = event, 0 = censored).
#' @param alpha A numeric value (default: 0.5) for elastic net mixing (0 = Ridge, 1 = LASSO).
#' @param nFPCs A numeric value specifying the number of principal components to use (default: 2).
#' @param lambda A numeric value specifying a fixed lambda for the elastic net. If `NULL`, lambda is selected using 10fold cross-validation.
#' @return A named numeric vector of coefficients for the selected principal components, with descriptive names.
#'
#' @importFrom dplyr rename inner_join arrange select mutate
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival Surv
#' @examples
#' \dontrun{
#'
#' data("ARIC_trained_models")  # Example FPCA trained models
#' surv_data <- data.frame(
#'   id = 1:100,
#'   time = rexp(100, rate = 0.1),
#'   status = sample(0:1, 100, replace = TRUE)
#' )

#'coefs <- get_cox_coef(
#'   trained_models = ARIC_trained_models,
#'   id = "id",
#'   surv_data = surv_data,
#'   time_to_event = "time",
#'   event = "status",
#'   alpha = 0.5,
#'   nFPCs = 2
#' )
#' print(coefs)
#' }
#'
#'
#'

#' @export
get_cox_coef <- function(trained_models, id = "id", surv_data, time_to_event, event, alpha = 0.5,nFPCs=2, lambda=NULL) {
  # Validate inputs
  if (!is.data.frame(surv_data)) {
    stop("The `surv_data` parameter must be a data frame.")
  }
  if (!all(c(time_to_event, event) %in% colnames(surv_data))) {
    stop("The `time_to_event` and `event` columns must exist in `surv_data`.")
  }

  # Extract IDs from `trained_models`
  trained_ids <- trained_models[["id"]]

  # Concatenate the first two `train_scores` (FPC1 and FPC2) across all models
  concatenated_scores <- data.frame(
    id = trained_ids,
    do.call(cbind, lapply(names(trained_models)[-length(trained_models)], function(protein_name) {
      # Access fpca_model for each protein
      fpca_model <- trained_models[[protein_name]]$fpca_model
      return(fpca_model$xiEst[, 1:nFPCs, drop = FALSE])
    } )
    )
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
  if (is.null(lambda)) {
    cv_el <- glmnet::cv.glmnet(concatenated_scores, y, family = "cox", alpha = alpha)
    el_model <- glmnet::glmnet(concatenated_scores, y, family = "cox", alpha = alpha,
                               lambda = cv_el$lambda.min)
  } else {
    el_model <- glmnet::glmnet(concatenated_scores, y, family = "cox", alpha = alpha,
                               lambda = lambda)
  }

  # Extract coefficients
  coefs <- as.vector(coef(el_model))
  protein_names <- names(trained_models)[-length(trained_models)]
  coef_names <- unlist(lapply(protein_names, function(protein) {
    paste0(protein, "_FPC", 1:nFPCs)
  }))
  names(coefs) <- coef_names

  # Return the coefficients
  return(coefs)
}
