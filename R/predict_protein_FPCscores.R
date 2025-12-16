#' Predict Functional Principal Component (FPC) Scores for Test Data
#'
#' This function predicts FPC scores for test data using trained FPCA models.
#' It validates input data, handles duplicate observations, filters out age values outside the training range,
#' and returns the requested FPC scores for each protein. The function will print messages to the console during execution.
#'
#' @param test_data A data frame containing test data with columns `id`, `age`, and protein data.
#' @param trained_models A named list of trained FPCA models. Default is `ARIC_trained_models`.
#' @param protein_indices A numeric vector specifying the column indices of protein data in `test_data`.
#' @param nFPCs An integer specifying the number of FPC scores to extract. Default is 2.
#' @param Seqid Logical. If `TRUE`, rename user columns to UniProt IDs. Default is `FALSE`.
#'
#' @return A data frame with rows corresponding to test samples, containing `id` and FPC scores for each protein.
#' The FPC scores are named as `protein_name_FPC1`, `protein_name_FPC2`, etc.
#'
#' @importFrom dplyr rename_with recode group_by filter pull n select arrange mutate
#' @importFrom tidyr pivot_wider
#' @importFrom fdapace FPCA
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ARIC_trained_models)
#' data(example.test.data)
#'
#' # Predict FPC scores
#' test_fpc_scores <- predict_protein_FPCscores(
#'   test_data = example.test.data,
#'   trained_models = ARIC_trained_models,
#'   protein_indices = c(3, 4),
#'   nFPCs = 2, Seqid = T
#' )
#' print(test_fpc_scores)
#' }
#' @export
#'
#'
predict_protein_FPCscores <- function(test_data, trained_models = ARIC_trained_models, protein_indices,
                                      nFPCs = 2, Seqid = FALSE) {
  # Validate input: Ensure required columns are present
  if (!all(c("id", "age") %in% colnames(test_data))) {
    stop("The test data must include 'id' and 'age' columns.")
  }

  # Ensure age is an integer
  test_data$age <- floor(test_data$age)

  # Validate protein indices
  if (any(protein_indices <= 0 | protein_indices > ncol(test_data))) {
    stop("Invalid protein_indices: Ensure indices match the column positions in the test data.")
  }

  # Ensure required packages are available
  required_packages <- c("fdapace")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "The following required packages are missing: ",
      paste(missing_packages, collapse = ", "),
      ". Please install them and try again."
    )
  }


  # Check duplicate age
  dup <- test_data %>%
    group_by(id, age) %>%
    filter(n() > 1) %>%
    pull(id)
  if (length(dup) > 1) {
    stop(sprintf(
      "Duplicated IDs with duplicate age values found: %s. These observations have been removed.",
      paste(dup, collapse = ", ")
    ))
  }

  # Get the age range from the first trained model
  first_model <- trained_models[[1]][["fpca_model"]]
  age_range <- range(first_model$workGrid)
  # Filter out observations with age outside the range
  out_of_range <- test_data$age < age_range[1] | test_data$age > age_range[2]
  if (any(out_of_range)) {
    warning(sprintf(
      "%d observations were removed because their age values were outside the training range [%f, %f].",
      sum(out_of_range), age_range[1], age_range[2]
    ))
  }
  test_data <- test_data[!out_of_range, ]

  # Rename columns if Seqid is TRUE
  if (Seqid) {
    print("Renaming your data columns to UniProt IDs")
    id_mapping_path <- system.file("extdata", "id_mapping.csv", package = "LPAI")
    if (!file.exists(id_mapping_path)) {
      stop("The id_mapping.csv file is missing from the package.")
    }
    id_mapping <- read.csv(id_mapping_path, stringsAsFactors = FALSE)



    # Rename columns
    seqid_old <- colnames(test_data[protein_indices])

    uniprot <- ifelse(seqid_old %in% id_mapping$seqid,
      yes = id_mapping$uniprot[match(seqid_old, id_mapping$seqid)],
      no = seqid_old
    )
    seqid_dup <- seqid_old[duplicated(uniprot)]
    test_data <- test_data %>% select(-all_of(seqid_dup))
    # Subset id_mapping to include only columns in test_data
    relevant_mapping <- id_mapping[id_mapping$seqid %in% colnames(test_data), ]

    test_data_new <- test_data %>%
      dplyr::rename_with(
        .cols = relevant_mapping$seqid,
        .fn = ~ dplyr::recode(., !!!setNames(relevant_mapping$uniprot, relevant_mapping$seqid))
      ) %>%
      as.data.frame()
  } else {
    test_data_new <- test_data
  }


  # Initialize a list to store FPC scores for each protein
  score_list <- list()

  # Loop through each protein index
  for (protein_index in protein_indices) {
    # Get protein name
    protein_name <- colnames(test_data_new)[protein_index]
    # Check if a trained model exists for this protein
    if (!protein_name %in% names(trained_models)) {
      if (is.na(protein_name)) {
        next
      }
      message(sprintf("%s not found in pre-trained FPCA models. Skipping.", names(test_data)[protein_index]))
      next
    }

    print(sprintf("Calculating FPC scores for protein %s ", protein_name))

    # Prepare test data for the current protein
    data_wide <- prepare_data(test_data_new, protein_index)
    # Create Ly and Lt lists for prediction
    Ly_list <- lapply(2:ncol(data_wide), function(x) na.omit(data_wide[[x]]))
    Lt_list <- lapply(2:ncol(data_wide), function(x) data_wide$age[!is.na(data_wide[[x]])])

    # Get the trained FPCA model
    fpca_model <- trained_models[[protein_name]][["fpca_model"]]

    # Predict FPCA scores
    predicted_scores <- predict(fpca_model, newLy = Ly_list, newLt = Lt_list)

    # Extract the first two FPC scores and name the columns
    fpc_scores <- predicted_scores$scores[, 1:nFPCs, drop = FALSE]
    colnames(fpc_scores) <- paste0(protein_name, "_FPC", 1:nFPCs)

    # Add the scores to the list
    score_list[[protein_name]] <- fpc_scores
  }


  # Combine all scores into a single data frame
  id <- colnames(data_wide)[-1]
  final_scores <- as.data.frame(do.call(cbind, score_list))
  final_scores <- cbind(id, final_scores)

  # Return the combined data frame
  return(final_scores)
}
