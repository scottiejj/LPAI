#' Train Functional Principal Component Analysis (FPCA) Models for Proteins
#'
#' This function trains FPCA models for specified protein columns in the dataset using parallel processing.
#' It validates input data, handles duplicates, and ensures compatibility of column indices.
#'
#' @param data A data frame containing columns `id`, `age`, and protein measurements.
#' @param protein_indices Numeric vector indicating the column indices of the protein data in `data`.
#' @param cores Integer specifying the number of cores for parallel computation (default: 1).
#' @return A named list where each element corresponds to a protein, containing:
#' \describe{
#'   \item{fpca_model}{The trained FPCA model.}
#'   \item{id}{A vector of unique sample IDs extracted from the input data.}
#' }
#'
#' @details
#' This function processes the input data and trains FPCA models for each specified protein column in parallel.
#' Each FPCA model captures the functional structure of protein data with respect to age.
#'
#' The function ensures that age values are rounded to integers before training the models.
#'
#' Parallel processing is implemented using the `foreach` and `doParallel` packages.
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr group_by filter pull n
#' @importFrom tidyr pivot_wider
#' @importFrom fdapace FPCA
#' @importFrom magrittr %>%

#' @examples
#' \dontrun{
#' data(example.train.data)
#' result <- train_FPCA_protein(example.train.data, protein_indices = 3:4, cores = 2)
#' print(result[[1]])
#' }
#'
#' @export
train_FPCA_protein <- function(data, protein_indices, cores = 1) {
  # Check for required packages
  required_packages <- c("fdapace", "foreach", "doParallel", "dplyr", "tidyr", "magrittr", "LPAI")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "The following required packages are missing: ",
      paste(missing_packages, collapse = ", "),
      ". Please install them and try again."
    )
  }

  # Validate input
  required_columns <- c("id", "age")
  missing_columns <- required_columns[!required_columns %in% colnames(data)]
  if (length(missing_columns) > 0) {
    stop(sprintf("The following required columns are missing from the data: %s. Please include them.",
                 paste(missing_columns, collapse = ", ")))
  }
  data$age <- floor(data$age) # Ensure age is an integer

  if (any(protein_indices <= 0 | protein_indices > ncol(data))) {
    stop("Invalid 'protein_indices': Ensure they are positive integers and within the range of column indices in the data.")
  }

  protein_names <- colnames(data)[protein_indices]

  # Check for duplicate age values per ID
  dup <- data %>%
    dplyr::group_by(id, age) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::pull(id)
  if (length(dup) > 0) {
    stop("Duplicate age values for the same ID found in the data. Please ensure each (id, age) pair is unique.")
  }

  # Set up parallel processing with doParallel
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- NULL  # Will run sequentially if cores = 1
  }

  log_file <- "protein_train_log.txt"
  sample_ids <- as.character(unique(data$id))

  results <- tryCatch({
    results <- foreach::foreach(
      protein_index = protein_indices,
      .packages = c("fdapace", "dplyr", "tidyr", "magrittr", "LPAI"),
      .export = c("prepare_data")
    ) %dopar% {
      # Get column name for the current index
      protein_name <- colnames(data)[protein_index]
      cat(sprintf("Processing protein %s \n", protein_name), file = log_file, append = TRUE)

      # Prepare data for this protein
      data_wide <- prepare_data(data, protein_index)

      # Create Ly and Lt lists for FPCA
      Ly_list <- lapply(2:ncol(data_wide), function(x) na.omit(data_wide[[x]]))
      Lt_list <- lapply(2:ncol(data_wide), function(x) data_wide$age[!is.na(data_wide[[x]])])

      # Train FPCA model
      fpca_model <- fdapace::FPCA(Ly = Ly_list, Lt = Lt_list, optns = list(dataType = "Sparse"))

      # Cleanup
      rm(Ly_list, Lt_list, data_wide)
      gc()

      # Return results for this protein
      list(
        fpca_model = fpca_model
      )
    }
    results
  },
  interrupt = function(e) {
    message("Training interrupted.")
    return(NULL)
  },
  finally = {
    # Ensure cluster is stopped in case of error or normal completion
    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }
  })



  # Set names for results using protein names
  names(results) <- protein_names
  results[["id"]] <- sample_ids

  cat(sprintf("Finished training."), file = log_file, append = TRUE)
  return(results)
}
