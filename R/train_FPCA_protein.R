#' Train FPCA Models for Proteins
#'
#' Trains Functional Principal Component Analysis (FPCA) models for specified protein columns in the dataset using parallel processing.
#'
#' @param data A data frame containing `id`, `age`, and protein columns.
#' @param protein_indices Numeric vector of column indices for the protein data in `data`.
#' @param cores Number of cores for parallel computation (default: 1).
#' @return A named list where each element corresponds to a protein, containing:
#' \describe{
#'   \item{fpca_model}{The trained FPCA model.}
#'   \item{train_scores}{Data frame of FPCA scores for each sample.}
#' }
#'
#' @importFrom foreach foreach
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multisession
#' @importFrom dplyr group_by filter pull n
#' @importFrom tidyr pivot_wider
#' @importFrom fdapace FPCA
#' @importFrom magrittr %>%
#' @importFrom doRNG %dorng%
#'
#' @examples
#' \dontrun{

#'data(example.train.data)
#'
#' result <- train_FPCA_protein(example.train.data, protein_indices = 3:4, cores = 2)
#' print(result[[1]])
#'}
#'
#' @export
#'
train_FPCA_protein <- function(data, protein_indices, cores = 1) {
  # Check for required packages
  required_packages <- c("fdapace", "foreach", "doFuture", "doRNG","LPAI")
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
  data$age<-floor(data$age) #Ensure age is an integer

  if (any(protein_indices <= 0 | protein_indices > ncol(data))) {
    stop("Invalid 'protein_indices': Ensure they are positive integers and within the range of column indices in the data.")
  }

  protein_names <- colnames(data)[protein_indices]

  #Check duplicate age

  dup<-data %>%
    group_by(id,age) %>%
    filter(n() > 1) %>%
    pull(id)
  if(length(dup)>1){
    stop("Duplicate age values for same ID found in the data. Please ensure each age value is unique.")
  }

  # Set up parallel processing with doFuture
  doFuture::registerDoFuture()
  future::plan(future::multisession, workers = cores)

  # Parallel loop to train FPCA models
  results <-foreach(
    protein_index = protein_indices,
    .packages = c("fdapace", "dplyr", "tidyr", "magrittr","LPAI"),
    .export = c("prepare_data"),.multicombine=TRUE
  ) %dorng% {
    # Get column name for the current index
    protein_name <- colnames(data)[protein_index]

    # Prepare data for this protein
    data_wide <- prepare_data(data, protein_index)

    # Create Ly and Lt lists for FPCA
    Ly_list <- lapply(2:ncol(data_wide), function(x) na.omit(data_wide[[x]]))
    Lt_list <- lapply(2:ncol(data_wide), function(x) data_wide$age[!is.na(data_wide[[x]])])
    sample_ids <- colnames(data_wide)[-1]

    # Train FPCA model
    fpca_model <- FPCA(Ly = Ly_list, Lt = Lt_list, optns = list(dataType="Sparse"))

    # FPCA scores
    scores <- data.frame(id = sample_ids, fpca_model$xiEst)

    # Remove unnecessary objects to free memory
    rm(Ly_list, Lt_list, data_wide)
    gc()

  # Return results for this protein
    list(
      fpca_model = fpca_model,
      train_scores = scores
    )
  }
  # Set names for results using protein name
  print((results))
  names(results)<-protein_names

  results <- lapply(results, function(x) {
    attr(x, "rng") <- NULL
    return(x)
  })

  return(results)
}
