#' Prepare Data for FPCA (Using Protein Index)
#'
#' This function processes data to convert it into a wide format for a specified protein column by its index.
#' The resulting wide-format data has `age` as rows and protein values for each `id` as columns.
#'
#' @param data A data frame containing training data with columns `id`, `age`, and protein data.
#' @param protein_index An integer specifying the index of the protein column to be converted.
#' @return A wide-format data frame where rows represent unique `age` values and columns represent protein values for each `id`.
#'
#' @importFrom dplyr select arrange n
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%

#' @export
prepare_data <- function(data, protein_index) {
  # Validate protein_index
  if (protein_index <= 0 || protein_index > ncol(data)) {
    stop("The specified protein_index is out of bounds")
  }

  # Select relevant columns using the index
  p=colnames(data)[protein_index]
  selected_columns <- c("id", "age", p)
  data_i <- data[,selected_columns] # Select id, age, and the specified protein column
  original_id_order <- as.character(unique(data$id))
  # Prepare the data in wide format
  data_wide <- data_i %>%
    tidyr::pivot_wider(names_from = id, values_from = all_of(p)) %>% # Convert to wide format
    dplyr::arrange(age) %>%
    dplyr::select(age, all_of(original_id_order)) # Arrange by age
  return(data_wide)
}
