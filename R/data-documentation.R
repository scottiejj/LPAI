#' ARIC Coefficients
#'
#' This dataset contains the coefficients or each protein's first 2 FPC trained from 4684 proteins from the ARIC study
#'
#' @format A named vector with coefficients for each protein's first 2 FPC
#' @details
#' The coefficients were extracted from the trained models and can be used for directly calculating LPAI.
#'
"aric_coef"



#' ARIC Trained Models
#'
#' This data contains the trained FPCA model for each protein based on the ARIC study data.
#'
#' @format A list containing trained FPCA models and training FPC scores for each protein:
#' @details
#' The models were trained using the FPCA method on ARIC study data.

"ARIC_trained_models"


#' Example training data
#'
#' This data contains example training data with columns id, age and proteins
#'
#' @format a dataframe

"example.train.data"


#' Example Test data
#'
#' This data contains example test data with columns id, age and proteins
#'
#' @format a dataframe

"example.test.data"

