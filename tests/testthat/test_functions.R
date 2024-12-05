# tests/testthat/test-functions.R
library(testthat)
library(future)
library(doFuture)
library(LPAI)  # Load your package

# Set future plan to sequential during testing
plan(sequential)
registerDoFuture()

test_that("train_protein works correctly", {

  age_range <- 50:85  # Age range for individuals
  n <- 45
  train.data <- data.frame(
    id = rep(1:15, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(age_range,size = 15,replace = F), times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(n, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(n, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )

  # Call train_protein
  result <- train_FPCA_protein(train.data, protein_indices = 3:4, cores = 2)

  # Check result structure
  expect_type(result, "list")
  expect_true("Q96PQ1" %in% names(result))
  expect_true("P25440" %in% names(result))

  # Check contents of each result
  expect_true("fpca_model" %in% names(result[["Q96PQ1"]]))
  expect_true("train_scores" %in% names(result[["P25440"]]))
})


test_that("get_cox_coef outputs correct coefficients", {
  # Simulated data

  age_range <- 50:85  # Age range for individuals
  n <- 45
  train.data <- data.frame(
    id = rep(1:15, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(age_range,size = 15,replace = F), times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(n, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(n, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )

  trained_models <- train_FPCA_protein(train.data, protein_indices = 3:4, cores = 2)

  surv_data <- data.frame(
    id = unique(train.data$id),
    time = rexp(15, rate = 0.1),
    status = sample(0:1, 15, replace = TRUE)
  )

  # Call get_cox_coef
  coefs <- get_cox_coef(trained_models, surv_data = surv_data, time_to_event = "time", event = "status")

  # Check coefficients
  expect_type(coefs, "double")
  expect_true("Q96PQ1_FPC1" %in% names(coefs))
  expect_true("P25440_FPC2" %in% names(coefs))

})


test_that("predict_protein_FPCscores returns correct scores", {

  age_range <- 50:85  # Age range for individuals
  n <- 45
  train.data <- data.frame(
    id = rep(1:15, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(age_range,size = 15,replace = F), times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(n, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(n, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )
  # Simulated data
  test.data <- data.frame(
    id = rep(20:29, each = 3),  # IDs, each individual has 3 observations
    age = rep(49:58, times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(10*3, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(10*3, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )
  # Mock trained models
  trained_models <- train_FPCA_protein(train.data, protein_indices = 3:4, cores = 2)
  # Call predict_protein_FPCscores
  predictions <- predict_protein_FPCscores(test.data, trained_models, protein_indices = 3:4, nFPCs = 2)

  # Check predictions structure
  expect_true("id" %in% colnames(predictions))
  expect_true("Q96PQ1_FPC1" %in% colnames(predictions))
  expect_true("Q96PQ1_FPC1" %in% colnames(predictions))
  expect_true("P25440_FPC1" %in% colnames(predictions))
  expect_true("P25440_FPC2" %in% colnames(predictions))
})



test_that("get_LPAI computes risk scores correctly", {

  age_range <- 50:85  # Age range for individuals
  n <- 45
  train.data <- data.frame(
    id = rep(1:15, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(age_range,size = 15,replace = F), times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(n, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(n, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )

  # Simulated data
  test.data <- data.frame(
    id = rep(20:29, each = 3),  # IDs, each individual has 3 observations
    age = rep(49:58, times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(10*3, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(10*3, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )
  # Mock trained models
  trained_models <- train_FPCA_protein(train.data, protein_indices = 3:4, cores = 2)
  # Call predict_protein_FPCscores
  test_fpc_scores <- predict_protein_FPCscores(test.data, trained_models, protein_indices = 3:4, nFPCs = 2)


  # Example Cox coefficients
  surv_data <- data.frame(
    id = unique(train.data$id),
    time = rexp(15, rate = 0.1),
    status = sample(0:1, 15, replace = TRUE)
  )

  # Call get_cox_coef
  coefs <- get_cox_coef(trained_models, surv_data = surv_data, time_to_event = "time", event = "status")


  # Call get_LPAI
  result <- get_LPAI(test_fpc_scores, coefs)

  # Check that the result has correct dimensions
  expect_equal(nrow(result), nrow(test_fpc_scores))
  expect_equal(ncol(result), 2)  # id and LPAI columns

  # Check that column names are correct
  expect_equal(colnames(result), c("id", "LPAI"))

  # Manually calculate LPAI for comparison
  expected_LPAI <- as.vector(as.matrix(test_fpc_scores[,-1]) %*% as.matrix(as.vector(coefs)))
  expect_equal(result$LPAI, expected_LPAI)

  # Check that no warnings or errors occur for valid input
  expect_warning(get_LPAI(test_fpc_scores, coefs), NA)
  expect_error(get_LPAI(test_fpc_scores, coefs), NA)
})


test_that("get_LPAI handles missing coefficients", {
  # Example FPC scores
  test_fpc_scores <- data.frame(
    id = 1:5,
    Protein1_FPC1 = c(0.5, 1.2, -0.8, 0.9, 0.4),
    Protein1_FPC2 = c(-0.2, 0.3, 1.5, -0.4, 0.6)
  )

  # Missing coefficients for Protein2
  cox_coefs <- c(
    Protein1_FPC1 = 0.3,
    Protein1_FPC2 = -0.1,
    Protein2_FPC1 = 0.5
  )

  # Expect error for missing coefficients
  expect_error(
    get_LPAI(test_fpc_scores, cox_coefs),
    "The following coefficients are missing in the FPC scores"
  )
})

test_that("get_LPAI handles extra FPC scores", {
  # Example FPC scores with an extra column
  test_fpc_scores <- data.frame(
    id = 1:5,
    Protein1_FPC1 = c(0.5, 1.2, -0.8, 0.9, 0.4),
    Protein1_FPC2 = c(-0.2, 0.3, 1.5, -0.4, 0.6),
    Extra_FPC = c(0.1, 0.2, 0.3, 0.4, 0.5)
  )

  # Example Cox coefficients
  cox_coefs <- c(
    Protein1_FPC1 = 0.3,
    Protein1_FPC2 = -0.1
  )

  # Expect a warning for extra FPC scores
  expect_warning(
    get_LPAI(test_fpc_scores, cox_coefs),
    "The following FPC scores are not used in the Cox model"
  )
})





test_that("predict_protein_FPC_scores works with Seqid argument", {
  # Load the id_mapping file from the installed package
  id_mapping_path <- system.file("extdata", "id_mapping.csv", package = "LPAI")
  expect_true(file.exists(id_mapping_path))  # Ensure the file exists


  age_range <- 50:85  # Age range for individuals
  n <- 45
  train.data <- data.frame(
    id = rep(1:15, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(age_range,size = 15,replace = F), times = 3),  # Ages for repeated measures
    P43320 = rnorm(n, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P04049 = rnorm(n, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )
  # Simulated data
  test.data <- data.frame(
    id = rep(20:29, each = 3),  # IDs, each individual has 3 observations
    age = rep(49:58, times = 3),  # Ages for repeated measures
    SeqId_10000_28 = rnorm(10*3, mean = 10, sd = 2),  # Simulated protein 1 measurements
    SeqId_10001_7 = rnorm(10*3, mean = 5, sd = 1)
  )
  # Mock trained models
  trained_models <- train_FPCA_protein(train.data, protein_indices = 3:4, cores = 2)
  # Call predict_protein_FPCscores
  # Test with Seqid
  result_with_seqid <- predict_protein_FPCscores(
    test_data = test.data,
    trained_models = trained_models,
    protein_indices = c(3, 4),
    nFPCs = 2,
    Seqid = TRUE
  )
  expect_true(all(c("P43320_FPC1", "P04049_FPC2") %in% colnames(result_with_seqid)))

  test.data.ex <- data.frame(
    id = rep(20:29, each = 3),  # IDs, each individual has 3 observations
    age = rep(49:58, times = 3),  # Ages for repeated measures
    SeqId_10000_28 = rnorm(10*3, mean = 10, sd = 2),  # Simulated protein 1 measurements
    SeqId_10001_7 = rnorm(10*3, mean = 5, sd = 1) , # Simulated protein 2 measurements
    SeqId_10615_18=rnorm(10*3, mean = 5, sd = 1)
  )

  expect_warning(
    result_with_seqid <- predict_protein_FPCscores(
      test_data = test.data.ex,
      trained_models = trained_models,
      protein_indices = c(3: 5),
      nFPCs = 2,
      Seqid = TRUE
    ),
    "Protein: SeqId_10615_18 not found in trained FPCA models. Skipping this protein."
  )

  })



