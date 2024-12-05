
# LPAI
LPAI provides tools for calculating the Longitudinal Proteomic Aging Index using Functional Principal Component Analysis (FPCA). 
It allows users to either directly calculate LPAI using pretrained weights from 185 FPCs of 181 unique proteins or train their own models on custom datasets.

## Installation

You can install the latest version of **LPAI** from GitHub using the `devtools` package:

```
# Install devtools if you haven't already
install.packages("devtools")

# Install LPAI from GitHub
devtools::install_github("scottiejj/LPAI")

```

## Detailed Steps

### Example: Calculating LPAI Using Pretrained Weights

This example demonstrates how to calculate LPAI directly using our pretrained weights from 4952 proteins measured in the ARIC study.
The protein data is measured using the SomaScan platform, expressed in relative fluorescence units (RFU), and is NOT log-transformed.

- Input data should Be in wide format.

- Contain columns with exact names: id, age.

- Have protein columns named using either SeqId or UniProt ID. SeqId: Format SeqId_XXXXX_XX. UniProt ID: Format QXXXXX.

``` 
library(LPAI)

# Load pre-trained FPCA models and Cox coefficients
data(ARIC_trained_models) # Load pre-trained FPCA models
data(aric_coef) # Load pre-trained Cox coefficients

# Create test data
test.data <- data.frame(
    id = rep(1:20, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(40:80, 20), times = 3),  # Ages for repeated measures
    SeqId_10000_28 = rnorm(20 * 3, mean = 10, sd = 2),  # Simulated protein 1 measurements
    SeqId_10001_7 = rnorm(20 * 3, mean = 7, sd = 1),   # Simulated protein 2 measurements
    SeqId_10615_18 = rnorm(20 * 3, mean = 8, sd = 1)   # Simulated protein 3 measurements
)

# Calculate FPC scores
test.fpc.scores <- predict_protein_FPCscores(
    test_data = test.data,
    trained_models = ARIC_trained_models,
    protein_indices = c(3:5),  #columns of protein
    Seqid = TRUE # protein columns are named in Seqid format
)

# Calculate LPAI
lpai <- get_LPAI(test_fpc_scores=test.fpc.scores, cox_coefs = aric_coef)
print(lpai)


```

### Example: training a new LPAI model
``` r
library(LPAI)
#Simulate training data 
age_range <- 50:85 
  n <- 45
  train.data <- data.frame(
    id = rep(1:15, each = 3),  # IDs, each individual has 3 observations
    age = rep(sample(age_range,size = 15,replace = F), times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(n, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(n, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )

# Simulate test data
  test.data <- data.frame(
    id = rep(20:29, each = 3),  # IDs, each individual has 3 observations
    age = rep(49:58, times = 3),  # Ages for repeated measures
    Q96PQ1 = rnorm(10*3, mean = 10, sd = 2),  # Simulated protein 1 measurements
    P25440 = rnorm(10*3, mean = 5, sd = 1)  # Simulated protein 2 measurements
  )
# Simulate Time to Event Data

  surv_data <- data.frame(
    id = unique(train.data$id),
    time = rexp(15, rate = 0.1),
    status = sample(0:1, 15, replace = TRUE)
  )


# Training
trained_models <- train_FPCA_protein(train.data, protein_indices = 3:4, cores = 2) #specify how many cores to use for parallel computation
coefs <- get_cox_coef(trained_models, surv_data = surv_data,id="id",time_to_event = "time", event = "status",alpha = 0.5,nFPCs=2) #compute coefficients assigned to each FPC

# Calculate FPC score in a test set 
test.fpc.scores <- predict_protein_FPCscores(test.data, trained_models, protein_indices = 3:4, nFPCs = 2)
# Calculate LPAI
lpai.test <- get_LPAI(test_fpc_scores=test.fpc.scores, cox_coefs = coefs)
print(lpai.test)

```



              

                          
                          
                          

                          
                          