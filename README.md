
# Longitudinal Proteomic Aging Index(LPAI)

LPAI is an R package for computing the Longitudinal Proteomic Aging Index using Functional Principal Component Analysis (FPCA). The package supports two workflows:

- Using the included pretrained model based on 181 proteins selected from a 4,684-protein panel to compute LPAI for new test data.
- Training your own model on longitudinal proteomics data: fit FPCA, then an elastic-net Cox model to derive a custom LPAI.

The pretrained model was developed from proteins measured in the ARIC study (Atherosclerosis Risk in Communities) using the SomaScan platform; values are expressed in relative fluorescence units (RFU) and log-transformed.

The package includes example data and ID mappings to help you get started quickly.

**Note**: Each individual must have at least two longitudinal protein measurements (two or more ages/visits) to compute LPAI.

## Installation

``` r
# Install from GitHub
install.packages("devtools")
devtools::install_github("scottiejj/LPAI")

# Recommended dependencies used internally
install.packages(c(
  "fdapace", "glmnet", "survival", "dplyr", "tidyr", "magrittr",
  "foreach", "doParallel"
))
```

## Data Requirements

- **Columns required**: **`id`, `age`** (**EXACT** naming), plus one column per protein (numeric).
- Long format input: each row is one observation (id at a given age), with protein columns alongside `id` and `age`.
- Ages are rounded down to integers internally; ensure $(id, age)$ pairs are unique.
- Per-person minimum: at least two measurements across age (≥2 visits) to enable FPCA scoring and LPAI computation.
- Protein names: either UniProt IDs (e.g., `Q96PQ1`) or SomaScan `SeqId_*` names. If your test data uses SeqIds, set `Seqid = TRUE` in `predict_protein_FPCscores()` to map to UniProt via `inst/extdata/id_mapping.csv`.

## Quick Start (Pretrained Models)

``` r
library(LPAI)

# Pretrained FPCA models and Cox coefficients
data(ARIC_trained_models)
data(example.test.data)
data(aric_coef)

# Predict FPC scores for a subset of proteins (columns 3 and 4 here)
test_fpc_scores <- predict_protein_FPCscores(
  test_data = example.test.data,
  trained_models = ARIC_trained_models,
  protein_indices = c(3, 4),
  Seqid = TRUE  # convert SeqId names in your data to UniProt
)

# Compute LPAI (linear predictor from Cox coefficients)
lpai <- get_LPAI(test_fpc_scores = test_fpc_scores, cox_coefs = aric_coef)
lpai

```

Tips:
- If your test data already uses UniProt IDs that match the pretrained models, set `Seqid = FALSE`.
- If you select proteins not present in the pretrained models, they are skipped with a message.
- Observations with ages outside the training age range are dropped with a warning.

## Train Your Own LPAI

``` r
library(LPAI)

# Simulate longitudinal training data (long format)
set.seed(1)
n <- 45
train.data <- data.frame(
  id = rep(1:15, each = 3),
  age = rep(sample(50:85, size = 15), times = 3),
  Q96PQ1 = rnorm(n, mean = 10, sd = 2),
  P25440 = rnorm(n, mean = 5, sd = 1)
)

# Simulate a test set in the same format
test.data <- data.frame(
  id = rep(20:29, each = 3),
  age = rep(49:58, times = 3),
  Q96PQ1 = rnorm(30, mean = 10, sd = 2),
  P25440 = rnorm(30, mean = 5, sd = 1)
)

# Survival data aligned to training IDs
surv_data <- data.frame(
  id = unique(train.data$id),
  time = rexp(15, rate = 0.1),
  status = sample(0:1, 15, replace = TRUE)
)

# 1) Train FPCA models (logs to protein_train_log.txt)
# Parallel computing: increase 'cores' to use multiple CPU cores via doParallel.

trained_models <- train_FPCA_protein(
  data = train.data,
  protein_indices = 3:4,
  cores = 2 # parallel computing; more cores -> faster training
)

# 2) Fit elastic-net Cox and get coefficients
coefs <- get_cox_coef(
  trained_models = trained_models,
  surv_data = surv_data,
  id = "id",
  time_to_event = "time",
  event = "status",
  alpha = 0.5,
  nFPCs = 2
)

# 3) Predict FPC scores in test data (no renaming needed here)
test.fpc.scores <- predict_protein_FPCscores(
  test_data = test.data,
  trained_models = trained_models,
  protein_indices = 3:4,
  nFPCs = 2,
  Seqid = FALSE
)

# 4) Compute LPAI in test set
lpai.test <- get_LPAI(test_fpc_scores = test.fpc.scores, cox_coefs = coefs)
lpai.test
```

## Interpreting LPAI

- LPAI is a linear predictor from a Cox model: $\text{LPAI} = \sum_i x_i \ \beta_i$, where $x_i$ are the first two FPC scores.
- Higher LPAI indicates higher estimated risk (hazard).
- Units are arbitrary; the scale depends on cohort. For comparison:
  - Compare individuals in the same cohort by differences in LPAI; a difference of $\Delta$ implies $\exp(\Delta)$ fold change in hazard.
  - You can center or standardize LPAI (e.g., z-score) to improve interpretability across cohorts.

## Troubleshooting

- Missing packages: install `fdapace`, `glmnet`, `survival`, `dplyr`, `tidyr`, `magrittr`, `foreach`, `doParallel`.
- Duplicate `(id, age)` rows: ensure one row per id per age; duplicates cause an error.
- Ages outside training range: such observations are dropped with a warning during prediction.
- Protein not in trained models: the protein is skipped; check names or use `Seqid = TRUE` to map.
- Invalid `protein_indices`: ensure indices refer to protein columns (after `id`, `age`).
- Missing `id_mapping.csv`: reinstall the package or check `inst/extdata/id_mapping.csv` is available.


## Citation

Rao, Z., S. Wang, A. Li, et al. 2025. "A Novel Longitudinal Proteomic Aging Index Predicts Mortality, Multimorbidity, and Frailty in Older Adults." Aging Cell e70317. https://doi.org/10.1111/acel.70317.

## License

See `LICENSE` in this repository.
                          










