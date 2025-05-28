# BiTSLS: Bidirectional Two-Stage Least Squares

A simple R package for estimating bidirectional causal effects using proxy variables.

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("Jiaqi-Min/BiTSLS")
```

## Usage

The `Bi_TSLS()` function estimates bidirectional causal effects between X and Y:

```r
library(BiTSLS)

# Prepare your data with required variables
data <- data.frame(
  X = ...,  # Treatment variable
  Y = ...,  # Outcome variable
  Z = ...,  # Negative control exposure
  W = ...,  # Negative control outcome
  # Additional covariates (At least one covariate)
)

# Run the estimation
result <- Bi_TSLS(data)

# View results
print(result)  # Effect of X on Y and Y on X
```

## Requirements

Your data must contain:
- `X`: Treatment/exposure variable (numeric)
- `Y`: Outcome variable (numeric)
- `Z`: Negative control exposure (numeric)
- `W`: Negative control outcome (numeric)
- Additional covariates are optional (At least one covariate)

## Sensitivity Analysis

You can test sensitivity to violations of the proxy structural conditions:

```r
# With sensitivity parameters
result <- Bi_TSLS(data, R_w = 0.1, R_z = -0.1)
```

## License

MIT License
