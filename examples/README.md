# Examples for clamp_binary()

This directory contains example scripts to demonstrate how `clamp_binary()` works with confounded binary treatments.

## Files in this Directory

### 1. `test_clamp_binary_simple.R`
A streamlined test script that's easy to run and understand.

**To run:**
```R
# From R console, navigate to the examples directory
setwd("/path/to/clamp/examples")
source("test_clamp_binary_simple.R")
```

### 2. `test_clamp_binary_example.R`
A comprehensive test with detailed output and diagnostics.

**To run:**
```R
# From R console, navigate to the examples directory
setwd("/path/to/clamp/examples")
source("test_clamp_binary_example.R")
```

### 3. `check_dependencies.R`
Check if all required R packages and source files are available.

**To run:**
```R
setwd("/path/to/clamp/examples")
source("check_dependencies.R")
```

## What the Test Does

### Data Generation

1. **Confounder**: Generates a 1-dimensional confounder `Z ~ N(0, 1)`

2. **Binary Treatments**: Creates 10 binary treatments `X1, ..., X10` where each treatment is confounded by Z:
   ```
   P(Xj = 1 | Z) = logistic(αj + βj * Z)
   ```
   - Each treatment has its own intercept αj ~ Uniform(-1, 1)
   - Each treatment has its own slope βj ~ Uniform(0.5, 2)

3. **Response**: Generates outcome `y` using only 2 true causal treatments (X3 and X7) plus the confounder:
   ```
   y = 2.0 + 1.5*Z + 3.0*X3 - 2.5*X7 + ε
   ```
   where ε ~ N(0, 1)

4. **Inverse Probability Weights**: Computes IPW for each treatment by:
   - Fitting logistic regression: `P(Xj = 1 | Z)`
   - Computing weights: `W[i,j] = 1/ps[i]` if `X[i,j]=1`, else `1/(1-ps[i])`
   - Truncating extreme weights at 1st and 99th percentiles

### What to Expect

If `clamp_binary()` works correctly, it should:

✓ **Detect X3 and X7** as causal variables (high PIP > 0.5)
✓ **Estimate effects** close to true values: X3 ≈ +3.0, X7 ≈ -2.5
✓ **Exclude X1, X2, X4, X5, X6, X8, X9, X10** (low PIP < 0.5)
✓ **Converge** within 100 iterations

### Example Output

```
=== RESULTS ===

Posterior Inclusion Probabilities:
  Variable   PIP  True
1       X3 0.982  TRUE
2       X7 0.967  TRUE
3       X5 0.034 FALSE
4       X2 0.015 FALSE
...

Posterior Mean Effects:
  Variable Estimated  True
1       X3     2.953   3.0
2       X7    -2.487  -2.5
3       X1     0.021   0.0
...

✓ SUCCESS: Both true causal variables detected!
```

## Understanding the Test

This test simulates a realistic causal inference scenario:

- **Confounding**: All treatments are associated with Z, but only X3 and X7 actually cause changes in y
- **Selection Bias**: Naive regression of y on X would show spurious associations for all treatments because they're all correlated with Z
- **IPW Correction**: By using inverse probability weights, `clamp_binary()` should correctly identify only the truly causal treatments

This is exactly the problem that methods like CLAMP are designed to solve!

## Troubleshooting

If the test fails to detect the true variables:

1. **Check convergence**: Look at `fit$converged` and `fit$niter`
2. **Increase iterations**: Set `max_iter = 200`
3. **Increase bootstrap samples**: Set `nboots = 100`
4. **Check weights**: Make sure IPW matrix W doesn't have extreme values
5. **Try different seed**: The random seed affects data generation
