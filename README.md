# Uni-cox
Univariable analysis for cox regression, with optional covariates

### Uni-variable cox regression

For genetic analysis, we would like to know which gene has effects on patients prognosis the most. However, the lineline package only handle with multi-variable analysis, which mixs up all genes. Therefore, implementing a package for signle gene regression is neccersary. 
With covariants:
$$
H(t) = H_0(t) \times \text{exp}(b_1x_1 + b_2x_2 + ...b_kx_k)
$$

### Usage
Run in the terminal
Here we uses rossi as template data(lifeline.dataset)
