# Uni-cox
Univariable analysis for cox regression, with optional covariates

### Uni-variable cox regression

For genetic analysis, we would like to know which gene has effects on patients prognosis the most. However, the lineline package only handle with multi-variable analysis, which mixs up all genes. Therefore, implementing a package for signle gene regression is neccersary. 
With covariants:
$$
H(t) = H_0(t) \times \text{exp}(b_1x_1 + b_2x_2 + ...b_kx_k)
$$

### Usage
1. Run in the terminal
Here we uses rossi as template data(lifeline.dataset)
2. Run with streamlit
```
streamlit run test_cox_app.py
```
See following example:
<img width="958" alt="image" src="https://user-images.githubusercontent.com/97223450/206990768-2ff36d4f-ce1d-4ae5-aed6-6a3c8e170ad8.png">

