# Notes on confounder correction for _cis_-eQTL calling

We identified and removed effects from non-genetic confounding factors
using [Half-sibling
regression](http://www.pnas.org/content/113/27/7391.long) framework.
Similar ideas were explored in RNA-seq data, e.g.,
[`RUV-g`](https://www.nature.com/articles/nbt.2931) where they use
control genes to correct confounders.

For each LD block:

1. Include histone peaks / genes whose putative QTLs could be located within this LD block.  Let's call them `Y1` and genotype matrix in this LD block `G`.

2. Find control peak signals `Y0` that qualify the following criteria:

    a. `Y0` are located out of this LD block

    b. `Y0` are significantly correlated with `Y1` (potentially due to confounding factor as most _cis_-regulatory programs do not span far beyond the LD).  We retain top signals with correlation `p < 1e-4`

    c. There is no apparent correlation between `Y0` and `G`.  We remove some of `Y0` with `p < 1e-4`.

3. Fit a regression model : `Y1 ~ (Y0 - Cδ)θ` where `C` are known demographic confounders such age and gender.

4. Take residuals : `Y1 <- Y1 - (Y0 - Cδ)θ`.

5. Call QTLs : `Y1 ~ G`

