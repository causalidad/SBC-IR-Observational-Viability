# IR Perturbative Sector Combined Likelihood with Pantheon+

## Components


This combined diagnostic includes:

- CMB Planck-lite spectral penalty
- DESI DR2 BAO-lite baseline
- f sigma8 growth validation
- Pantheon+ full covariance baseline

Because the current current infrared-sector implementation is primarily perturbative/spectral
rather than background-modifying, BAO and Pantheon+ are treated as consistency
baselines at this stage.

## Combined comparison

```text
         model  chi2_CMB  chi2_BAO  chi2_growth  chi2_Pantheon  chi2_total  k_params  N_eff   AIC_proxy   BIC_proxy
LCDM_reference   0.00000  32.99668     1.039661    1458.160871 1492.197212         0   1650 1492.197212 1492.197212
SBC_F3_current   2.12626  32.99668     1.039661    1458.160871 1494.323472         3   1650 1500.323472 1516.549064
```


## Differences

Delta chi2 = 2.126260
Delta AIC  = 8.126260
Delta BIC  = 24.351852


## Interpretation

LCDM remains the statistical reference model and is not outperformed by the infrared-sector extension
in the current proxy comparison.

However, the infrared-sector extension introduces only a controlled additional CMB spectral penalty
while preserving:

- BAO consistency,
- Pantheon+ full-covariance background consistency,
- f sigma8 growth consistency,
- high-ell TT/EE/TE acoustic structure.

This supports SBC/F3 as a viable weak IR-safe phenomenological deformation,
not yet as a statistically preferred replacement for LCDM.
