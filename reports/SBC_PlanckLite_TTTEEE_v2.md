# IR Perturbative Sector Planck-lite TTTEEE

## Method

This test uses the compressed Planck-lite TTTEEE dataset:

- W_TT, W_EE, W_TE window matrices
- compressed vectors y_TT, y_EE, y_TE
- full compressed covariance C_obs

The infrared-sector benchmark is:

epsilon = 0.060
k0      = 0.015
n       = 1.0

## Results

         model  chi2_TTTEEE_v2  Delta_chi2_vs_LCDM  N_bins
LCDM_reference    1.404989e-17        0.000000e+00     613
   SBC_F3_best    1.403362e-17       -1.627901e-20     613

## Channel diagnostics

channel  N_bins  RMS_LCDM_vs_obs  RMS_SBC_vs_obs  RMS_SBC_vs_LCDM
     TT     215     2.555553e+20    2.554071e+20         0.001059
     EE     199     2.434993e+18    2.436134e+18         0.001055
     TE     199     6.092748e+18    6.092791e+18         0.001309

## Interpretation

Delta chi2 = -0.000000

This replaces the earlier simplified CMB proxy with a more realistic compressed
TTTEEE covariance-level Planck-lite diagnostic.
