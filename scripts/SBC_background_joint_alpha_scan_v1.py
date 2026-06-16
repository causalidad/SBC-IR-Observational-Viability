import numpy as np
import pandas as pd
from pathlib import Path
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
import matplotlib.pyplot as plt

outdir = Path("likelihood/results")
plotdir = Path("likelihood/plots")
repdir = Path("likelihood/reports")
outdir.mkdir(parents=True, exist_ok=True)
plotdir.mkdir(parents=True, exist_ok=True)
repdir.mkdir(parents=True, exist_ok=True)

# Cosmology
c_km_s = 299792.458
h = 0.6736
H0 = 100.0 * h
omega_b = 0.02237
omega_cdm = 0.1200
Omega_m = (omega_b + omega_cdm) / h**2
Omega_L = 1.0 - Omega_m
rd = 147.09

epsilon = 0.060
k0 = 0.015
n = 1.0

def E_lcdm(z):
    return np.sqrt(Omega_m * (1 + z)**3 + Omega_L)

def F_sbc(z, alpha):
    return -alpha * epsilon * np.exp(-k0 * (1 + z)**n)

def E_sbc(z, alpha):
    return E_lcdm(z) * (1 + F_sbc(z, alpha))

def DM(z, Efunc):
    integral, _ = quad(lambda zp: 1.0 / Efunc(zp), 0, z, epsabs=1e-8, epsrel=1e-8)
    return (c_km_s / H0) * integral

def DH(z, Efunc):
    return c_km_s / (H0 * Efunc(z))

def DV(z, Efunc):
    return (z * DM(z, Efunc)**2 * DH(z, Efunc))**(1/3)

def DL(z, Efunc):
    return (1 + z) * DM(z, Efunc)

def mu(z, Efunc):
    return 5*np.log10(DL(z, Efunc)) + 25

# BAO data
bao_mean = "data/bao/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt"
bao_covp = "data/bao/bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt"

bao = pd.read_csv(
    bao_mean,
    sep=r"\s+",
    comment="#",
    header=None,
    names=["z","value","observable"]
)

bao_cov = np.loadtxt(bao_covp)
bao_icov = np.linalg.inv(bao_cov)
bao_data = bao["value"].to_numpy()

def predict_bao(Efunc):
    pred = []
    for _, row in bao.iterrows():
        z = float(row["z"])
        obs = row["observable"]

        if obs == "DM_over_rs":
            pred.append(DM(z, Efunc) / rd)
        elif obs == "DH_over_rs":
            pred.append(DH(z, Efunc) / rd)
        elif obs == "DV_over_rs":
            pred.append(DV(z, Efunc) / rd)
        else:
            raise ValueError(obs)

    return np.array(pred)

def chi2_bao(Efunc):
    delta = predict_bao(Efunc) - bao_data
    return float(delta @ bao_icov @ delta)

# Pantheon data
sn_path = "data/pantheon_plus/DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat"
sn_covp = "data/pantheon_plus/DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov"

df_all = pd.read_csv(sn_path, sep=r"\s+")
mask = df_all["IS_CALIBRATOR"].to_numpy() == 0
sn = df_all[mask].copy()

raw_cov = np.loadtxt(sn_covp)
if raw_cov.ndim == 1:
    n_cov = int(raw_cov[0])
    cov_all = raw_cov[1:].reshape((n_cov, n_cov))
else:
    cov_all = raw_cov

sn_cov = cov_all[np.ix_(mask, mask)]
cho = cho_factor(sn_cov, lower=True, check_finite=False)

z_sn = sn["zHD"].to_numpy()
mu_obs = sn["MU_SH0ES"].to_numpy()
one = np.ones_like(mu_obs)

def chi2_pantheon(Efunc):
    mu_raw = np.array([mu(float(zi), Efunc) for zi in z_sn])
    r0 = mu_raw - mu_obs

    Cinv_r0 = cho_solve(cho, r0, check_finite=False)
    Cinv_one = cho_solve(cho, one, check_finite=False)

    offset = - (one @ Cinv_r0) / (one @ Cinv_one)

    res = r0 + offset
    Cinv_res = cho_solve(cho, res, check_finite=False)

    return float(res @ Cinv_res), offset

alphas = np.linspace(0.0, 0.20, 81)

chi2_bao_lcdm = chi2_bao(E_lcdm)
chi2_sn_lcdm, offset_lcdm = chi2_pantheon(E_lcdm)

rows = []

for alpha in alphas:
    Efunc = lambda zz, a=alpha: E_sbc(zz, a)

    cb = chi2_bao(Efunc)
    cs, off = chi2_pantheon(Efunc)

    rows.append({
        "alpha_bg": alpha,
        "chi2_BAO": cb,
        "Delta_chi2_BAO": cb - chi2_bao_lcdm,
        "chi2_Pantheon": cs,
        "Delta_chi2_Pantheon": cs - chi2_sn_lcdm,
        "offset_Pantheon": off,
        "chi2_joint": cb + cs,
        "Delta_chi2_joint": (cb + cs) - (chi2_bao_lcdm + chi2_sn_lcdm),
    })

out = pd.DataFrame(rows)
best = out.loc[out["chi2_joint"].idxmin()]

csv_out = outdir / "SBC_background_joint_alpha_scan_v1.csv"
out.to_csv(csv_out, index=False)

plt.figure()
plt.plot(out["alpha_bg"], out["Delta_chi2_joint"], label="BAO+Pantheon")
plt.plot(out["alpha_bg"], out["Delta_chi2_BAO"], label="BAO")
plt.plot(out["alpha_bg"], out["Delta_chi2_Pantheon"], label="Pantheon")
plt.axhline(0, ls=":")
plt.axvline(best["alpha_bg"], ls=":")
plt.xlabel("alpha_bg")
plt.ylabel("Delta chi2")
plt.title("BAO + Pantheon+ background scan")

plt.legend()
plt.grid(True, ls=":")
plt.tight_layout()
plt.savefig(plotdir / "SBC_background_joint_alpha_scan_v1.png", dpi=200)
plt.close()

report = f"""# IR Perturbative Sector Background Joint Alpha Scan

## Setup

epsilon = {epsilon}
k0      = {k0}
n       = {n}
alpha range = 0.0 to 0.20

## LCDM reference

chi2_BAO_LCDM = {chi2_bao_lcdm:.6f}
chi2_Pantheon_LCDM = {chi2_sn_lcdm:.6f}
chi2_joint_LCDM = {chi2_bao_lcdm + chi2_sn_lcdm:.6f}

## Best joint point

{best.to_string()}

## Interpretation

This scan determines whether the phenomenological background deformation can improve Pantheon+
without being excluded by BAO.

A viable background branch must keep Delta chi2_BAO small while allowing any
Pantheon+ improvement.
"""

md_out = repdir / "SBC_background_joint_alpha_scan_v1.md"
md_out.write_text(report)

print("Saved:")
print(csv_out)
print(plotdir / "SBC_background_joint_alpha_scan_v1.png")
print(md_out)
print()
print("BEST:")
print(best)
