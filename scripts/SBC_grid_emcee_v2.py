import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
import emcee

# =========================
# Load grid
# =========================

df = pd.read_csv("output/grid_suppressing_full/grid_full_results.csv")

print("Columns:", df.columns.tolist())
print("Grid shape:", df.shape)

points = df[["epsilon", "k0", "n"]].to_numpy()

# =========================
# Interpolators
# =========================

interp_sigma8 = LinearNDInterpolator(points, df["sigma8"].to_numpy())
interp_TT = LinearNDInterpolator(points, df["TT_RMS"].to_numpy())
interp_EE = LinearNDInterpolator(points, df["EE_RMS"].to_numpy())
interp_TE = LinearNDInterpolator(points, df["TE_robust_RMS"].to_numpy())
interp_PP = LinearNDInterpolator(points, df["PP_RMS"].to_numpy())

# Optional columns
interp_Pk = None
interp_S8 = None

if "Pk_RMS" in df.columns:
    interp_Pk = LinearNDInterpolator(points, df["Pk_RMS"].to_numpy())

if "S8" in df.columns:
    interp_S8 = LinearNDInterpolator(points, df["S8"].to_numpy())

# =========================
# Targets / tolerances
# =========================

sigma8_target = 0.8228
sigma8_err = 0.0005

S8_target = 0.776
S8_err = 0.017

# Proxy CMB tolerance
cmb_tol = 1.0e-3

# P(k) tolerance, if available
pk_tol = 1.0e-2

# =========================
# Likelihood
# =========================

def safe_interp(interp, eps, k0, n):
    v = interp(eps, k0, n)
    if v is None:
        return None
    v = float(np.asarray(v))
    if not np.isfinite(v):
        return None
    return v

def log_prior(theta):
    eps, k0, n = theta

    if 0.06 <= eps <= 0.12 and 0.015 <= k0 <= 0.04 and 0.0 <= n <= 1.0:
        return 0.0

    return -np.inf

def log_likelihood(theta):
    eps, k0, n = theta

    TT = safe_interp(interp_TT, eps, k0, n)
    EE = safe_interp(interp_EE, eps, k0, n)
    TE = safe_interp(interp_TE, eps, k0, n)
    PP = safe_interp(interp_PP, eps, k0, n)
    sigma8 = safe_interp(interp_sigma8, eps, k0, n)

    if any(v is None for v in [TT, EE, TE, PP, sigma8]):
        return -np.inf

    chi2 = 0.0

    # CMB residual proxy
    chi2 += (TT / cmb_tol)**2
    chi2 += (EE / cmb_tol)**2
    chi2 += (TE / cmb_tol)**2
    chi2 += (PP / cmb_tol)**2

    # Growth amplitude
    chi2 += ((sigma8 - sigma8_target) / sigma8_err)**2

    # Matter power proxy if available
    if interp_Pk is not None:
        Pk = safe_interp(interp_Pk, eps, k0, n)
        if Pk is None:
            return -np.inf
        chi2 += (Pk / pk_tol)**2

    # S8 if available
    if interp_S8 is not None:
        S8 = safe_interp(interp_S8, eps, k0, n)
        if S8 is None:
            return -np.inf
        chi2 += ((S8 - S8_target) / S8_err)**2

    return -0.5 * chi2

def log_prob(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll

# =========================
# Run MCMC
# =========================

np.random.seed(42)

ndim = 3
nwalkers = 32
nsteps = 300

initial = np.array([0.08, 0.02, 0.5])
scale = np.array([0.01, 0.004, 0.15])

p0 = initial + scale * np.random.randn(nwalkers, ndim)

# force inside prior
p0[:, 0] = np.clip(p0[:, 0], 0.0601, 0.1199)
p0[:, 1] = np.clip(p0[:, 1], 0.0151, 0.0399)
p0[:, 2] = np.clip(p0[:, 2], 0.0001, 0.9999)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
sampler.run_mcmc(p0, nsteps, progress=True)
# MCMC diagnostic: integrated autocorrelation time
try:
    tau = sampler.get_autocorr_time(tol=0)
    n_eff = nwalkers * nsteps / tau

    tau_report = pd.DataFrame({
        "parameter": ["epsilon", "k0", "n"],
        "tau_integrated": tau,
        "effective_samples": n_eff
    })

    tau_report.to_csv(
        "likelihood/reports/SBC_grid_emcee_v2_autocorr.csv",
        index=False
    )

    print("\nAUTOCORRELATION TIME:")
    print(tau_report)

except Exception as e:
    print("\nAutocorrelation diagnostic could not be estimated reliably:")
    print(e)
samples = sampler.get_chain(discard=80, flat=True)
logp = sampler.get_log_prob(discard=80, flat=True)

out = pd.DataFrame(samples, columns=["epsilon", "k0", "n"])
out["log_prob"] = logp


os.makedirs("likelihood/chains", exist_ok=True)
os.makedirs("likelihood/plots", exist_ok=True)

chain_file = "likelihood/chains/SBC_grid_emcee_v2_chain.csv"
out.to_csv(chain_file, index=False)

# =========================
# Plots
# =========================

for col in ["epsilon", "k0", "n"]:
    plt.figure()
    plt.hist(out[col], bins=40)
    plt.xlabel(col)
    plt.ylabel("samples")
    plt.title(f"Posterior proxy v2: {col}")
    plt.grid(True, ls=":")
    plt.tight_layout()
    plt.savefig(f"likelihood/plots/SBC_grid_emcee_v2_{col}.png", dpi=200)
    plt.close()

plt.figure()
plt.scatter(out["k0"], out["epsilon"], c=out["log_prob"], s=8)
plt.xlabel("k0")
plt.ylabel("epsilon")
plt.title("Proxy posterior v2: epsilon vs k0")
plt.colorbar(label="log_prob")
plt.grid(True, ls=":")
plt.tight_layout()
plt.savefig("likelihood/plots/SBC_grid_emcee_v2_eps_k0.png", dpi=200)
plt.close()

plt.figure()
plt.plot(out["log_prob"].to_numpy(), lw=0.7)
plt.xlabel("sample")
plt.ylabel("log probability")
plt.title("IR perturbative sector MCMC log probability trace")
plt.grid(True, ls=":")
plt.tight_layout()
plt.savefig("likelihood/plots/SBC_grid_emcee_v2_logprob_trace.png", dpi=200)
plt.close()

print("Saved:")
print(chain_file)
print("likelihood/plots/SBC_grid_emcee_v2_*.png")

print("\nSUMMARY:")
print(out.describe())

best = out.loc[out["log_prob"].idxmax()]
print("\nBEST POSTERIOR POINT v2:")
print(best)

