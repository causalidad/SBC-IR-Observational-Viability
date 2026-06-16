import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

basedir = Path("output/grid_suppressing_full")
outdir = Path("likelihood/plots")
resdir = Path("likelihood/results")

outdir.mkdir(parents=True, exist_ok=True)
resdir.mkdir(parents=True, exist_ok=True)

best_file = basedir / "suppress_eps0p060_k0p015_n1p0_cl_lensed.dat"
lcdm_file = basedir / "lcdm_cl_lensed.dat"

best = np.loadtxt(best_file)
lcdm = np.loadtxt(lcdm_file)

ell = best[:,0]
ell_lcdm = lcdm[:,0]

# Interpolar LCDM a la misma malla ell del archivo SBC
lcdm_interp = np.zeros((len(ell), lcdm.shape[1]))
lcdm_interp[:,0] = ell
for j in range(1, lcdm.shape[1]):
    lcdm_interp[:,j] = np.interp(ell, ell_lcdm, lcdm[:,j])

lcdm = lcdm_interp

names = ["TT","EE","TE"]

best_cols = {
    "TT": 1,
    "EE": 2,
    "TE": 3,
    "PP": 5,
}

lcdm_cols = {
    "TT": 1,
    "EE": 2,
    "TE": 3,
    "PP": 4,
}

rows = []

def rms(x):
    return float(np.sqrt(np.mean(x**2)))

for name in names:

    y_best = best[:, best_cols[name]]
    y_lcdm = lcdm[:, lcdm_cols[name]]

    if name == "TE":
        norm = np.sqrt(np.abs(lcdm[:,1] * lcdm[:,2])) + 1e-30
        res = (y_best - y_lcdm) / norm
    else:
        res = (y_best - y_lcdm) / (np.abs(y_lcdm) + 1e-30)

    for lo, hi in [(2,30),(30,600),(600,2500)]:

        m = (ell >= lo) & (ell <= hi)

        rows.append({
            "spectrum": name,
            "ell_min": lo,
            "ell_max": hi,
            "RMS_best_vs_LCDM": rms(res[m]),
            "max_abs": float(np.max(np.abs(res[m])))
        })

    m = (ell >= 2) & (ell <= 2500)

    plt.figure()
    plt.plot(ell[m], res[m])

    plt.xscale("log")
    plt.xlabel(r"$\ell$")
    plt.ylabel("normalized residual")
    plt.title(f"IR perturbative sector vs LCDM: {name}")
    plt.grid(True, which="both", ls=":")
    plt.tight_layout()


    plt.savefig(
        outdir / f"SBC_best_vs_LCDM_{name}_v1.png",
        dpi=200
    )

    plt.close()

df = pd.DataFrame(rows)

outfile = resdir / "SBC_best_vs_LCDM_real_cls_v1.csv"
df.to_csv(outfile, index=False)

print(df)

print("\nSaved:")
print(outfile)


