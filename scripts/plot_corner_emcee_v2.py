import pandas as pd
import matplotlib.pyplot as plt

try:
    import corner
except ImportError:
    raise SystemExit("Falta instalar corner: pip install corner")

chain_path = "likelihood/chains/SBC_grid_emcee_v2_chain.csv"
out_path = "figures/SBC_corner_posterior_v1.png"

df = pd.read_csv(chain_path)
samples = df[["epsilon", "k0", "n"]].to_numpy()

fig = corner.corner(
    samples,
    labels=[r"$\epsilon$", r"$k_0$", r"$n$"],
    quantiles=[0.16, 0.50, 0.84],
    show_titles=True,
    title_fmt=".3f",
    levels=(0.68, 0.95),
    bins=40,
    smooth=1.0,
    fill_contours=True,
    plot_datapoints=True
)

fig.savefig(out_path, dpi=300, bbox_inches="tight")
print("Saved:", out_path)
