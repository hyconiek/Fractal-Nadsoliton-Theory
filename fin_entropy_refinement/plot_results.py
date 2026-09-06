"""Render the actual ST8547 counterexample and correlation calculations."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import analysis

root = Path(__file__).resolve().parent
data = json.loads((root/"results.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.8), constrained_layout=True)
t = np.linspace(-7, 0, 500)
x = 2.**t
axes[0].plot(t, x*analysis.periodic_ddphi(x), label="dyadic-compatible adversary", color="#b34d29")
axes[0].axhline(1., color="#184e77", linestyle="--", label="Shannon")
axes[0].set(xlabel=r"$\log_2 x$", ylabel=r"$x\,\phi''(x)$",
            title="Binary refinement leaves a periodic freedom")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=.2)
sweep = data["log_mean_refinement"]["sweep"]
eps = [s["epsilon"] for s in sweep]
axes[1].plot(eps, [s["entropy_gap"] for s in sweep], "o-", color="#184e77",
             label="hidden horizontal dissipation")
axes[1].axhline(0, color="#888888", linestyle="--",
                label="change in coarse heat law: zero")
axes[1].set(xlabel="correlation parameter", ylabel="dimensionless entropy dissipation",
            title="Same coarse dynamics, different fine dissipation")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=.2)
fig.savefig(root/"findings.pdf")
fig.savefig(root/"findings.png", dpi=170)
plt.close(fig)
