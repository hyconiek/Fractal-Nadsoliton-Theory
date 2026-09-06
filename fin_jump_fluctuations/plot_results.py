"""Publication figures from the actual ST8548 action and count formulas."""
from pathlib import Path
import json
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import analysis as a

root=Path(__file__).resolve().parent
d=json.loads((root/"results.json").read_text())["strict_edge"]
forward,reverse=d["a"],d["b"]
traffic=forward+reverse
currents=np.linspace(-2.5,2.5,450)*traffic
exact=np.array([a.current_cost(j,forward,reverse) for j in currents])
quadratic=np.array([a.quadratic_cost(j,forward,reverse) for j in currents])
fig,ax=plt.subplots(1,2,figsize=(10.4,3.9),constrained_layout=True)
ax[0].plot(currents/traffic,exact/traffic,label="Markov jump action",color="#184e77")
ax[0].plot(currents/traffic,quadratic/traffic,"--",label="quadratic log-mean action",color="#b34d29")
ax[0].axvline((forward-reverse)/traffic,color=".65",lw=.8)
ax[0].set(xlabel="current / nominal traffic",ylabel="action / nominal traffic",
          title="Same minimum, different fluctuation cost")
ax[0].legend(fontsize=8)
F=np.linspace(1,4,350)
ax[1].plot(F,1/np.tanh(F/2),label="Markov jumps",color="#184e77")
ax[1].plot(F,2/F,"--",label="quadratic prescription",color="#b34d29")
ax[1].axhline(1,color=".4",ls=":",label="integer-jump lower bound")
ax[1].axvline(math.log(12),color=".65",lw=.8)
ax[1].scatter([math.log(12),math.log(12)],
              [d["shot_variance_to_absolute_drift"],d["quadratic_variance_to_absolute_drift"]],
              c=["#184e77","#b34d29"],zorder=4)
ax[1].set(xlabel="absolute log population ratio",ylabel="local variance / absolute mean",
          title="A microscopic counting consistency test")
ax[1].legend(fontsize=8)
for x in ax:
    x.grid(alpha=.2)
fig.savefig(root/"findings.pdf")
fig.savefig(root/"findings.png",dpi=180)
plt.close(fig)
