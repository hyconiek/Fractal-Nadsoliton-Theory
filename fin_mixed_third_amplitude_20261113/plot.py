"""Plot the actual coefficient estimates; no fitted limit is presented as data."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent
o=json.loads((H/'comparison.json').read_text());old=json.loads((H.parent/'fin_mixed_scaling_phase_20261113/comparison.json').read_text())['runs']
kx=o['theory'][2]
fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True)
for ax,s in zip(axes,['0.25','0.5']):
    for phase,style in [(0.,'-'),(np.pi/2,'--')]:
        rr=sorted([r for r in old+o['runs'] if r['phase']==phase and r['dx']==.1],key=lambda r:r['epsilon'])
        for i,color in [(0,'#1768ac'),(1,'#c35b10')]:
            ax.plot([r['epsilon']**2 for r in rr],[100*(r['fits'][s]['1'][i]/kx-1) for r in rr],style,marker='o',color=color,label=f"{['I','J'][i]}, phase {'0' if phase==0 else 'pi/2'}")
    fine=next(r for r in o['runs'] if r['dx']==.05)
    for i,color in [(0,'#1768ac'),(1,'#c35b10')]:ax.scatter([.05**2],[100*(fine['fits'][s]['1'][i]/kx-1)],marker='x',s=70,color=color)
    ax.set_title(f'Slow-time window [{s}, 1]');ax.set_xlabel('Squared amplitude');ax.grid(alpha=.2)
axes[0].set_ylabel('Mixed-coefficient deviation from theory (%)');axes[1].legend(fontsize=8)
fig.suptitle('Full nonlinear PDE: three amplitudes; crosses mark refined phase-zero run')
fig.tight_layout();fig.savefig(H/'amplitude_comparison.svg');plt.close(fig)
