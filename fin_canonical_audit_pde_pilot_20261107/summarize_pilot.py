"""Summarize reproducible full-PDE runs; retain separate sampled trajectories."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent
raw=json.loads((H/'pde_results.json').read_text());runs=raw['runs']
assert len(runs)==9,'Wait for all runs before summarizing'
rows=[];summaries=[]
for i,v in enumerate(runs):
    hist=v['history'];first,last=hist[0],hist[-1]
    out={k:x for k,x in v.items() if k!='history'}
    out.update({'run_id':i,'final_I':last['I'],'initial_I':first['I'],
        'bound_energy_loss':np.sqrt(v['w'])*(first['I']-last['I'])+np.sqrt(3)*(first['J']-last['J']),
        'continuum_energy_gain':last['continuum_quadratic_energy']-first['continuum_quadratic_energy'],
        'interaction_energy_change':last['interaction_and_chart_energy']-first['interaction_and_chart_energy']})
    summaries.append(out)
    for j,q in enumerate(hist):
        if j%4==0 or j==len(hist)-1:rows.append({'run_id':i,**q})
with (H/'trajectories.csv').open('w') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
base=summaries[5]
metrics={'time_halving_slope_change':abs(summaries[6]['finite_window_inverse_action_slope']-base['finite_window_inverse_action_slope']),
 'domain_enlargement_slope_change':abs(summaries[7]['finite_window_inverse_action_slope']-base['finite_window_inverse_action_slope']),
 'space_halving_slope_change_at_same_dt':abs(summaries[8]['finite_window_inverse_action_slope']-summaries[6]['finite_window_inverse_action_slope'])}
checks={'all_finite':all(np.isfinite(v['final_I']) for v in summaries),
 'energy_drift_below_4e_minus6':max(v['relative_energy_drift'] for v in summaries)<4e-6,
 'time_refinement_slope_below_1e_minus6':metrics['time_halving_slope_change']<1e-6,
 'domain_refinement_slope_below_1e_minus8':metrics['domain_enlargement_slope_change']<1e-8,
 'spatial_refinement_slope_below_1e_minus4':metrics['space_halving_slope_change_at_same_dt']<1e-4,
 'positive_continuum_gain':all(v['continuum_energy_gain']>0 for v in summaries)}
# These are declared numerical quality gates, not theorem acceptance tests.
(H/'pilot_summary.json').write_text(json.dumps({'status':'PILOT_ONLY','runs':summaries,'refinement_metrics':metrics,'quality_checks':checks},indent=2)+'\n')
fig,ax=plt.subplots(1,2,figsize=(10,3.8))
for idx in [1,3,5]:
    v=runs[idx];hist=v['history'];time=np.array([q['t'] for q in hist])*v['epsilon']**2
    ax[0].plot(time,[q['I']/hist[0]['I'] for q in hist],label=f"PDE epsilon={v['epsilon']}")
v=runs[5];hist=v['history'];ax[0].plot(np.array([q['t'] for q in hist])*.01,[q['leading_ode_I']/hist[0]['I'] for q in hist],'k--',label='Leading action ODE (epsilon=0.1)')
ax[0].set(xlabel='Slow time epsilon^2 t',ylabel='Canonical action I / I(0)');ax[0].legend(fontsize=7)
t=np.array([q['t'] for q in hist]);E0=hist[0]['total_energy']
ax[1].plot(t,[q['continuum_quadratic_energy']/E0 for q in hist],label='Continuum quadratic energy')
ax[1].plot(t,[(np.sqrt(1.5)*(hist[0]['I']-q['I'])+np.sqrt(3)*(hist[0]['J']-q['J']))/E0 for q in hist],label='Canonical bound energy loss')
ax[1].plot(t,[(q['interaction_and_chart_energy']-hist[0]['interaction_and_chart_energy'])/E0 for q in hist],label='Interaction / chart change')
ax[1].set(xlabel='Model time (epsilon=0.1)',ylabel='Fraction of initial total energy');ax[1].legend(fontsize=7)
for a in ax:a.grid(alpha=.25)
fig.suptitle('Full nonlinear PDE pilot: fixed w=1.5, g4=g6=1; no damping fit in evolution')
fig.tight_layout();fig.savefig(H/'pde_comparison.svg');plt.close(fig)
print(json.dumps({'refinements':metrics,'checks':checks,'fine_runs':[summaries[i] for i in [1,3,5,8]]},indent=2))
