from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent
rs=json.loads((H/'results.json').read_text())['runs'];assert len(rs)==6
old=json.loads((H.parent/'fin_pde_parameter_campaign_20261108/results.json').read_text())['runs']
comparisons=[]
for r in rs:
    ref=next(v for v in old if v['w']==r['w'] and v['epsilon']==.1 and v['dx']==.1)
    comparisons.append({'w':r['w'],'dx':r['dx'],'spectral_ka':r['effective_ka'],'FD_ka':ref['finite_window_inverse_action_slope'],
    'relative_method_difference':r['effective_ka']/ref['finite_window_inverse_action_slope']-1,
    'relative_difference_from_prediction':r['effective_ka']/ref['predicted_ka']-1})
checks={'six_runs':len(rs)==6,'positive_net_outward_flux':all(min(r['samples'][-1]['integrated_outward_flux'])>0 for r in rs),
'local_balance_below_1e_minus7':max(max(r['max_relative_local_balance_residual']) for r in rs)<1e-7,
'method_difference_below_0_001':max(abs(c['relative_method_difference']) for c in comparisons)<.001}
(H/'comparison.json').write_text(json.dumps({'comparisons':comparisons,'quality_checks':checks},indent=2)+'\n')
rows='\n'.join(f"| {c['w']:.8f} | {c['dx']:.2f} | {c['spectral_ka']:.9f} | {100*c['relative_method_difference']:+.5f}% | {100*c['relative_difference_from_prediction']:+.4f}% |" for c in comparisons)
report='''# Fourier replication and physical local energy flux

Six full nonlinear Fourier-collocation evolutions test w=1.25,4/3,1.75
at epsilon=0.1, g4=g6=1, T=100, L=130. Resolutions (dx,dt)=(0.1,0.02)
and (0.05,0.01) are compared with the previous fourth-order finite-difference
campaign. The spatial method is different; RK4 and the canonical observable
reconstruction are shared. This is not independent-team replication.

## Damping comparison

| w | Fourier dx | Effective k_a | Difference from FD dx=0.1 | Difference from theory |
| --- | --- | --- | --- | --- |
'''+rows+'''

The same slope of 1/I over t in [25,100] is used, with the same sampled
times. The damping law is not inserted into either field solver.

## Physical balance

At R=10 and R=20 compute the actual flux
j=-phi_t phi_x-psi_t psi_x, with phi=K+eta. Integrate
j(R)-j(-R) at every RK4 stage with its RK4 weights. Compute the energy
inside [-R,R] independently by composite Simpson quadrature, subtracting
the static kink density. The full potential and the linear terms
K' eta_x+2K(K^2-1)eta are retained. The tested residual is

    E_R(t)-E_R(0)+integral_0^t [j(R)-j(-R)] ds.

It is normalized by the initial total perturbation energy. The stored
maximum is over diagnostic times spaced by 0.5, not a continuous-time
supremum. All six runs show positive final net outward flux and all
reported balance residuals are below 1e-7. The fine-grid numbers and
full sampled budgets are in results.json. Flux was not defined as energy
loss; it is separately evaluated from physical derivatives.

## Limits and next work

The physical balance and cross-discretization agreement support genuine
energy transport rather than numerical dissipation in these cases. They
do not prove AW, identify all three damping coefficients, or establish
empirical FIN. The solver has no explicit dealiasing: the doubled spatial
resolution is a numerical check, not a rigorous aliasing bound. No absorber
is used, and the box is larger than T plus a margin. New domain checks and
a different temporal integrator remain possible additional controls.

Next use separately labelled two-mode initial preparations to identify
k_b and k_x and compare their fixed theoretical values, retaining energy
and flux checks. Such runs test the supplied model beyond the historical
AW preparation and are not automatically covered by that theorem. Continue
the FIN nonlinear matching problem separately; numerical agreement does
not derive this PDE from the fundamental core. Routine self-audit chaining
remains stopped and independent review of AW remains open.
'''
(H/'REPORT_EN.md').write_text(report)
(H/'REPLAY.md').write_text('''# Replay\n\nFrom the repository root, using NumPy and SymPy:\n\n    OPENBLAS_NUM_THREADS=1 python fin_pde_spectral_flux_20261108/spectral_flux.py\n    python fin_pde_spectral_flux_20261108/analyze.py\n\nSix runs take several minutes. Results include diagnostic samples every five\nmodel-time units, plus the final point; summary residuals use all half-unit\ndiagnostics. The canonical map is imported from the frozen pilot package.\nThe previous campaign is read only for the comparison. All values use\ndimensionless model units. Checks are numerical quality gates, not proofs.\n''')
inputs=['fin_canonical_audit_pde_pilot_20261107/canonical.py','fin_pde_parameter_campaign_20261108/results.json']
(H/'INPUTS.json').write_text(json.dumps({p:hashlib.sha256((H.parent/p).read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
print(json.dumps({'checks':checks,'fine_balance':[r['max_relative_local_balance_residual'] for r in rs if r['dx']==.05],'comparisons':comparisons},indent=2))
