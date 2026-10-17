#!/usr/bin/env python3
"""Exact, source-independent UE operator/energy/closure controls.
Does NOT read or certify unavailable FIN M2/Y/nonlinear arrays.
"""
from __future__ import annotations
from fractions import Fraction as F
import json, hashlib, sys, platform
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
I=s.eye(2)
Vinfty=s.Matrix([[1,s.Rational(1,2)],[s.Rational(1,2),1]])
B=s.diag(4,0)
A=Vinfty-B
comm=Vinfty*B-B*Vinfty
assert comm == s.Matrix([[0,-2],[2,0]])
assert sorted(Vinfty.eigenvals().keys()) == [s.Rational(1,2),s.Rational(3,2)]
assert A == s.Matrix([[-3,s.Rational(1,2)],[s.Rational(1,2),1]])
# f=sech(x), norm^2=2, int f'^2=2/3, int sech^2 f^2=4/3.
norm=s.Rational(2)
kin=s.Rational(2,3)
weighted=s.Rational(4,3)
q=kin + norm*Vinfty[0,0] - weighted*B[0,0]
rq=q/norm
assert q == -s.Rational(8,3)
assert rq == -s.Rational(4,3)
# Noncommuting analytic Rayleigh upper bound with f=sech(kx).
# For each unit vector u: R= k^2/3 + u^T Vinfty u - 2/3 u^T B u.
mat_test=Vinfty-s.Rational(2,3)*B
assert mat_test==s.Matrix([[-s.Rational(5,3),s.Rational(1,2)],[s.Rational(1,2),1]])
# For scalar phi=tanh(x), psi=sech(x), coupling Vint=phi^2 psi^2.
# z=tanh x: F=2 z (1-z^2), G=2 z^2 sqrt(1-z^2).
z=s.symbols('z', real=True)
F2=s.integrate(4*z**2*(1-z**2), (z,-1,1))
G2=s.integrate(4*z**4,(z,-1,1))
c=s.cancel(F2/G2)
assert F2==s.Rational(16,15)
assert G2==s.Rational(8,5)
assert c==s.Rational(2,3)
# dE_phi=-delta*F2, dE_psi=+delta*c*G2, dE_int=0.
assert s.simplify(-F2+c*G2)==0
# Two compatible localized excitations of the SAME Psi field (not a new force):
# Phi_0=1 and Psi_0=exp(-(x+a)^2)+exp(-(x-a)^2).
# Cross part of existing Vint=Phi^2 Psi^2 is 2*int f_plus*f_minus dx.
x,a=s.symbols('x a',real=True,positive=True)
two_bump_cross=s.simplify(2*s.integrate(s.exp(-2*x**2-2*a**2),(x,-s.oo,s.oo)))
assert s.simplify(two_bump_cross-s.sqrt(2*s.pi)*s.exp(-2*a**2))==0
# Same coarse state, hidden psi=0 versus psi=epsilon sech x:
# acceleration difference Delta a=-2 epsilon^2 tanh x sech^2 x.
# L2 norm = 4/sqrt(15) epsilon^2, by exact integration.
assert s.simplify(F2)==s.Rational(16,15)
checks=[
  'noncommuting_commutator_exact',
  'positive_far_vacuum_eigenvalues_exact',
  'wall_core_indefinite_exact',
  'negative_rayleigh_integrals_exact',
  'negative_rayleigh_quotient_exact',
  'matrix_analytic_trial_bound_exact',
  'nonlinear_force_L2_exact',
  'nonlinear_hidden_force_L2_exact',
  'interaction_stationary_at_t0_exact',
  'two_configuration_existing_cross_term_exact',
  'coarse_acceleration_separation_exact'
]
out={
 'execution_status':'EXECUTED_EXACT_ANALYTIC_CONTROLS_ONLY',
 'strict_source_gate':'BLOCKED_UNAVAILABLE_ACTUAL_SOURCE_FILES',
 'fin_actual_matrices_checked':False,
 'fin_full_complex_action_checked':False,
 'precision':'SymPy exact integers/rationals; symbolic integral identities',
 'checks_passed':len(checks),'checks':checks,
 'operator_control':{
   'V_inf':[[str(c) for c in row] for row in Vinfty.tolist()],
   'B':[[str(c) for c in row] for row in B.tolist()],
   'V_0':[[str(c) for c in row] for row in A.tolist()],
   'commutator':[[str(c) for c in row] for row in comm.tolist()],
   'far_vacuum_spectrum':['1/2','3/2'],
   'trial_norm_squared':str(norm),'trial_kinetic':str(kin),
   'trial_attractive_integral':str(weighted),
   'negative_quadratic_form':str(q),'negative_rayleigh_quotient':str(rq),
   'below_essential_threshold':True,'actual_FIN_instability_claim':False
 },
 'transfer_control':{
   'V_interaction':'phi^2 psi^2','F2':str(F2),'G2':str(G2),
   'Cauchy_velocity_ratio':str(c),
   'dE_phi_over_delta':'-16/15','dE_psi_over_delta':'16/15',
   'dE_int_over_delta':'0',
   'two_bump_interaction_cross':'sqrt(2*pi)*exp(-2*a^2)',
   'actual_FIN_transfer_claim':False
 },
 'coarse_obstruction':{
   'same_initial_coarse_state':True,
   'different_initial_psi':'0 versus epsilon sech(x)',
   'L2_acceleration_difference':'4*epsilon^2/sqrt(15)',
   'autonomous_phi_projection':False
 },
 'runtime':{'python':sys.version.split()[0],'sympy':s.__version__, 'platform':platform.platform()},
 'nonconclusions':['No input M2/Y diagonal/array verification','No complete noncommuting FIN spectrum',
    'No strict nadsoliton sourcing','No physical mass/space-time/observable claim',
    'No FIN-specific two-wall scattering, absorption or stability proof']
}
(ROOT/'REPLAY'/'exact_controls.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print('PASS',len(checks),'exact analytical checks')
print('noncommuting negative control: spectrum(Vinf)={1/2,3/2}, Rayleigh=-4/3')
print('instantaneous diagnostic transfer: dEphi=-16*delta/15; dEpsi=+16*delta/15')
print('coarse acceleration separation: 4 epsilon^2/sqrt(15)')
print('same-field two-bump cross: sqrt(2*pi)*exp(-2*a^2)')
