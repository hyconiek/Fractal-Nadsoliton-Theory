#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any
import sympy as sp
import numpy as np
import scipy.integrate as si

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
OUT=GEN/'p1950_s900_strict_renormalization_exact_integration_probe.json'
MD=GEN/'p1950_s900_strict_renormalization_exact_integration_theorem.md'

def strict_kernel(d,omega,phi,beta,eta):
    return sp.cos(omega*d+phi)/(1+beta*d**eta)



def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required backend artifact: {path}")
    return json.loads(path.read_text(encoding='utf-8'))


def _direct_profile_exprs(
    direct_profiles: dict[str, Any],
    d: sp.Symbol,
    K: sp.Expr,
    Kd: sp.Expr,
    Kdd: sp.Expr,
) -> dict[str, sp.Expr] | None:
    profiles = direct_profiles.get('profiles') or {}
    required = ('R2', 'Ric2', 'Riem2', 'GB')
    if not all(key in profiles for key in required):
        return None

    expr_by_channel: dict[str, sp.Expr] = {}
    base_locals = {
        'd': d,
        'K': K,
        'Kd': Kd,
        'Kdd': Kdd,
        'Abs': sp.Abs,
        'sqrt': sp.sqrt,
        'sin': sp.sin,
        'cos': sp.cos,
        'exp': sp.exp,
        'log': sp.log,
    }
    for channel in required:
        raw = profiles[channel]
        expr_text = raw.get('expression') if isinstance(raw, dict) else raw
        if not isinstance(expr_text, str) or not expr_text.strip():
            return None
        local_exprs = {**base_locals, **expr_by_channel}
        expr_by_channel[channel] = sp.sympify(expr_text, locals=local_exprs)
    return expr_by_channel


def _backend_channel_defs(d: sp.Symbol, K: sp.Expr, Kd: sp.Expr, Kdd: sp.Expr) -> tuple[list[tuple[str, str, sp.Expr]], dict[str, Any]]:
    p1848 = _load_json(GEN / 'p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json')
    projection = ((p1848.get('gravity_counterterm_projection_pack') or {}).get('strict_gravity_projection') or {})
    if not projection:
        raise ValueError('Backend strict_gravity_projection is empty; cannot bind delta_c_gr_i channels.')

    channel_label = {'R^2': 'R2', 'Ricci^2': 'Ric2', 'Riemann^2': 'Riem2', 'GaussBonnet': 'GB'}
    direct_profiles = (p1848.get('gravity_operator_profiles_B1') or {})
    direct_expr_by_channel = _direct_profile_exprs(direct_profiles, d, K, Kd, Kdd)
    if direct_expr_by_channel is None:
        # Fallback is retained only to keep the obstruction witness executable if
        # the backend profile object disappears.  The gate below remains OPEN
        # whenever the direct backend profile export is incomplete.
        direct_profiles_present = False
        profile_source_kind = 'fallback_surrogate_profiles_missing_backend_direct_export'
        profile_by_channel = {
            'R2': sp.simplify(K**2),
            'Ric2': sp.simplify(Kd**2),
            'Riem2': sp.simplify(sp.Abs(K * Kdd)),
            'GB': sp.simplify((K**2) / (1 + d)),
        }
    else:
        direct_profiles_present = True
        profile_source_kind = 'p1848_direct_gravity_operator_profiles_B1'
        profile_by_channel = direct_expr_by_channel

    defs: list[tuple[str, str, sp.Expr]] = []
    for operator_name, delta_label in projection.items():
        if operator_name not in channel_label:
            raise ValueError(f'Unsupported backend operator in strict_gravity_projection: {operator_name}')
        channel = channel_label[operator_name]
        defs.append((delta_label, channel, profile_by_channel[channel]))

    defs.sort(key=lambda row: row[0])
    component_pack = p1848.get('gravity_componentwise_variation_pack') or {}
    backend_meta = {
        'projection_labels_source': 'p1848::gravity_counterterm_projection_pack.strict_gravity_projection',
        'tensor_operator_source': 'p1848::gravity_componentwise_variation_pack',
        'tensor_operator_exports_present': all(k in component_pack for k in ['H_R2_munu','H_Ric2_munu','H_Riem2_munu','H_GB_munu']),
        'tensor_operator_exports': {
            'R2': component_pack.get('H_R2_munu'),
            'Ric2': component_pack.get('H_Ric2_munu'),
            'Riem2': component_pack.get('H_Riem2_munu'),
            'GB': component_pack.get('H_GB_munu'),
        },
        'direct_operator_profiles_present': direct_profiles_present,
        'direct_operator_profiles_source': 'p1848::gravity_operator_profiles_B1',
        'direct_operator_profile_source_kind': profile_source_kind,
        'direct_operator_profile_generation_rule': direct_profiles.get('profile_generation_rule'),
        'direct_operator_profile_expressions': {
            channel: sp.sstr(expr)
            for _, channel, expr in defs
        },
    }
    return defs, backend_meta


def _eval_grid(fn: Any, xs: np.ndarray) -> np.ndarray:
    values = np.asarray(fn(xs), dtype=float)
    if values.shape == ():
        return np.full_like(xs, float(values), dtype=float)
    return values.astype(float, copy=False)
def main():
    omega=sp.Float('0.18575'); phi=sp.Float('0.16250'); beta=sp.Float('1.0'); eta=sp.Rational(9,5)
    d=sp.symbols('d', positive=True)
    eps=sp.Symbol('eps', positive=True)
    K=strict_kernel(d,omega,phi,beta,eta)
    Kd=sp.diff(K,d)
    Kdd=sp.diff(Kd,d)
    integrand=sp.lambdify(d,sp.simplify(K**2),'numpy')
    epsv=1e-6
    m0=si.quad(lambda x: float(integrand(x)),epsv,1.0,epsabs=1e-12,epsrel=1e-12,limit=400)[0]
    m1=si.quad(lambda x: float(x*integrand(x)),epsv,1.0,epsabs=1e-12,epsrel=1e-12,limit=400)[0]
    m2=si.quad(lambda x: float((x**2)*integrand(x)),epsv,1.0,epsabs=1e-12,epsrel=1e-12,limit=400)[0]
    # Channel labels are pulled from backend tensor projection pack (p1848)
    # and then bound to strict analytical profiles for direct B1 replay.
    channel_defs, backend_meta = _backend_channel_defs(d, K, Kd, Kdd)
    chan_funcs=[sp.lambdify(d,ch_expr,'numpy') for _,_,ch_expr in channel_defs]
    # delta_c_gr_i rows are generated from projected channel basis moments on B1.
    xs=np.linspace(epsv,1.0,12001,dtype=float)
    k2=_eval_grid(integrand,xs)
    basis=[_eval_grid(f,xs) for f in chan_funcs]
    # backend-generated operator matrix A_ij = <b_i,b_j>_B1  (Gram matrix)
    A_num=np.zeros((4,4),dtype=float)
    for i in range(4):
        for j in range(4):
            A_num[i,j]=float(si.simpson((basis[i]*basis[j]),x=xs))
    A=sp.Matrix(A_num)
    # divergence target vector from projected strict UV moments on tensor channels
    b_num=np.array([
        float(si.simpson((basis[0]*k2),x=xs)),
        float(si.simpson((basis[1]*k2),x=xs)),
        float(si.simpson((basis[2]*k2),x=xs)),
        float(si.simpson((basis[3]*k2),x=xs)),
    ],dtype=float)
    singular_values=np.linalg.svd(A_num,compute_uv=False)
    max_singular=float(np.max(singular_values)) if singular_values.size else 0.0
    min_singular=float(np.min(singular_values)) if singular_values.size else 0.0
    rank_tolerance=max(float(A_num.shape[0]),float(A_num.shape[1]))*max_singular*np.finfo(float).eps*100.0
    operator_profile_gram_rank=int(np.sum(singular_values>rank_tolerance)) if max_singular > 0.0 else 0
    operator_profile_gram_nullity=int(A_num.shape[1]-operator_profile_gram_rank)
    operator_profile_gram_condition_number=None if min_singular<=rank_tolerance else float(max_singular/min_singular)
    rcond=(rank_tolerance/max_singular) if max_singular > 0.0 else None
    sol_np, *_ = np.linalg.lstsq(A_num,b_num,rcond=rcond)
    names=['a_R2','a_Ric2','a_Riem2','a_GB']
    delta_names=[x[0] for x in channel_defs]
    channel_names=[x[1] for x in channel_defs]
    coeffs={n:float(v) for n,v in zip(names,sol_np)}
    residual_np=A_num@sol_np-b_num
    res_norm=float(np.linalg.norm(residual_np,ord=2))
    term_rows=[]
    for i in range(A.rows):
        lhs=float(np.dot(A_num[i,:],sol_np))
        rhs=float(b_num[i])
        term_rows.append({
            'delta_c_gr_i':delta_names[i],
            'channel':channel_names[i],
            'lhs_counterterm_combo':lhs,
            'rhs_divergence_target':rhs,
            'signed_residual':lhs-rhs,
            'residual_abs':abs(lhs-rhs),
            'row_operator':{
                'a_R2':float(A[i,0]),
                'a_Ric2':float(A[i,1]),
                'a_Riem2':float(A[i,2]),
                'a_GB':float(A[i,3]),
            }
        })
    term_max_abs=max(r['residual_abs'] for r in term_rows)

    A_norm=float(np.linalg.norm(A_num,ord=2))
    b_norm=float(np.linalg.norm(b_num,ord=2))
    x_norm=float(np.linalg.norm(sol_np,ord=2))
    machine_eps=float(np.finfo(float).eps)
    backward_error_bound=max(1.0,A_norm*x_norm+b_norm)*machine_eps*50.0
    simpson_grid_refined=np.linspace(epsv,1.0,24001,dtype=float)
    k2_ref=_eval_grid(integrand,simpson_grid_refined)
    basis_ref=[_eval_grid(f,simpson_grid_refined) for f in chan_funcs]
    A_num_ref=np.zeros((4,4),dtype=float)
    for i in range(4):
        for j in range(4):
            A_num_ref[i,j]=float(si.simpson((basis_ref[i]*basis_ref[j]),x=simpson_grid_refined))
    b_num_ref=np.array([
        float(si.simpson((basis_ref[0]*k2_ref),x=simpson_grid_refined)),
        float(si.simpson((basis_ref[1]*k2_ref),x=simpson_grid_refined)),
        float(si.simpson((basis_ref[2]*k2_ref),x=simpson_grid_refined)),
        float(si.simpson((basis_ref[3]*k2_ref),x=simpson_grid_refined)),
    ],dtype=float)
    quadrature_discretization_bound=max(float(np.max(np.abs(A_num_ref-A_num))), float(np.max(np.abs(b_num_ref-b_num))))
    residual_tolerance=max(backward_error_bound,quadrature_discretization_bound,1e-18)
    quadrature_scale=max(float(np.max(np.abs(A_num))), float(np.max(np.abs(b_num))), 1.0)
    quadrature_relative_discretization_bound=quadrature_discretization_bound/quadrature_scale
    quadrature_relative_stability_tolerance=1e-6
    missing_direct_profiles_count = 0 if backend_meta['direct_operator_profiles_present'] else 4
    rank_gate_pass = operator_profile_gram_rank == len(channel_defs)
    quadrature_gate_pass = quadrature_relative_discretization_bound <= quadrature_relative_stability_tolerance
    residual_gate_pass = res_norm<=residual_tolerance and term_max_abs<=residual_tolerance
    direct_profile_gate_pass = missing_direct_profiles_count==0
    closed_local_witness = residual_gate_pass and direct_profile_gate_pass and rank_gate_pass and quadrature_gate_pass
    if direct_profile_gate_pass is False:
        fail_trace = f'missing_direct_operator_profiles_count={missing_direct_profiles_count} > 0'
    elif rank_gate_pass is False:
        fail_trace = (
            f'operator_profile_gram_rank={operator_profile_gram_rank} < {len(channel_defs)}; '
            f'nullity={operator_profile_gram_nullity}; condition_number=rank_deficient'
        )
    elif quadrature_gate_pass is False:
        fail_trace = (
            f'quadrature_relative_discretization_bound={quadrature_relative_discretization_bound:.6e} '
            f'> {quadrature_relative_stability_tolerance:.6e}'
        )
    elif res_norm>residual_tolerance:
        fail_trace = f'linear_residual_l2={res_norm:.6e} > residual_tolerance={residual_tolerance:.6e}'
    elif term_max_abs>residual_tolerance:
        fail_trace = f'termwise_residual_abs_max={term_max_abs:.6e} > residual_tolerance={residual_tolerance:.6e}'
    else:
        fail_trace = ''
    witness={
      'schema_version':'1.0.0',
      'theorem_target':'STRICT_B1_ONE_LOOP_DIVERGENCE_EXACT_ALGEBRAIC_CANCELLATION',
      'strict_lane_assumptions':['strict_kernel_only','no_legacy_transfer','background_family_B1_backend_projection_channels'],
      'domain':{
          'background_family':'B1',
          'epsilon_uv':epsv,
          'kernel_params':{'omega':float(omega),'phi':float(phi),'beta':float(beta),'eta':float(eta)},
          'backend_counterterm_basis':[{'delta_label':dl,'channel':cn} for dl,cn,_ in channel_defs],
          'backend_matrix_generation':'gram_projection_on_B1_of_tensor_channel_basis',
          'backend_channel_source':'p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json::gravity_counterterm_projection_pack.strict_gravity_projection',
          'backend_tensor_operator_source':'p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json::gravity_componentwise_variation_pack',
          'backend_tensor_operator_exports_present':backend_meta['tensor_operator_exports_present'],
          'backend_direct_operator_profiles_source':backend_meta['direct_operator_profiles_source'],
          'backend_direct_operator_profiles_present':backend_meta['direct_operator_profiles_present'],
          'backend_direct_operator_profile_source_kind':backend_meta['direct_operator_profile_source_kind'],
          'backend_direct_operator_profile_generation_rule':backend_meta['direct_operator_profile_generation_rule'],
      },
      'computed_rows':[
          {'moment':'m0','value':float(m0)},
          {'moment':'m1','value':float(m1)},
          {'moment':'m2','value':float(m2)},
          *term_rows
      ],
      'aggregate_metrics':{
          **coeffs,
          'linear_residual_l2':res_norm,
          'termwise_residual_abs_max':term_max_abs,
          'residual_tolerance':residual_tolerance,
          'backward_error_bound':backward_error_bound,
          'quadrature_discretization_bound':quadrature_discretization_bound,
          'quadrature_scale':quadrature_scale,
          'quadrature_relative_discretization_bound':quadrature_relative_discretization_bound,
          'quadrature_relative_stability_tolerance':quadrature_relative_stability_tolerance,
          'backend_tensor_operator_exports':backend_meta['tensor_operator_exports'],
          'backend_direct_operator_profile_expressions':backend_meta['direct_operator_profile_expressions'],
          'missing_direct_operator_profiles_count':missing_direct_profiles_count,
          'operator_profile_gram_rank':operator_profile_gram_rank,
          'operator_profile_gram_nullity':operator_profile_gram_nullity,
          'operator_profile_gram_singular_values':[float(x) for x in singular_values],
          'operator_profile_gram_rank_tolerance':rank_tolerance,
          'operator_profile_gram_condition_number':operator_profile_gram_condition_number,
      },
      'pass_fail_criteria':{
          'linear_residual_l2_le_residual_tolerance':residual_tolerance,
          'termwise_residual_abs_max_le_residual_tolerance':residual_tolerance,
          'missing_direct_operator_profiles_count_eq_0':0,
          'operator_profile_gram_rank_eq_channel_count':len(channel_defs),
          'quadrature_relative_discretization_bound_le':quadrature_relative_stability_tolerance,
      },
      'verdict':'CLOSED_NUMERICAL_WITNESS_TASK1' if closed_local_witness else 'OPEN_OBSTRUCTION_WITH_TRACE',
      'fail_trace':fail_trace
    }
    GEN.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(witness,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    md=f"""# P1950 S900 Strict Renormalization Exact Integration Theorem

Target: `{witness['theorem_target']}`

We define strict kernel

`K_strict(d)=cos(omega d+phi)/(1+beta d^eta)` with
`omega={float(omega)}, phi={float(phi)}, beta={float(beta)}, eta={float(eta)}`.

B1 proxy UV moments on `[epsilon,1]`, `epsilon=1e-6`:
- m0={float(m0):.16e}
- m1={float(m1):.16e}
- m2={float(m2):.16e}

Solve exact linear cancellation system `A a = b` for counterterm coefficients
`a=(a_R2,a_Ric2,a_Riem2,a_GB)`.

Solution:
- a_R2={coeffs['a_R2']:.16e}
- a_Ric2={coeffs['a_Ric2']:.16e}
- a_Riem2={coeffs['a_Riem2']:.16e}
- a_GB={coeffs['a_GB']:.16e}

Residual norm `||A a - b||_2 = {res_norm:.3e}`.
Termwise max residual `max_i |(A a - b)_i| = {term_max_abs:.3e}`.
Residual tolerance `{residual_tolerance:.3e}` is computed from backward-error and Simpson-refinement bounds.
Quadrature relative discretization bound `{quadrature_relative_discretization_bound:.3e}`
against stability tolerance `{quadrature_relative_stability_tolerance:.3e}`.

Direct backend profile source:
`{backend_meta['direct_operator_profile_source_kind']}`.

Backend profile expressions:
- R2: `{backend_meta['direct_operator_profile_expressions'].get('R2')}`
- Ric2: `{backend_meta['direct_operator_profile_expressions'].get('Ric2')}`
- Riem2: `{backend_meta['direct_operator_profile_expressions'].get('Riem2')}`
- GB: `{backend_meta['direct_operator_profile_expressions'].get('GB')}`

Gram rank: `{operator_profile_gram_rank}` of `{len(channel_defs)}`.
Gram nullity: `{operator_profile_gram_nullity}`.

Verdict: `{witness['verdict']}`.
Fail trace: `{witness['fail_trace']}`.

Interpretation: the direct tensor-profile export is now consumed by the same
witness+gate pipeline.  If the Gauss-Bonnet tensor identity makes the B1 scalar
profile Gram matrix rank-deficient, the local cancellation residual may still
be small, but the four-channel coefficient identification is not theorem-grade
unique on this B1 scalar projection.
"""
    MD.write_text(md,encoding='utf-8')
    print(OUT)
    print(MD)

if __name__=='__main__':
    main()
