#!/usr/bin/env python3
"""FIN ST106--ST117: projection selection, recovery, and adversarial closure.

All physical interpretations remain conditional.  Exact claims concern the
declared finite-dimensional operators, algebras, intervals, and synthetic
probability models only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import minimize_scalar, root

from fin_programs_497_506_next_research import iv_bounds
from fin_st01_st15_research import N, random_orthogonal_fixing_uniform, strict_operator
from fin_st28_st45_research import saturation_energy_gradient_hessians
from fin_st46_st57_research import carrier_probability_table
from fin_st58_st69_research import interval_left_product, interval_matvec
from fin_st70_st81_research import reflection_expansion


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST106_ST117_Results.json"
SUMMARY = ROOT / "FIN_ST106_ST117_Summary.csv"
CERT108 = ROOT / "FIN_ST108_Transcendental_Fold_Certificate.json"
CERT115 = ROOT / "FIN_ST115_Independent_Nuisance_Global_Certificate.json"
PACK117 = ROOT / "FIN_ST117_Adversarial_Projection_Alternatives.json"
FIG_DIR = ROOT / "FIN_ST106_ST117_Figures"
SEED = 20260820


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plus_minus(n: int = N) -> tuple[np.ndarray, np.ndarray]:
    p = np.zeros((2 * n, n)); f = np.zeros((2 * n, n))
    for i in range(n):
        p[i, i] = p[i + n, i] = 1 / math.sqrt(2)
        f[i, i] = 1 / math.sqrt(2); f[i + n, i] = -1 / math.sqrt(2)
    return p, f


def st106_naturality_classification() -> dict:
    return {
        "program": "ST106",
        "object": "Naturality and Symmetry Classification for the Observable Algebra",
        "complex_dimensions": {
            "full_M24": 24**2,
            "layer_swap_commutant_M12_plus_M12": 2 * 12**2,
            "declared_coarse_visible_fine_blind_M12_plus_C": 12**2 + 1,
            "fully_translation_covariant_circulant_commutant": 24,
        },
        "theorem": (
            "The commutant of the layer swap is M12 direct-sum M12 in the plus/minus basis. "
            "It does not select M12 direct-sum C.  Requiring invariance of the fine block under every U(12) "
            "conjugation collapses that block to a scalar by Schur's lemma and yields M12 direct-sum C, but "
            "this requirement is exactly a supplied fine-blindness premise, not a consequence of strict A."
        ),
        "status": "proven_symmetry_classification_and_fine_blindness_premise_identified",
        "boundary": "Naturality does not derive an apparatus, a preferred sector, or a physical projection.",
    }


def st107_lift_facial_geometry(a: np.ndarray) -> dict:
    # The feasible B-set is PSD plus one-sided diagonal and two-sided off-diagonal bounds.
    lower_active = N + N * (N - 1) // 2
    recession_rays = []
    for i in range(N):
        direction = np.zeros((N, N)); direction[i, i] = 1.0
        recession_rays.append({"site": i, "rank": 1, "trace": 1.0})
    return {
        "program": "ST107",
        "object": "Faces and Recession Geometry of the Complete Lift Spectrahedron",
        "ambient_real_symmetric_dimension": N * (N + 1) // 2,
        "active_independent_constraints_at_B_equals_A": lower_active,
        "B_equals_A_is_a_vertex": True,
        "recession_cone": "nonnegative diagonal orthant",
        "recession_extreme_rays": recession_rays,
        "theorem": (
            "Every recession direction D has D_ij=0 for i!=j because all off-diagonal coordinates have finite "
            "two-sided bounds.  Positive-semidefinite recession then gives D_ii>=0.  Hence the recession cone is "
            "the nonnegative diagonal orthant, with exactly twelve rank-one extreme rays E_ii.  At B=A all 78 "
            "independent lower coordinate bounds are active, so B=A is a vertex of the complete feasible set."
        ),
        "unresolved": "A complete enumeration of all bounded faces and vertices was not obtained.",
        "status": "proven_recession_cone_and_one_exact_vertex_partial_facial_classification",
        "boundary": "This convex geometry supplies adversarial lifts, not a law selecting one lift.",
    }


def iv(value: float | str | tuple[float, float]):
    if isinstance(value, tuple):
        return mp.iv.mpf([str(value[0]), str(value[1])])
    if isinstance(value, (float, np.floating)):
        return mp.iv.mpf([str(float(np.nextafter(value, -np.inf))), str(float(np.nextafter(value, np.inf)))])
    return mp.iv.mpf(str(value))


def transcendental_strict_intervals():
    """Exact-decimal parameter intervals for the declared strict profile."""
    mp.iv.dps = 70
    omega, phi, eta = iv("0.18575"), iv("0.16250"), iv(9) / iv(5)
    wd = {d: mp.iv.cos(omega * d + phi) / (1 + iv(d) ** eta) for d in range(1, 7)}
    s = 2 * sum((wd[d] for d in range(1, 6)), iv(0)) + wd[6]
    aa = [[iv(0) for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            aa[i][j] = s if i == j else -wd[min((i-j) % N, (j-i) % N)]
    return aa, wd, s


def st108_transcendental_fold(a: np.ndarray) -> dict:
    previous = json.loads((ROOT / "FIN_ST70_ST81_Results.json").read_text(encoding="utf-8"))["ST77"]
    expansion = reflection_expansion(); selected = np.arange(7)
    x0 = np.r_[previous["fold_state_reduced"], previous["fold_kappa"], previous["fold_null_vector"]].astype(float)

    def system(candidate: np.ndarray) -> np.ndarray:
        q, kappa, vector = candidate[:7], candidate[7], candidate[8:]
        state = expansion @ q
        _, gradient, hessian, _ = saturation_energy_gradient_hessians(kappa * a, state)
        matrix = hessian[np.ix_(selected, np.arange(N))] @ expansion
        return np.r_[gradient[selected], matrix @ vector, 0.5 * (vector @ vector - 1.0)]

    x0 = root(system, x0, method="lm", options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14, "maxiter": 20000}).x
    aa, wd, s_iv = transcendental_strict_intervals()

    def interval_system_jacobian(radius: float, point_a: bool = False):
        q = [iv((x0[i]-radius, x0[i]+radius)) for i in range(7)]
        kappa = iv((x0[7]-radius, x0[7]+radius))
        v = [iv((x0[8+i]-radius, x0[8+i]+radius)) for i in range(7)]
        u, w = [], []
        for site in range(N):
            index = site if site <= 6 else N-site
            u.append(q[index]); w.append(v[index])
        amat = [[iv(a[i, j]) for j in range(N)] for i in range(N)] if point_a else aa
        au = [sum((amat[i][j]*u[j] for j in range(N)), iv(0)) for i in range(N)]
        aw = [sum((amat[i][j]*w[j] for j in range(N)), iv(0)) for i in range(N)]
        h, rdiag, drdu = [], [], []
        for item in u:
            rho = item**2; den = 1 + rho/2
            qfun = rho/den; qp = den**-2; qpp = -(den**-3); qppp = iv("1.5")*den**-4
            hh = -qfun*qp + iv("0.075"); hp = -(qp**2 + qfun*qpp); hpp = -(3*qp*qpp + qfun*qppp)
            h.append(hh); rdiag.append(2*hh + 4*rho*hp); drdu.append(2*item*(6*hp + 4*rho*hpp))
        g = [kappa*au[i] + 2*u[i]*h[i] for i in selected]
        g += [kappa*aw[i] + rdiag[i]*w[i] for i in selected]
        g += [iv("0.5")*(sum((item**2 for item in v), iv(0))-1)]
        jac = [[iv(0) for _ in range(15)] for _ in range(15)]
        for i in range(7):
            for col in range(7):
                total = iv(0)
                for site in range(N):
                    if (site if site <= 6 else N-site) == col:
                        total += kappa*amat[i][site]
                if i == col: total += rdiag[i]
                jac[i][col] = total; jac[7+i][8+col] = total
            jac[i][7] = au[i]; jac[7+i][i] = drdu[i]*w[i]; jac[7+i][7] = aw[i]; jac[14][8+i] = v[i]
        glo = np.array([iv_bounds(x)[0] for x in g]); ghi = np.array([iv_bounds(x)[1] for x in g])
        jlo = np.array([[iv_bounds(x)[0] for x in row] for row in jac]); jhi = np.array([[iv_bounds(x)[1] for x in row] for row in jac])
        return glo, ghi, jlo, jhi

    # The inverse is built from the binary64 midpoint Jacobian; all proof boxes use transcendental A intervals.
    _, _, pjl, pju = interval_system_jacobian(0.0, point_a=True)
    c = np.linalg.inv(0.5*(pjl+pju))
    g0lo, g0hi, _, _ = interval_system_jacobian(0.0)
    attempts, accepted = [], None
    for radius in [1e-8, 3e-9, 1e-9, 3e-10]:
        _, _, jlo, jhi = interval_system_jacobian(radius)
        cglo, cghi = interval_matvec(c, c, g0lo, g0hi)
        ylo = np.nextafter(x0-cghi, -np.inf); yhi = np.nextafter(x0-cglo, np.inf)
        cjlo, cjhi = interval_left_product(c, jlo, jhi)
        mlo, mhi = -cjhi, -cjlo
        for i in range(15):
            mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf); mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
        dlo = np.full(15, -radius); dhi = np.full(15, radius)
        mdlo, mdhi = interval_matvec(mlo, mhi, dlo, dhi)
        klo = np.nextafter(ylo+mdlo, -np.inf); khi = np.nextafter(yhi+mdhi, np.inf)
        margin = float(min(np.min(klo-(x0-radius)), np.min((x0+radius)-khi)))
        row = {"radius": radius, "minimum_strict_inclusion_margin": margin, "included": margin > 0,
               "maximum_Krawczyk_half_width": float(np.max((khi-klo)/2))}
        attempts.append(row)
        if margin > 0 and accepted is None: accepted = row
    kernel_widths = {str(d): iv_bounds(wd[d])[1]-iv_bounds(wd[d])[0] for d in wd}
    certificate = {
        "declared_exact_parameters": {"omega": "0.18575", "phi": "0.16250", "eta": "9/5"},
        "row_sum_interval": list(iv_bounds(s_iv)), "kernel_interval_widths": kernel_widths,
        "center": x0.tolist(), "attempts": attempts, "accepted": accepted,
        "arithmetic": "mpmath.iv 70-digit transcendental intervals and nextafter-enclosed binary64 linear algebra",
        "theorem_scope": (
            "Uniform parametric Krawczyk inclusion for every strict A enclosed by evaluation of the declared exact-decimal "
            "formula.  For each A in that interval family there is one root in the common state box; this does not assert "
            "that all A share one identical root."
        ),
        "proof_boundary": "Source-code interval theorem, not a proof-assistant replay and not a physical bifurcation measurement.",
    }
    CERT108.write_text(json.dumps(native(certificate), indent=2, sort_keys=True), encoding="utf-8")
    return {"program":"ST108", "object":"Transcendental-Kernel Fold Certificate", "certificate_file":CERT108.name,
            "certificate_sha256":sha(CERT108), **certificate,
            "status":"proven_uniform_interval_fold_for_declared_transcendental_kernel" if accepted else "blocked_no_interval_inclusion",
            "boundary":certificate["proof_boundary"]}


def st109_maximal_recoverable_code() -> dict:
    alpha = 0.37
    return {
        "program":"ST109", "object":"Maximal Quantum Code Recoverable after Coarse/Fine Pinching",
        "maximum_code_dimension":N, "achieved_dimension":N, "example_alpha":alpha,
        "encoding":"V psi = sqrt(alpha) P psi + sqrt(1-alpha) F U psi, for arbitrary unitary U on C^12",
        "knill_laflamme":"Q Pc Q = alpha Q and Q Pf Q = (1-alpha) Q",
        "theorem":(
            "The pinching channel has errors Pc and Pf.  The displayed matched-copy code obeys the Knill--Laflamme "
            "conditions and is recovered by measuring the sector and applying the inverse branch isometry.  Its dimension "
            "is twelve.  Any code with both nonzero branches has injective restrictions into each twelve-dimensional "
            "sector, hence dimension at most twelve.  The bound is attained and noncommuting reference families are recoverable."
        ),
        "status":"proven_maximal_twelve_dimensional_recoverable_code",
        "boundary":"The code, alpha, U, sector instrument, and recovery are supplied operational resources.",
    }


def st110_equivariant_measure_nonselection() -> dict:
    probabilities = np.ones(12)/12
    return {
        "program":"ST110", "object":"Equivariant Branch Measure and Nonselection Theorem",
        "unique_invariant_probability":probabilities.tolist(), "entropy_nats":float(-np.sum(probabilities*np.log(probabilities))),
        "theorem":(
            "C12 acts transitively on the twelve branches, so invariance forces equal mass and normalization forces 1/12. "
            "The unique invariant probability law is uniform.  It produces no canonical branch: every deterministic "
            "equivariant selector from a singleton invariant input would have to be fixed by the transitive action, and none is."
        ),
        "status":"proven_unique_uniform_measure_and_deterministic_nonselection",
        "boundary":"A random realization additionally needs a stochastic event and record; QW-2191 remains open.",
    }


def st111_antisymmetric_memory() -> dict:
    gamma, nu = 0.35, 1.7
    frequencies = np.linspace(0,5,41)
    # S=gamma I, Omega=[[0,-nu],[nu,0]], collocated input/output b=e1.
    c = np.array([[gamma,-nu],[nu,gamma]]); b = np.array([1.0,0.0])
    rows=[]
    for omega in frequencies:
        response = b @ np.linalg.inv(1j*omega*np.eye(2)+c) @ b
        rows.append({"omega":float(omega),"real_response":float(np.real(response)),"absolute_response":float(abs(response))})
    return {
        "program":"ST111", "object":"One-Antisymmetric-Block Memory Realization",
        "poles":[{"real":-gamma,"imag":nu},{"real":-gamma,"imag":-nu}],
        "minimum_sampled_real_response":min(r["real_response"] for r in rows), "rows":rows,
        "theorem":(
            "For x_dot=-(S+Omega)x+Bu, y=B^T x with S positive definite and Omega^T=-Omega, the storage "
            "V=||x||^2/2 satisfies V_dot=-x^T S x+u^T y because x^T Omega x=0.  The antisymmetric block creates "
            "oscillatory phase and nonreciprocal internal motion but no active gain in the collocated passive realization."
        ),
        "status":"proven_antisymmetric_memory_remains_passive",
        "boundary":"Active gain requires an indefinite dissipative block, pump, noncollocated port, or nonequilibrium source not derived here.",
    }


def binary_gibbs_probability(a: np.ndarray, beta: float, q: float) -> float:
    # For the dyadic family B(0)=A+2D with D=diag(row sum weights); use the exact construction from the prior result.
    from fin_st28_st45_research import dyadic_lift
    p, f = plus_minus(); b0 = f.T @ dyadic_lift(a,0.0) @ f
    zc = float(np.sum(np.exp(-beta*np.linalg.eigvalsh(a))))
    zf = float(np.sum(np.exp(-beta*np.linalg.eigvalsh(b0))))
    return zc/(zc+math.exp(-2*beta*q)*zf)


def bernoulli_chernoff(p0: float,p1: float,delta: float=0.0):
    a=delta+(1-2*delta)*p0; b=delta+(1-2*delta)*p1
    opt=minimize_scalar(lambda s:a**s*b**(1-s)+(1-a)**s*(1-b)**(1-s),bounds=(0,1),method="bounded")
    return -math.log(opt.fun),float(opt.x)


def st112_optimal_design(a: np.ndarray) -> dict:
    q0,q1,delta_max=.2,.7,.05
    betas=np.linspace(.05,5.0,500)
    rows=[]
    for beta in betas:
        p0=binary_gibbs_probability(a,float(beta),q0); p1=binary_gibbs_probability(a,float(beta),q1)
        info,s=bernoulli_chernoff(p0,p1,delta_max)
        rows.append({"beta":float(beta),"p_q0":p0,"p_q1":p1,"worst_detector_delta":delta_max,"Chernoff_information":info,"s":s})
    best=max(rows,key=lambda x:x["Chernoff_information"])
    p_at_q=best["p_q0"]
    fisher=4*best["beta"]**2*(1-2*delta_max)**2*p_at_q**2*(1-p_at_q)**2/((delta_max+(1-2*delta_max)*p_at_q)*(1-delta_max-(1-2*delta_max)*p_at_q))
    return {
        "program":"ST112","object":"Synthetic Fisher/Chernoff Design under Temperature and Detector Nuisance",
        "hypotheses":{"q0":q0,"q1":q1},"beta_grid":[float(betas[0]),float(betas[-1]),len(betas)],
        "detector_bit_flip_upper_bound":delta_max,"best_grid_design":best,"Fisher_information_at_q0_best_beta":float(fisher),
        "analytic_Fisher_formula":"4 beta^2 (1-2 delta)^2 p^2(1-p)^2/[p_delta(1-p_delta)]",
        "status":"strong_numerical_design_with_exact_fisher_formula",
        "boundary":"The optimum is grid-conditional and synthetic; beta, detector model, preparation, and physical scale are supplied.",
    }


def st113_spin_refinement_category() -> dict:
    depth=5
    return {
        "program":"ST113","object":"Z2 Spin Twists on the Directed Refinement Category",
        "depth":depth,"untwisted_compatible_sequences":2,"independently_twisted_sequences":2**(depth+1),
        "recursion":"h_(n+1)=t_n h_n^2=t_n for h_n,t_n in Z2",
        "theorem":(
            "Refinement maps are noninvertible and therefore form a directed category, not a groupoid.  Untwisted degree-two "
            "pullback leaves only two holonomy histories, both periodic after the first step.  With independent Z2 twists, "
            "every subsequent holonomy equals the newly supplied twist, so all 2^(L+1) histories occur through depth L."
        ),
        "status":"proven_finite_refinement_category_twist_classification",
        "boundary":"Each nontrivial fine holonomy is new boundary data; no strict spin source is derived.",
    }


def st114_locality_algebra() -> dict:
    return {
        "program":"ST114","object":"Observable Algebra Generated by Vertex Locality",
        "complex_dimensions":{"full_vertex_locality_plus_connected_generator":24**2,"layer_swap_fixed":2*12**2,"target_fine_blind":12**2+1},
        "theorem":(
            "All vertex projectors together with any connected weighted adjacency generate every matrix unit along graph "
            "paths and hence the full algebra M24.  Restricting only by layer swap gives M12 direct-sum M12.  Neither "
            "locality nor swap symmetry yields the desired M12 direct-sum C; fine blindness must be imposed separately."
        ),
        "status":"proven_locality_selects_too_large_an_algebra",
        "boundary":"Locality cannot presently derive the projection algebra used by the operational model.",
    }


def st115_independent_nuisance_certificate(a: np.ndarray) -> dict:
    rng=np.random.default_rng(20260818+87); qmat=random_orthogonal_fixing_uniform(rng)
    p=carrier_probability_table(a,np.eye(N),transported=False)/N
    q=carrier_probability_table(a,qmat,transported=False)/N
    u=np.full_like(p,1/p.size); emax=.12
    pp=(1-emax)*p+emax*u; qq=(1-emax)*q+emax*u
    opt=minimize_scalar(lambda s:float(np.sum(pp**s*qq**(1-s))),bounds=(0,1),method="bounded",options={"xatol":1e-15})
    s0=float(opt.x); rad=2e-7
    mp.iv.dps=60; si=iv((s0-rad,s0+rad))
    def arrays(e,p0):
        return [[iv((float(np.nextafter(x,-np.inf)),float(np.nextafter(x,np.inf)))) for x in row] for row in ((1-e)*p0+e*u)]
    ppi,qqi=arrays(emax,p),arrays(emax,q)
    derivative_s=iv(0); second_s=iv(0); dep=iv(0); deq=iv(0)
    for i in range(N):
        for j in range(N):
            term=ppi[i][j]**si*qqi[i][j]**(1-si); lr=mp.iv.log(ppi[i][j]/qqi[i][j])
            derivative_s += term*lr; second_s += term*lr**2
            dep += si*iv(float(u[i,j]-p[i,j]))*ppi[i][j]**(si-1)*qqi[i][j]**(1-si)
            deq += (1-si)*iv(float(u[i,j]-q[i,j]))*ppi[i][j]**si*qqi[i][j]**(-si)
    # Sign change is evaluated at the two scalar endpoints.
    def ds_at(s):
        total=iv(0)
        for i in range(N):
            for j in range(N): total += ppi[i][j]**iv(s)*qqi[i][j]**(1-iv(s))*mp.iv.log(ppi[i][j]/qqi[i][j])
        return list(iv_bounds(total))
    certificate={
        "nuisance_box":[0,emax,0,emax],"corner":[emax,emax],"s_center":s0,"s_interval":[s0-rad,s0+rad],
        "dB_ds_at_left":ds_at(s0-rad),"dB_ds_at_right":ds_at(s0+rad),
        "B_second_s_interval":list(iv_bounds(second_s)),
        "corner_partial_epsilon_P_interval":list(iv_bounds(dep)),"corner_partial_epsilon_Q_interval":list(iv_bounds(deq)),
        "corner_Bhattacharyya_minimum":float(opt.fun),"worst_case_Chernoff_information":float(-math.log(opt.fun)),
        "proof":(
            "For fixed s, the Chernoff coefficient is jointly concave in the two distributions, hence in independent "
            "contaminations.  Its pointwise infimum F over s is concave.  Strict positivity of B_ss and the certified "
            "sign change isolate a unique minimizing s at the upper corner.  The envelope derivatives in both nuisance "
            "coordinates are strictly positive there.  The supporting-hyperplane inequality for concave F then proves "
            "F(eP,eQ)<=F(0.12,0.12) throughout the rectangle; therefore Chernoff information is globally minimized there."
        ),
        "scope":"Frozen synthetic ST87 probability tables with independently bounded uniform contamination.",
    }
    signs=(certificate["dB_ds_at_left"][1]<0<certificate["dB_ds_at_right"][0]
           and certificate["B_second_s_interval"][0]>0 and certificate["corner_partial_epsilon_P_interval"][0]>0
           and certificate["corner_partial_epsilon_Q_interval"][0]>0)
    certificate["accepted"]=bool(signs)
    CERT115.write_text(json.dumps(native(certificate),indent=2,sort_keys=True),encoding="utf-8")
    return {"program":"ST115","object":"Global Independent-Nuisance Chernoff Certificate","certificate_file":CERT115.name,
            "certificate_sha256":sha(CERT115),**certificate,
            "status":"proven_global_independent_box_worst_corner" if signs else "blocked_interval_sign_certificate_failed",
            "boundary":"This closes the declared synthetic rectangle, not unknown laboratory nuisance or model misspecification."}


def st116_cgp_vs_thermal_operations(a: np.ndarray) -> dict:
    vals=np.linalg.eigvalsh(a)
    return {
        "program":"ST116","object":"Typing Obstruction for CGP versus Thermal Operations",
        "system_spectrum":vals.tolist(),
        "attempted_target":"construct a covariant Gibbs-preserving channel outside thermal operations",
        "result":(
            "No valid strict-separation witness was promoted.  CGP is a predicate on the system channel, its Hamiltonian, "
            "and beta.  Membership in thermal operations additionally depends on a declared bath-Hamiltonian class, exact "
            "energy-conserving dilation, ancillary dimensions, and whether catalysts or limits are admitted.  These are absent, "
            "so 'outside TO' is not a well-typed proposition on the current FIN artifacts."
        ),
        "status":"proven_current_artifact_typing_obstruction_strict_separation_remains_open",
        "boundary":"This is not a proof that CGP equals TO and not a physical thermal theorem.",
    }


def st117_adversarial_lifts(a: np.ndarray) -> dict:
    p,f=plus_minus(); rows=[]
    base=p@a@p.T+f@a@f.T
    rng=np.random.default_rng(SEED+117); x=rng.normal(size=(N,N)); x=(x+x.T)/2
    for site,t in [(0,.2),(3,.5),(7,1.0),(11,2.0)]:
        b=a.copy(); b[site,site]+=t; lift=p@a@p.T+f@b@f.T
        rows.append({"site":site,"t":t,"coarse_generator_residual":float(np.linalg.norm(p.T@lift@p-a)),
                     "intertwining_residual":float(np.linalg.norm(lift@p-p@a)),
                     "fine_site_expectation_shift":float(f[:,site]@lift@f[:,site]-f[:,site]@base@f[:,site]),
                     "minimum_eigenvalue":float(np.min(np.linalg.eigvalsh(lift))),
                     "translation_commutator_norm":float(np.linalg.norm(np.roll(np.roll(lift,1,0),1,1)-lift)),
                     "coarse_test_residual":float(abs(np.trace(x@(p.T@lift@p))-np.trace(x@a)))})
    packet={"family":"L_(i,t)=PAP*+F(A+t E_ii)F*, t>=0","rows":rows,
            "theorem":(
                "Each ray starts at the ST107 vertex B=A, stays in the complete graph-lift fiber, and preserves every "
                "conditioned coarse functional-calculus experiment exactly.  For t>0 it breaks cyclic translation and is "
                "outside the one-parameter circulant dyadic family.  The supplied fine-site observable detects exactly t."
            ),
            "scope":"Synthetic adversarial alternatives; sites on one cyclic orbit are operationally identical under coarse access."}
    PACK117.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8")
    return {"program":"ST117","object":"Adversarial Projection Alternatives from Extreme Rays","packet_file":PACK117.name,
            "packet_sha256":sha(PACK117),**packet,"status":"proven_adversarial_noncirculant_lift_family",
            "boundary":"The alternatives falsify coarse uniqueness, not FIN physics; fine apparatus remains supplied."}


def make_figures(results:dict)->None:
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4)); dims=results["ST106"]["complex_dimensions"]; ax.bar(["M24","swap","fine-blind"],[dims["full_M24"],dims["layer_swap_commutant_M12_plus_M12"],dims["declared_coarse_visible_fine_blind_M12_plus_C"]]); ax.set(ylabel="complex algebra dimension",title="ST106 symmetry does not select fine blindness"); fig.tight_layout(); fig.savefig(FIG_DIR/"st106_algebra_dimensions.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); ax.bar(range(12),[1]*12); ax.set(xlabel="site",ylabel="extreme-ray trace",title="ST107 twelve recession extreme rays"); fig.tight_layout(); fig.savefig(FIG_DIR/"st107_recession_rays.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=results["ST108"]["attempts"]; ax.semilogx([r["radius"] for r in rows],[r["minimum_strict_inclusion_margin"] for r in rows],"o-"); ax.axhline(0,color="black",lw=.8); ax.set(xlabel="root-box radius",ylabel="inclusion margin",title="ST108 transcendental-kernel certificate"); fig.tight_layout(); fig.savefig(FIG_DIR/"st108_krawczyk.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); ax.bar(range(12),results["ST110"]["unique_invariant_probability"]); ax.set(xlabel="branch",ylabel="invariant probability",title="ST110 invariant measure is uniform, not selective"); fig.tight_layout(); fig.savefig(FIG_DIR/"st110_uniform_measure.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=results["ST111"]["rows"]; ax.plot([r["omega"] for r in rows],[r["real_response"] for r in rows]); ax.axhline(0,color="black",lw=.8); ax.set(xlabel="frequency",ylabel="Re transfer",title="ST111 antisymmetric memory remains passive"); fig.tight_layout(); fig.savefig(FIG_DIR/"st111_passive_antisymmetric.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); best=results["ST112"]["best_grid_design"]; ax.bar(["q0","q1"],[best["p_q0"],best["p_q1"]]); ax.set(ylim=(0,1),ylabel="detected coarse probability",title=f"ST112 best synthetic beta={best['beta']:.3f}"); fig.tight_layout(); fig.savefig(FIG_DIR/"st112_design.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); cert=results["ST115"]; ax.bar(["dF/deP","dF/deQ"],[cert["corner_partial_epsilon_P_interval"][0],cert["corner_partial_epsilon_Q_interval"][0]]); ax.axhline(0,color="black",lw=.8); ax.set(ylabel="certified lower bound",title="ST115 upper corner supports global worst case"); fig.tight_layout(); fig.savefig(FIG_DIR/"st115_nuisance_gradient.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=results["ST117"]["rows"]; ax.plot([r["t"] for r in rows],[r["fine_site_expectation_shift"] for r in rows],"o-",label="fine response"); ax.plot([r["t"] for r in rows],[r["coarse_generator_residual"] for r in rows],"s--",label="coarse residual"); ax.set(xlabel="adversarial ray t",ylabel="response",title="ST117 coarse-indistinguishable alternatives"); ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR/"st117_adversarial.png",dpi=190); plt.close(fig)


def main()->None:
    _,a,_=strict_operator()
    out:dict[str,Any]={"metadata":{"programs":"ST106-ST117","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"projection_hypothesis":"counterfactual"}}
    out["ST106"]=st106_naturality_classification(); out["ST107"]=st107_lift_facial_geometry(a)
    out["ST108"]=st108_transcendental_fold(a); out["ST109"]=st109_maximal_recoverable_code()
    out["ST110"]=st110_equivariant_measure_nonselection(); out["ST111"]=st111_antisymmetric_memory()
    out["ST112"]=st112_optimal_design(a); out["ST113"]=st113_spin_refinement_category()
    out["ST114"]=st114_locality_algebra(); out["ST115"]=st115_independent_nuisance_certificate(a)
    out["ST116"]=st116_cgp_vs_thermal_operations(a); out["ST117"]=st117_adversarial_lifts(a)
    out["recommended_next_programs"]=[
        {"id":"ST118","priority":1,"study":"derive or obstruct a fine-blind observable algebra from an explicit operational equivalence relation"},
        {"id":"ST119","priority":2,"study":"enumerate bounded vertices/faces of the 78-dimensional lift spectrahedron with symmetry reduction"},
        {"id":"ST120","priority":3,"study":"export a dependency-minimal exact replay of the ST108 transcendental interval certificate"},
        {"id":"ST121","priority":4,"study":"classify all maximal recoverable twelve-dimensional codes up to sector unitaries"},
        {"id":"ST122","priority":5,"study":"test whether a strict-derived invariant state can induce a nonuniform branch measure without QW-2191 replay"},
        {"id":"ST123","priority":6,"study":"classify the minimal nonpassive port/source extension and its thermodynamic cost"},
        {"id":"ST124","priority":7,"study":"interval-certify the continuous beta optimum and robust finite-count design of ST112"},
        {"id":"ST125","priority":8,"study":"classify coherent twist systems under non-binary and mixed refinement maps"},
        {"id":"ST126","priority":9,"study":"replace vertex locality by causal-net or code-subspace locality and retest algebra selection"},
        {"id":"ST127","priority":10,"study":"extend ST115 to nonuniform and correlated detector contamination"},
        {"id":"ST128","priority":11,"study":"declare a finite bath/dilation class and decide CGP outside TO by exact semidefinite feasibility"},
        {"id":"ST129","priority":12,"study":"minimax-optimize fine observables against the ST117 adversarial lift cone"},
    ]
    out["central_verdict"]=(
        "The projection layer is not selected by symmetry or vertex locality: swap symmetry leaves M12 plus M12 and locality generates M24. "
        "The desired M12 plus C algebra requires an explicit fine-blindness premise.  The complete lift cone has twelve unbounded adversarial "
        "rays.  Nonetheless, the declared transcendental strict kernel now has a uniform interval fold certificate, pinching supports maximal "
        "twelve-dimensional recoverable codes, and the full independent synthetic nuisance box has a globally certified worst corner."
    )
    out["epistemic_boundary"]="No QW-2191 closure, strict physical projection, apparatus, dimensional scale, laboratory evidence, legacy-role transfer, Standard Model, gravity, L_total, or ToE closure is claimed."
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h); w.writerow(["program","object","status"])
        for k in range(106,118): w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out)
    print(json.dumps({"results":RESULTS.name,"programs":12,"figures":8},indent=2))


if __name__=="__main__": main()
