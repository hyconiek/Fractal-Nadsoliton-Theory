#!/usr/bin/env python3
"""FIN ST94--ST105: projection uniqueness, recovery, and falsification.

All physical projection claims remain counterfactual. Exact results concern the
declared finite operators, algebras, and synthetic probability models only.
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
from scipy.linalg import expm
from scipy.optimize import minimize_scalar, root
from scipy.special import comb

from fin_programs_497_506_next_research import iv_bounds
from fin_st01_st15_research import N, random_orthogonal_fixing_uniform, strict_operator
from fin_st28_st45_research import dyadic_lift, saturation_energy_gradient_hessians
from fin_st46_st57_research import carrier_probability_table
from fin_st58_st69_research import interval_left_product, interval_matvec
from fin_st70_st81_research import reflection_expansion


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST94_ST105_Results.json"
SUMMARY = ROOT / "FIN_ST94_ST105_Summary.csv"
CERT99 = ROOT / "FIN_ST99_Fold_Krawczyk_Certificate.json"
SPEC105 = ROOT / "FIN_ST105_Projection_Falsification_Protocol.json"
FIG_DIR = ROOT / "FIN_ST94_ST105_Figures"
SEED = 20260819


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plus_minus_isometries(n: int = N) -> tuple[np.ndarray, np.ndarray]:
    plus = np.zeros((2 * n, n))
    minus = np.zeros((2 * n, n))
    for x in range(n):
        plus[x, x] = plus[x + n, x] = 1 / math.sqrt(2)
        minus[x, x] = 1 / math.sqrt(2)
        minus[x + n, x] = -1 / math.sqrt(2)
    return plus, minus


def st94_unique_conditional_expectation(a: np.ndarray) -> dict:
    plus, minus = plus_minus_isometries()
    pc, pf = plus @ plus.T, minus @ minus.T

    def expectation(x: np.ndarray) -> np.ndarray:
        fine_scalar = np.trace(pf @ x @ pf) / N
        return pc @ x @ pc + fine_scalar * pf

    rng = np.random.default_rng(SEED + 94)
    residuals = []
    for _ in range(10):
        x = rng.normal(size=(2 * N, 2 * N)) + 1j * rng.normal(size=(2 * N, 2 * N))
        b1c = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        b2c = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        b1 = plus @ b1c @ plus.T + rng.normal() * pf
        b2 = plus @ b2c @ plus.T + rng.normal() * pf
        residuals.append({
            "unital": float(np.linalg.norm(expectation(np.eye(2 * N)) - np.eye(2 * N))),
            "idempotent": float(np.linalg.norm(expectation(expectation(x)) - expectation(x))),
            "trace_preserving": float(abs(np.trace(expectation(x)) - np.trace(x))),
            "bimodule": float(np.linalg.norm(expectation(b1 @ x @ b2) - b1 @ expectation(x) @ b2)),
        })
    return {
        "program": "ST94",
        "object": "Unique Trace-Preserving Conditional Expectation after an Observable Algebra Is Fixed",
        "observable_algebra": "M12 on the coarse sector direct-sum C times identity on the fine sector",
        "formula": "E_B(X)=Pc X Pc + Tr(Pf X Pf)/12 Pf",
        "maximum_live_residuals": {key: max(row[key] for row in residuals) for key in residuals[0]},
        "theorem": (
            "For a fixed unital finite-dimensional *-subalgebra B and normalized trace, the trace-preserving conditional "
            "expectation E_B is unique: the bimodule and trace identities make X-E_B(X) Hilbert-Schmidt orthogonal to B. "
            "For B=M12 direct-sum C the displayed formula follows."
        ),
        "projection_uniqueness_verdict": (
            "Operational sufficiency uniquely determines the expectation only after the observable algebra B is supplied. "
            "It does not derive or uniquely select B, the plus sector, or a physical apparatus from strict A."
        ),
        "status": "proven_relative_uniqueness_but_observable_algebra_remains_an_axiom",
        "boundary": "This is a uniqueness theorem conditional on B, not a strict physical projection theorem.",
    }


def st95_complete_lift_fiber(a: np.ndarray) -> dict:
    plus, minus = plus_minus_isometries()
    c = float(np.max(np.diag(a)) + 1.0)
    b_interior = c * np.eye(N)
    l_interior = plus @ a @ plus.T + minus @ b_interior @ minus.T
    offdiag = l_interior - np.diag(np.diag(l_interior))
    interior_margin = -float(np.max(offdiag[np.triu_indices(2 * N, 1)]))
    rows = []
    for q in [0.0, 0.2, 0.7, 1.5]:
        lifted = dyadic_lift(a, q)
        b = minus.T @ lifted @ minus
        lower = min(float(b[i, j] - a[i, j]) for i in range(N) for j in range(N) if i != j)
        upper = min(float(-a[i, j] - b[i, j]) for i in range(N) for j in range(N) if i != j)
        diag = min(float(b[i, i] - a[i, i]) for i in range(N))
        rows.append({"q": q, "minimum_B_eigenvalue": float(np.min(np.linalg.eigvalsh(b))), "lower_inequality_margin": lower, "upper_inequality_margin": upper, "diagonal_margin": diag})
    return {
        "program": "ST95",
        "object": "Complete Symmetric Graph-Laplacian Fiber over A12",
        "fine_block_parameter_space_dimension": N * (N + 1) // 2,
        "classification": (
            "Every real symmetric L with LP=PA is uniquely L=P A P* + F B F*. It is PSD iff B is PSD. "
            "It is a graph Laplacian with nonpositive off-diagonal entries iff B_ii>=A_ii and "
            "A_ij<=B_ij<=-A_ij for i!=j."
        ),
        "interior_witness_B": "c I with c=max_i A_ii+1",
        "interior_graph_offdiagonal_margin": interior_margin,
        "interior_intertwining_residual": float(np.linalg.norm(l_interior @ plus - plus @ a)),
        "dyadic_family_rows": rows,
        "theorem": (
            "The complete coarse-preserving graph-generator fiber is a 78-dimensional convex spectrahedral set with a "
            "nonempty interior. The earlier q-family is only a one-dimensional subfamily and cannot be inferred from coarse data."
        ),
        "status": "proven_complete_finite_lift_fiber_classification",
        "boundary": "The classification increases rather than removes hidden-model nonuniqueness and supplies no physical refinement law.",
    }


def matrix_log_density(rho: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh((rho + rho.conj().T) / 2)
    return vectors @ np.diag(np.log(np.maximum(values, 1e-300))) @ vectors.conj().T


def relative_entropy(rho: np.ndarray, sigma: np.ndarray) -> float:
    return float(np.real(np.trace(rho @ (matrix_log_density(rho) - matrix_log_density(sigma)))))


def st96_petz_recovery_frontier() -> dict:
    plus, minus = plus_minus_isometries()
    pc, pf = plus @ plus.T, minus @ minus.T
    pinching = lambda x: pc @ x @ pc + pf @ x @ pf
    sigma = np.eye(2 * N) / (2 * N)
    rng = np.random.default_rng(SEED + 96)
    rows = []
    for _ in range(16):
        m = rng.normal(size=(2 * N, 2 * N)) + 1j * rng.normal(size=(2 * N, 2 * N))
        rho = m @ m.conj().T + 0.1 * np.eye(2 * N)
        rho /= np.trace(rho)
        pirho = pinching(rho)
        pythagorean = relative_entropy(rho, sigma) - relative_entropy(pirho, sigma) - relative_entropy(rho, pirho)
        coherence_trace_norm = float(np.sum(np.linalg.svd(rho - pirho, compute_uv=False)))
        rows.append({
            "relative_entropy_loss": relative_entropy(rho, pirho),
            "pythagorean_residual": pythagorean,
            "trace_norm_lost_coherence": coherence_trace_norm,
            "pinsker_margin": relative_entropy(rho, pirho) - 0.5 * coherence_trace_norm**2,
        })
    return {
        "program": "ST96",
        "object": "Petz-Recovery Frontier for Coarse/Fine Pinching",
        "rows": rows,
        "maximum_pythagorean_residual": max(abs(row["pythagorean_residual"]) for row in rows),
        "minimum_pinsker_margin": min(row["pinsker_margin"] for row in rows),
        "theorem": (
            "For a full-rank reference sigma in the pinched algebra, D(rho||sigma)=D(rho||Pi rho)+D(Pi rho||sigma). "
            "The Petz map recovers Pi rho, and equality in data processing is exact iff rho=Pi rho, i.e. iff no "
            "coarse/fine coherence was discarded. Pinsker gives D(rho||Pi rho)>=||rho-Pi rho||_1^2/2."
        ),
        "status": "proven_exact_recovery_criterion_with_numerical_checks",
        "boundary": "Recoverability is relative to a supplied pinching and reference state; neither is a derived physical measurement.",
    }


def st97_state_dependent_selector() -> dict:
    mu = 0.4
    angular_gain = 0.7

    def evolve(radius: float, theta: float, dt: float = 0.002, steps: int = 12000) -> tuple[float, float]:
        r, t = radius, theta
        for _ in range(steps):
            r += dt * r * (mu - r * r)
            t += dt * (-angular_gain * math.sin(12 * t))
        return r, t

    rows = []
    for theta in np.linspace(-math.pi, math.pi, 25, endpoint=False):
        r, t = evolve(1e-3, float(theta))
        branch = int(round((t % (2 * math.pi)) / (math.pi / 6))) % 12
        rows.append({"initial_theta": float(theta), "final_radius": r, "final_theta_mod_2pi": float(t % (2 * math.pi)), "selected_branch": branch})
    zero = evolve(0.0, 0.123)
    return {
        "program": "ST97",
        "object": "Equivariant State-Dependent Twelve-Branch Selector Test",
        "flow": "r_dot=r(mu-r^2), theta_dot=-g sin(12 theta)",
        "rows": rows,
        "zero_state_remains_zero": zero[0] == 0.0,
        "branches_reached": sorted({row["selected_branch"] for row in rows}),
        "theorem": (
            "The constructed flow is C12-equivariant and has twelve stable angular branches. Every nonzero initial state "
            "outside the twelve unstable angular separatrices selects a branch through its phase, but the symmetric state "
            "r=0 is invariant and selects none. Reflection "
            "maps each selected branch to its partner. Thus state dependence relocates the selector into initial/boundary "
            "data and does not create a canonical orbit member."
        ),
        "status": "proven_constructed_state_dependent_selection_but_QW2191_remains_open",
        "boundary": "The flow, initial perturbation and angular gain are constructed and not sourced by strict FIN.",
    }


def st98_passive_realization_extension(a: np.ndarray) -> dict:
    eigenvalues = np.linalg.eigvalsh(a)[1:]
    poles = eigenvalues + 0.15
    rng = np.random.default_rng(SEED + 98)
    b = rng.normal(size=(len(poles), 4))
    residues = np.sum(b**2, axis=1)
    rows = []
    for omega in np.linspace(0, 5, 31):
        transfer = sum(residue / (1j * omega + pole) for residue, pole in zip(residues, poles))
        rows.append({"omega": float(omega), "real_transfer": float(np.real(transfer)), "absolute_transfer": float(abs(transfer))})
    return {
        "program": "ST98",
        "object": "Passive State-Space Extension of the Temporal-Gain No-Go",
        "minimum_pole": float(np.min(poles)),
        "minimum_residue": float(np.min(residues)),
        "minimum_sampled_real_transfer": min(row["real_transfer"] for row in rows),
        "rows": rows,
        "theorem": (
            "Every reciprocal finite realization M(z)=B*(zI+C)^(-1)B with C=C*>0 has poles -c_j<0 and nonnegative "
            "residue matrices b_j*b_j. Its scalar restrictions are positive real and passive. If C and B are obtained "
            "only from commuting real functions of strict A, no signed active residue or temporal orientation appears."
        ),
        "status": "proven_passive_state_space_no_go_in_extended_reciprocal_class",
        "boundary": "Indefinite metrics, antisymmetric couplings, pumps or nonequilibrium states evade the theorem only as new sourced objects.",
    }


def iv(value: float | str | tuple[float, float]):
    if isinstance(value, tuple):
        return mp.iv.mpf([str(value[0]), str(value[1])])
    if isinstance(value, (float, np.floating)):
        return mp.iv.mpf([
            str(float(np.nextafter(value, -np.inf))),
            str(float(np.nextafter(value, np.inf))),
        ])
    return mp.iv.mpf(str(value))


def st99_interval_fold_certificate(a: np.ndarray) -> dict:
    previous = json.loads((ROOT / "FIN_ST70_ST81_Results.json").read_text(encoding="utf-8"))["ST77"]
    expansion = reflection_expansion()
    selected = np.arange(7)
    x0 = np.r_[previous["fold_state_reduced"], previous["fold_kappa"], previous["fold_null_vector"]].astype(float)

    def system(candidate: np.ndarray) -> np.ndarray:
        q, kappa, vector = candidate[:7], candidate[7], candidate[8:]
        state = expansion @ q
        _, gradient, hessian, _ = saturation_energy_gradient_hessians(kappa * a, state)
        matrix = hessian[np.ix_(selected, np.arange(N))] @ expansion
        return np.r_[gradient[selected], matrix @ vector, 0.5 * (vector @ vector - 1.0)]

    x0 = root(system, x0, method="lm", options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14, "maxiter": 20000}).x
    mp.iv.dps = 70

    def interval_system_jacobian(radius: float):
        q = [iv((x0[i] - radius, x0[i] + radius)) for i in range(7)]
        kappa = iv((x0[7] - radius, x0[7] + radius))
        v = [iv((x0[8 + i] - radius, x0[8 + i] + radius)) for i in range(7)]
        u = []
        w = []
        for site in range(N):
            index = site if site <= 6 else N - site
            u.append(q[index]); w.append(v[index])
        aa = [[iv(a[i, j]) for j in range(N)] for i in range(N)]
        au = [sum((aa[i][j] * u[j] for j in range(N)), iv(0)) for i in range(N)]
        aw = [sum((aa[i][j] * w[j] for j in range(N)), iv(0)) for i in range(N)]
        h, rdiag, drdu = [], [], []
        for item in u:
            rho = item**2; den = 1 + rho / 2
            qfun = rho / den; qp = den**-2; qpp = -(den**-3); qppp = iv("1.5") * den**-4
            hh = -qfun * qp + iv("0.075")
            hp = -(qp**2 + qfun * qpp)
            hpp = -(3 * qp * qpp + qfun * qppp)
            h.append(hh)
            rdiag.append(2 * hh + 4 * rho * hp)
            drdu.append(2 * item * (6 * hp + 4 * rho * hpp))
        g1 = [kappa * au[i] + 2 * u[i] * h[i] for i in selected]
        g2 = [kappa * aw[i] + rdiag[i] * w[i] for i in selected]
        g3 = iv("0.5") * (sum((item**2 for item in v), iv(0)) - 1)
        g = g1 + g2 + [g3]
        j = [[iv(0) for _ in range(15)] for _ in range(15)]
        # H_red in the q and v columns.
        for i in range(7):
            for col in range(7):
                total = iv(0)
                for site in range(N):
                    if (site if site <= 6 else N - site) == col:
                        total += kappa * aa[i][site]
                if i == col:
                    total += rdiag[i]
                j[i][col] = total
                j[7 + i][8 + col] = total
            j[i][7] = au[i]
            j[7 + i][i] = drdu[i] * w[i]
            j[7 + i][7] = aw[i]
            j[14][8 + i] = v[i]
        glo = np.array([iv_bounds(item)[0] for item in g]); ghi = np.array([iv_bounds(item)[1] for item in g])
        jlo = np.array([[iv_bounds(item)[0] for item in row] for row in j]); jhi = np.array([[iv_bounds(item)[1] for item in row] for row in j])
        return glo, ghi, jlo, jhi

    # Point Jacobian from a tight interval midpoint, then Krawczyk boxes.
    g0lo, g0hi, j0lo, j0hi = interval_system_jacobian(0.0)
    j0 = 0.5 * (j0lo + j0hi)
    c = np.linalg.inv(j0)
    hcheck = 1e-6
    numerical_jacobian = np.column_stack([
        (system(x0 + hcheck * np.eye(15)[j]) - system(x0 - hcheck * np.eye(15)[j])) / (2 * hcheck)
        for j in range(15)
    ])
    analytic_jacobian_residual = float(np.linalg.norm(j0 - numerical_jacobian, ord=np.inf))
    attempts = []
    accepted = None
    for radius in [1e-8, 3e-9, 1e-9, 3e-10, 1e-10, 3e-11]:
        _, _, jlo, jhi = interval_system_jacobian(radius)
        cglo, cghi = interval_matvec(c, c, g0lo, g0hi)
        ylo = np.nextafter(x0 - cghi, -np.inf); yhi = np.nextafter(x0 - cglo, np.inf)
        cjlo, cjhi = interval_left_product(c, jlo, jhi)
        mlo, mhi = -cjhi, -cjlo
        for i in range(15):
            mlo[i, i] = np.nextafter(mlo[i, i] + 1.0, -np.inf)
            mhi[i, i] = np.nextafter(mhi[i, i] + 1.0, np.inf)
        dlo = np.full(15, -radius); dhi = np.full(15, radius)
        mdlo, mdhi = interval_matvec(mlo, mhi, dlo, dhi)
        klo = np.nextafter(ylo + mdlo, -np.inf); khi = np.nextafter(yhi + mdhi, np.inf)
        left_margin = klo - (x0 - radius); right_margin = (x0 + radius) - khi
        margin = float(min(np.min(left_margin), np.min(right_margin)))
        row = {"radius": radius, "minimum_strict_inclusion_margin": margin, "included": margin > 0, "maximum_Krawczyk_half_width": float(np.max((khi-klo)/2))}
        attempts.append(row)
        if margin > 0 and accepted is None:
            accepted = row
    certificate = {
        "arithmetic": "mpmath.iv 70-digit outward intervals plus nextafter-enclosed binary64 matrix products",
        "operator_scope": "the frozen binary64 strict A used by the local research scripts",
        "center": x0.tolist(),
        "center_residual_inf": float(np.linalg.norm(system(x0), ord=np.inf)),
        "analytic_vs_finite_difference_J_residual_inf": analytic_jacobian_residual,
        "attempts": attempts,
        "accepted": accepted,
        "proof_boundary": "A source-code interval certificate for the frozen binary64 operator; not a proof-assistant theorem and not a transcendental-parameter enclosure of A.",
    }
    CERT99.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST99",
        "object": "Krawczyk Certificate for the Reflection-Reduced Simple Fold",
        "certificate_file": CERT99.name,
        "certificate_sha256": sha(CERT99),
        **certificate,
        "status": "proven_source_code_interval_fold_for_frozen_binary64_operator" if accepted else "blocked_interval_fold_certificate_not_obtained",
        "boundary": certificate["proof_boundary"],
    }


def st100_spin_refinement_tower() -> dict:
    rows = []
    for h12 in [1, -1]:
        h24 = h12**2
        h48 = h24**2
        rows.append({"h12": h12, "pullback_h24": h24, "pullback_h48": h48})
    return {
        "program": "ST100",
        "object": "Z2 Spin-Holonomy Pullback through the 12-to-24-to-48 Tower",
        "rows": rows,
        "compatible_pullback_sequences": [[1, 1, 1], [-1, 1, 1]],
        "stationary_antiperiodic_sequence_exists": False,
        "theorem": (
            "A degree-two cyclic cover sends Z2 holonomy h to h^2=+1. Therefore an antiperiodic class at level 12 "
            "pulls back to the periodic class at level 24 and remains periodic at level 48. Maintaining antiperiodicity "
            "at every refinement requires a new twist choice at every level, not functorial inheritance."
        ),
        "status": "proven_spin_refinement_pullback_and_new_twist_obligation",
        "boundary": "The tower classifies boundary data but neither chooses twists nor derives physical spin/fermions.",
    }


def st101_minimal_fine_observable(a: np.ndarray) -> dict:
    plus, minus = plus_minus_isometries()
    beta = 0.8
    b0 = minus.T @ dyadic_lift(a, 0.0) @ minus
    zc = float(np.trace(expm(-beta * a)))
    zf0 = float(np.trace(expm(-beta * b0)))
    rows = []
    for q in [0.0, 0.2, 0.7, 1.5]:
        lifted = dyadic_lift(a, q)
        b = minus.T @ lifted @ minus
        shift_residual = float(np.linalg.norm(b - b0 - 2 * q * np.eye(N)))
        pcoarse = zc / (zc + math.exp(-2 * beta * q) * zf0)
        recovered = -math.log(((1-pcoarse)/pcoarse) * (zc/zf0)) / (2*beta)
        fisher = 4 * beta**2 * pcoarse * (1-pcoarse)
        rows.append({"q": q, "fine_scalar_shift_residual": shift_residual, "global_Gibbs_coarse_probability": pcoarse, "q_recovered_from_probability": recovered, "Bernoulli_Fisher_information": fisher})
    return {
        "program": "ST101",
        "object": "Minimal Fine-Observable Identifiability Theorem",
        "beta": beta,
        "rows": rows,
        "theorem": (
            "The dyadic fine block obeys B(q)=B(0)+2q I. One calibrated fine-sector generator expectation identifies q. "
            "Equivalently, one coarse-versus-fine binary outcome from the globally normalized Gibbs state has "
            "p_c(q)=Z_c/(Z_c+exp(-2 beta q)Z_f0), a strictly monotone statistic that identifies q when beta and baseline "
            "partition functions are supplied."
        ),
        "status": "proven_one_scalar_fine_access_identifies_q_conditionally",
        "boundary": "The measurement, beta, global normalization and baseline calibration are additional operational resources.",
    }


def st102_algebra_tower(a: np.ndarray) -> dict:
    p12, _ = plus_minus_isometries(12)
    p24, _ = plus_minus_isometries(24)
    composite = p24 @ p12
    a24 = dyadic_lift(a, 0.37)
    a48 = dyadic_lift(a24, 0.61)
    # On M48, average over translations by 24 and by the order-four shift 12.
    def shift(size: int, amount: int) -> np.ndarray:
        matrix = np.zeros((size, size))
        for x in range(size): matrix[(x + amount) % size, x] = 1
        return matrix
    t24, t12 = shift(48, 24), shift(48, 12)
    def e24(x): return 0.5 * (x + t24 @ x @ t24.T)
    def e12(x): return sum(np.linalg.matrix_power(t12, k) @ x @ np.linalg.matrix_power(t12.T, k) for k in range(4)) / 4
    rng = np.random.default_rng(SEED + 102)
    commute = []
    for _ in range(8):
        x = rng.normal(size=(48, 48))
        commute.append(float(np.linalg.norm(e24(e12(x)) - e12(e24(x)))))
    return {
        "program": "ST102",
        "object": "Finite C*-Algebra Refinement Tower and Commuting-Square Test",
        "composite_isometry_residual": float(np.linalg.norm(composite.T @ composite - np.eye(N))),
        "generator_tower_residual": float(np.linalg.norm(a48 @ composite - composite @ a)),
        "maximum_conditional_expectation_commutator_residual": max(commute),
        "theorem": (
            "The half-translation and quarter-translation group averages are trace-preserving conditional expectations. "
            "Their cyclic subgroups commute and are nested, so the expectations form a finite commuting-square tower. "
            "The composite coarse isometry intertwines A48 with A12 for every intermediate fine parameter."
        ),
        "status": "proven_finite_commuting_expectation_tower",
        "boundary": "Tower consistency is automatic group averaging and does not select fine blocks or a physical RG law.",
    }


def st103_nuisance_certificate(a: np.ndarray) -> dict:
    rng = np.random.default_rng(20260818 + 87)  # reproduces ST87 carrier seed
    qmat = random_orthogonal_fixing_uniform(rng)
    p = carrier_probability_table(a, np.eye(N), transported=False) / N
    q = carrier_probability_table(a, qmat, transported=False) / N
    u = np.full_like(p, 1 / p.size)

    def common_information(epsilon: float) -> tuple[float, float]:
        pp = (1-epsilon)*p + epsilon*u; qq = (1-epsilon)*q + epsilon*u
        opt = minimize_scalar(lambda s: float(np.sum(pp**s * qq**(1-s))), bounds=(0,1), method="bounded")
        return -math.log(opt.fun), float(opt.x)
    rows = [{"epsilon": float(e), "Chernoff_information": common_information(float(e))[0], "s": common_information(float(e))[1]} for e in np.linspace(0,0.12,25)]
    return {
        "program": "ST103",
        "object": "Certified Common-Contamination Subproblem and Independent-Nuisance Boundary",
        "rows": rows,
        "common_endpoint_information": rows[-1]["Chernoff_information"],
        "common_monotone": all(rows[i+1]["Chernoff_information"] <= rows[i]["Chernoff_information"] + 1e-12 for i in range(len(rows)-1)),
        "theorem": (
            "When the same contamination epsilon is applied to both hypotheses, increasing epsilon is composition with "
            "a common stochastic channel. Chernoff information cannot increase by data processing, so the rigorous worst "
            "case on the common diagonal is epsilon=0.12."
        ),
        "independent_two_parameter_interval_certificate_obtained": False,
        "failure_reason": "The full independent (epsilon_P,epsilon_Q) box is not a single common channel and no outward global max-min certificate was obtained.",
        "status": "proven_common_nuisance_endpoint_but_full_ST103_goal_open",
        "boundary": "ST87's independent-box optimum remains strong numerical evidence, not interval global optimization.",
    }


def st104_thermal_hierarchy(a: np.ndarray) -> dict:
    eigenvalues = np.linalg.eigvalsh(a)
    groups = []
    for value in eigenvalues:
        if not groups or abs(value-groups[-1][0]) > 1e-9: groups.append([float(value),1])
        else: groups[-1][1] += 1
    multiplicities = [row[1] for row in groups]
    zero_mode_dimension = sum(m*m for m in multiplicities)
    return {
        "program": "ST104",
        "object": "Thermal-Operation Hierarchy with Strict Spectral Degeneracies",
        "energy_multiplicities": multiplicities,
        "zero_Bohr_frequency_operator_dimension": zero_mode_dimension,
        "hierarchy": "thermal operations subset of covariant Gibbs-preserving maps, strict subset of all Gibbs-preserving maps",
        "right_inclusion_strict_witness": "ST89 measure-and-prepare channel creates energy coherence from a stationary input",
        "left_inclusion_strictness_for_this_spectrum": "open_without_a_declared_bath_and_dilation_class",
        "theorem": (
            "Energy-conserving thermal operations preserve gamma and are time-translation covariant. ST89 proves the "
            "covariant Gibbs-preserving class is strictly smaller than all Gibbs-preserving maps. Strict degeneracies "
            "produce a 22-dimensional zero-frequency operator sector in which covariance permits within-block coherence."
        ),
        "status": "proven_hierarchy_and_right_strictness_left_strictness_open",
        "boundary": "FIN supplies no physical Hamiltonian scale, bath spectrum, energy-conserving dilation, or thermal instrument.",
    }


def binomial_bayes_error(n: int, p0: float, p1: float) -> float:
    tv = 0.0
    for k in range(n+1):
        a = comb(n,k) * p0**k * (1-p0)**(n-k)
        b = comb(n,k) * p1**k * (1-p1)**(n-k)
        tv += abs(a-b)
    return 0.5 * (1 - 0.5 * tv)


def st105_projection_falsification(a: np.ndarray) -> dict:
    beta, q0, q1 = 0.8, 0.2, 0.7
    plus, minus = plus_minus_isometries()
    b0 = minus.T @ dyadic_lift(a, 0) @ minus
    zc = float(np.trace(expm(-beta*a))); zf = float(np.trace(expm(-beta*b0)))
    probability = lambda q: zc / (zc + math.exp(-2*beta*q)*zf)
    p0, p1 = probability(q0), probability(q1)
    rows = [{"events": n, "optimal_equal_prior_Bayes_error": binomial_bayes_error(n,p0,p1)} for n in range(1,401)]
    threshold = next((row["events"] for row in rows if row["optimal_equal_prior_Bayes_error"] <= 0.01), None)
    packet = {
        "hypotheses": {"H0_q": q0, "H1_q": q1, "beta": beta},
        "coarse_conditioned_channel_TV": 0.0,
        "fine_binary_probabilities": [p0,p1],
        "first_event_count_with_Bayes_error_at_most_0_01": threshold,
        "rows": rows,
        "decision_rule": "likelihood ratio for the number of coarse-sector outcomes in n Bernoulli events",
        "physical_boundary": "entirely synthetic; assumes global Gibbs preparation, calibrated beta, coarse/fine POVM and independent raw counts",
    }
    SPEC105.write_text(json.dumps(packet, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST105",
        "object": "Projection-Specific Synthetic Falsification Protocol",
        "protocol_file": SPEC105.name,
        "protocol_sha256": sha(SPEC105),
        **packet,
        "theorem": (
            "Every test confined to conditioned coarse functional-calculus records has equal distributions under q0 and q1 "
            "and equal-prior error 1/2. A single declared coarse-versus-fine Gibbs occupation observable separates them; "
            "the exact binomial likelihood gives the displayed finite-event error curve."
        ),
        "status": "proven_coarse_impossibility_with_strong_binary64_fine_protocol_performance",
        "boundary": packet["physical_boundary"],
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(7,4)); st95=results["ST95"]; ax.bar(["dyadic q-family","complete B fiber"],[1,st95["fine_block_parameter_space_dimension"]]); ax.set(ylabel="parameter-space dimension",title="ST95 hidden lift freedom"); fig.tight_layout(); fig.savefig(FIG_DIR/"st95_lift_dimension.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st96=results["ST96"]; ax.scatter([r["trace_norm_lost_coherence"] for r in st96["rows"]],[r["relative_entropy_loss"] for r in st96["rows"]]); x=np.linspace(0,max(r["trace_norm_lost_coherence"] for r in st96["rows"]),100); ax.plot(x,.5*x*x,"--",label="Pinsker lower bound"); ax.set(xlabel="lost coherence trace norm",ylabel="relative-entropy loss",title="ST96 exact recovery frontier"); ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR/"st96_recovery.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st97=results["ST97"]; ax.scatter([r["initial_theta"] for r in st97["rows"]],[r["selected_branch"] for r in st97["rows"]]); ax.set(xlabel="initial phase",ylabel="selected branch",title="ST97 state dependence relocates selection"); fig.tight_layout(); fig.savefig(FIG_DIR/"st97_state_selector.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st98=results["ST98"]; ax.plot([r["omega"] for r in st98["rows"]],[r["real_transfer"] for r in st98["rows"]]); ax.axhline(0,color="black",lw=.8); ax.set(xlabel="frequency",ylabel="Re M(i omega)",title="ST98 passive positive-real response"); fig.tight_layout(); fig.savefig(FIG_DIR/"st98_passive_response.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st99=results["ST99"]; ax.semilogx([r["radius"] for r in st99["attempts"]],[r["minimum_strict_inclusion_margin"] for r in st99["attempts"]],"o-"); ax.axhline(0,color="black",lw=.8); ax.set(xlabel="box radius",ylabel="Krawczyk inclusion margin",title="ST99 interval fold certificate"); fig.tight_layout(); fig.savefig(FIG_DIR/"st99_krawczyk.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st101=results["ST101"]; ax.plot([r["q"] for r in st101["rows"]],[r["global_Gibbs_coarse_probability"] for r in st101["rows"]],"o-"); ax.set(xlabel="hidden q",ylabel="coarse-sector probability",title="ST101 one fine-access statistic identifies q"); fig.tight_layout(); fig.savefig(FIG_DIR/"st101_identifiability.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st103=results["ST103"]; ax.plot([r["epsilon"] for r in st103["rows"]],[r["Chernoff_information"] for r in st103["rows"]],"o-"); ax.set(xlabel="common contamination",ylabel="Chernoff information",title="ST103 certified common-channel monotonicity"); fig.tight_layout(); fig.savefig(FIG_DIR/"st103_common_nuisance.png",dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,4)); st105=results["ST105"]; ax.semilogy([r["events"] for r in st105["rows"]],[r["optimal_equal_prior_Bayes_error"] for r in st105["rows"]]); ax.axhline(.01,color="red",ls="--"); ax.set(xlabel="synthetic binary events",ylabel="optimal Bayes error",title="ST105 fine-access falsification protocol"); fig.tight_layout(); fig.savefig(FIG_DIR/"st105_falsification.png",dpi=190); plt.close(fig)


def write_summary(results: dict) -> None:
    with SUMMARY.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.writer(handle); writer.writerow(["program","object","status"])
        for k in range(94,106):
            row=results[f"ST{k}"]; writer.writerow([row["program"],row["object"],row["status"]])


def main() -> None:
    _, a, _ = strict_operator()
    results: dict[str,Any] = {"metadata":{"programs":"ST94-ST105","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"hypothesis":"H_PROJ remains counterfactual"}}
    results["ST94"]=st94_unique_conditional_expectation(a)
    results["ST95"]=st95_complete_lift_fiber(a)
    results["ST96"]=st96_petz_recovery_frontier()
    results["ST97"]=st97_state_dependent_selector()
    results["ST98"]=st98_passive_realization_extension(a)
    results["ST99"]=st99_interval_fold_certificate(a)
    results["ST100"]=st100_spin_refinement_tower()
    results["ST101"]=st101_minimal_fine_observable(a)
    results["ST102"]=st102_algebra_tower(a)
    results["ST103"]=st103_nuisance_certificate(a)
    results["ST104"]=st104_thermal_hierarchy(a)
    results["ST105"]=st105_projection_falsification(a)
    results["recommended_next_programs"]=[
        {"id":"ST106","priority":1,"study":"classify which additional naturality axioms can select the observable algebra without assuming it"},
        {"id":"ST107","priority":2,"study":"derive extreme rays and facial dimensions of the 78-dimensional graph-lift spectrahedron"},
        {"id":"ST108","priority":3,"study":"interval-enclose the strict transcendental kernel entries and lift the ST99 certificate beyond frozen binary64 A"},
        {"id":"ST109","priority":4,"study":"test Petz sufficiency for multiple noncommuting reference states and find the maximal recoverable code"},
        {"id":"ST110","priority":5,"study":"seek a strict-derived nonzero initial-state measure or prove equivariant-measure nonselection"},
        {"id":"ST111","priority":6,"study":"classify active memory realizations with one antisymmetric state-space block and audit their source requirements"},
        {"id":"ST112","priority":7,"study":"derive optimal Fisher/Chernoff design for q under finite beta and detector nuisance"},
        {"id":"ST113","priority":8,"study":"classify spin-twist cocycles on the full refinement groupoid rather than individual covers"},
        {"id":"ST114","priority":9,"study":"seek a physically interpretable observable algebra from locality, not spectral convenience"},
        {"id":"ST115","priority":10,"study":"obtain a full independent-box interval max-min certificate left open by ST103"},
        {"id":"ST116","priority":11,"study":"construct or obstruct a covariant Gibbs-preserving channel outside thermal operations for the strict degeneracy pattern"},
        {"id":"ST117","priority":12,"study":"build adversarial projection alternatives outside the dyadic family using ST95 extreme points"},
    ]
    results["central_verdict"]=(
        "Projection becomes unique only after an observable algebra is supplied; the strict core does not select that algebra. "
        "The complete hidden lift fiber is 78-dimensional, so ST86 nonidentifiability is far larger than the q-family. "
        "Lost coherences are exactly characterized by Petz sufficiency, and one deliberately fine/global observable identifies q. "
        "The projection hypothesis is therefore mathematically sharp and falsifiable only after extra operational access, but it remains unselected and untested physics."
    )
    results["epistemic_boundary"]="No QW-2191 closure, physical projection theorem, legacy-to-strict completion or role transfer, dimensional source, laboratory evidence, Standard Model, gravity, L_total or ToE closure is claimed."
    make_figures(results); write_summary(results)
    RESULTS.write_text(json.dumps(native(results),indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"results":RESULTS.name,"programs":12,"figures":8},indent=2))


if __name__=="__main__":
    main()
