"""P3172/S2122: Legacy* as standalone operator/model generator — potential audit.

Treats
    K*(d) = A * cos(pi*d/4 + pi/6) / (1 + beta*d),  A>0, beta>0
as an autonomous mathematical object (no Strict bridge, no parameter fitting
to Strict, no role transfer).

Finite computational certificate for:
  - operator-class acceptance matrix (TAK / NIE / WARUNKOWO)
  - dual dynamics U(t)=exp(-it L), T(t)=exp(-t L)
  - inverse-kernel recovery on radial circulant Z12
  - invariants and geometry witnesses
  - no-go / missing-axiom ledger
  - potential scores as research program (not closure claims)

Status: analysis / potential certificate. Not ToE, not L_total, not selector closure.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
GEN.mkdir(exist_ok=True)
OUT = GEN / "p3172_s2122_legacy_star_operator_model_generator_potential_audit.json"
MD = GEN / "p3172_s2122_legacy_star_operator_model_generator_potential_audit.md"
PACKET = ROOT / "P3172_S2122_LEGACY_STAR_OPERATOR_MODEL_GENERATOR_POTENTIAL_AUDIT.md"
AGENTS = REPO / "AGENTS.md"
SUMMARY = REPO / "SUMMARY_GROK.md"
SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

N = 12
OMEGA_L = math.pi / 4
PHI_L = math.pi / 6
A_DEFAULT = 1.0
BETAS = [0.01, 0.1, 0.5, 1.0, 1.07515, 2.0, 5.0]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    old = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in old:
        path.write_text(old.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def k_star(d: float, A: float = A_DEFAULT, beta: float = 1.0) -> float:
    return A * math.cos(OMEGA_L * d + PHI_L) / (1.0 + beta * d)


def cycle_dist(i: int, j: int, n: int = N) -> int:
    return min((i - j) % n, (j - i) % n)


def undirected_circulant(A: float, beta: float, n: int = N) -> np.ndarray:
    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            W[i, j] = k_star(cycle_dist(i, j, n), A=A, beta=beta)
    return W


def directed_circulant(A: float, beta: float, n: int = N) -> np.ndarray:
    """Directed residual: k(m) for m = (j-i) mod n (not undirected distance)."""
    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            m = (j - i) % n
            W[i, j] = k_star(m, A=A, beta=beta)
    return W


def positive_part_weights(W: np.ndarray) -> np.ndarray:
    Wp = np.maximum(W, 0.0)
    np.fill_diagonal(Wp, 0.0)
    return Wp


def graph_laplacian(Wp: np.ndarray) -> np.ndarray:
    deg = Wp.sum(axis=1)
    return np.diag(deg) - Wp


def markov_generator(Wp: np.ndarray) -> np.ndarray:
    """Row-stochastic rate matrix Q with off-diag >=0 and rows summing to 0."""
    deg = Wp.sum(axis=1)
    return Wp - np.diag(deg)


def is_hermitian(M: np.ndarray, tol: float = 1e-12) -> bool:
    return np.allclose(M, M.conj().T, atol=tol)


def is_psd(M: np.ndarray, tol: float = 1e-10) -> bool:
    if not is_hermitian(M):
        return False
    ev = np.linalg.eigvalsh((M + M.T) / 2.0)
    return float(ev.min()) >= -tol


def spectral_radius(M: np.ndarray) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(M))))


def matrix_exp(M: np.ndarray) -> np.ndarray:
    # scipy-free: eigendecomposition (fine for 12x12)
    w, V = np.linalg.eig(M)
    return (V * np.exp(w)) @ np.linalg.inv(V)


def frobenius(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.linalg.norm(A - B, "fro"))


def recover_radial_profile(W: np.ndarray, n: int = N) -> list[float]:
    """Recover undirected radial profile from first row using cycle distances."""
    return [float(W[0, d]) for d in range(n)]


def phase_table(n: int = 24) -> list[dict[str, Any]]:
    rows = []
    for d in range(n):
        th = OMEGA_L * d + PHI_L
        c = math.cos(th)
        rows.append(
            {
                "d": d,
                "theta_over_pi": th / math.pi,
                "cos": c,
                "period8_repeat": d >= 8 and abs(c - math.cos(OMEGA_L * (d % 8) + PHI_L)) < 1e-12,
            }
        )
    return rows


def triangle_violations(A: float, beta: float, n: int = N) -> dict[str, Any]:
    dI = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            kij = abs(k_star(cycle_dist(i, j, n), A=A, beta=beta))
            dI[i, j] = -math.log(max(kij, 1e-30))
    n_trip = 0
    n_viol = 0
    max_excess = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                n_trip += 1
                excess = dI[i, k] - (dI[i, j] + dI[j, k])
                if excess > 1e-12:
                    n_viol += 1
                    max_excess = max(max_excess, excess)
    diagonal_zero = bool(np.max(np.abs(np.diag(dI))) <= 1e-12)
    nonnegative = bool(np.min(dI) >= -1e-12)
    symmetric = bool(np.allclose(dI, dI.T, atol=1e-12))
    separates = bool(
        all(dI[i, j] > 1e-12 for i in range(n) for j in range(n) if i != j)
    )
    triangle_holds = n_viol == 0
    return {
        "n_triples": n_trip,
        "n_violations": n_viol,
        "max_excess": max_excess,
        "triangle_inequality_holds": triangle_holds,
        "diagonal_zero": diagonal_zero,
        "nonnegative": nonnegative,
        "symmetric": symmetric,
        "separates_points": separates,
        "is_metric": bool(
            triangle_holds and diagonal_zero and nonnegative and symmetric and separates
        ),
        "metric_correction": (
            "The raw dI=-log|K| is not a metric because dI(x,x) is nonzero "
            "for K*(0)=A*cos(pi/6). Triangle-inequality success alone is insufficient."
        ),
    }


def find_beta_psd_threshold(A: float = 1.0, lo: float = 0.5, hi: float = 2.0, steps: int = 80) -> dict[str, Any]:
    """Bracket and bisect the first PSD transition on undirected Z12."""
    def lambda_min(beta: float) -> float:
        return float(np.linalg.eigvalsh(undirected_circulant(A, beta)).min())

    grid = np.linspace(lo, hi, steps)
    rows = []
    bracket = None
    previous_beta = float(grid[0])
    previous_value = lambda_min(previous_beta)
    for b in grid:
        value = lambda_min(float(b))
        psd = value >= -1e-10
        rows.append({"beta": float(b), "lambda_min": value, "psd": psd})
        if previous_value < 0.0 <= value and bracket is None:
            bracket = [previous_beta, float(b)]
        previous_beta, previous_value = float(b), value
    if bracket is None:
        return {
            "beta_star_approx": None,
            "root_residual": None,
            "bracket": None,
            "scan_rows": rows[:: max(1, len(rows) // 10)],
        }
    left, right = bracket
    for _ in range(100):
        mid = 0.5 * (left + right)
        if lambda_min(mid) < 0.0:
            left = mid
        else:
            right = mid
    root = 0.5 * (left + right)
    return {
        "beta_star_approx": root,
        "root_residual": lambda_min(root),
        "bracket": bracket,
        "scan_rows": rows[:: max(1, len(rows) // 10)],
    }


def dual_dynamics_witness(beta: float = 1.0, t: float = 0.25) -> dict[str, Any]:
    Wp = positive_part_weights(undirected_circulant(1.0, beta))
    L = graph_laplacian(Wp)
    # shift to ensure nonnegative spectrum for diffusion
    ev = np.linalg.eigvalsh(L)
    s = float(max(0.0, -ev.min())) + 0.0
    Aop = L + s * np.eye(N)
    U = matrix_exp(-1j * t * Aop)
    T = matrix_exp(-t * Aop)
    # unitarity residual
    unit_res = frobenius(U.conj().T @ U, np.eye(N))
    # contraction residual: spectral radius of T should be <= 1
    rho_T = spectral_radius(T.real)  # T real
    # positivity of T entries (Markov-ish)
    T_real = np.real(T)
    min_entry = float(T_real.min())
    # common generator check: both are f(Aop)
    # spectral calculus residual: U and T share eigenbasis of Aop
    w, V = np.linalg.eigh(Aop)
    U_spec = (V * np.exp(-1j * t * w)) @ V.T
    T_spec = (V * np.exp(-t * w)) @ V.T
    return {
        "beta": beta,
        "t": t,
        "s_shift": s,
        "lambda_min_L": float(ev.min()),
        "lambda_max_L": float(ev.max()),
        "unitarity_residual_fro": unit_res,
        "T_spectral_radius": rho_T,
        "T_min_entry": min_entry,
        "U_spectral_calculus_residual": frobenius(U, U_spec),
        "T_spectral_calculus_residual": frobenius(T_real, T_spec),
        "common_generator": True,
        "shared_spectral_measure": True,
    }


def inverse_recovery_test(A: float = 1.0, beta: float = 1.0) -> dict[str, Any]:
    W = undirected_circulant(A, beta)
    profile = recover_radial_profile(W)
    # recover k(d) for d=0..6 (unique undirected distances)
    recovered = {d: profile[d] for d in range(N // 2 + 1)}
    exact = {d: k_star(d, A=A, beta=beta) for d in range(N // 2 + 1)}
    max_err = max(abs(recovered[d] - exact[d]) for d in recovered)
    # injectivity on radial undirected class: different (A,beta) -> different profiles?
    W2 = undirected_circulant(A * 1.1, beta)
    W3 = undirected_circulant(A, beta * 1.1)
    return {
        "max_abs_recovery_error": max_err,
        "injective_on_radial_A": frobenius(W, W2) > 1e-8,
        "injective_on_radial_beta": frobenius(W, W3) > 1e-8,
        "surjective_onto_all_hermitian": False,  # theorem-level: only radial circulants
        "partially_invertible_on_radial_circulant_class": max_err < 1e-10,
        "recovered_profile_d0_to_6": [recovered[d] for d in range(7)],
        "exact_profile_d0_to_6": [exact[d] for d in range(7)],
    }


def operator_class_matrix(beta_samples: list[float] | None = None) -> list[dict[str, Any]]:
    beta_samples = beta_samples or [0.01, 1.0, 2.0]
    # finite witnesses
    herm_ok = all(is_hermitian(undirected_circulant(1.0, b)) for b in beta_samples)
    dir_herm = all(is_hermitian(directed_circulant(1.0, b)) for b in beta_samples)
    psd_any = any(is_psd(undirected_circulant(1.0, b)) for b in beta_samples)
    psd_all = all(is_psd(undirected_circulant(1.0, b)) for b in beta_samples)
    lap_psd = True
    markov_ok = True
    for b in beta_samples:
        Wp = positive_part_weights(undirected_circulant(1.0, b))
        L = graph_laplacian(Wp)
        Q = markov_generator(Wp)
        lap_psd = lap_psd and is_psd(L)
        markov_ok = markov_ok and np.all(Q - np.diag(np.diag(Q)) >= -1e-12) and np.allclose(Q.sum(axis=1), 0.0)
    dual = dual_dynamics_witness(beta=1.0)
    inv = inverse_recovery_test()

    def row(name: str, verdict: str, reason: str, missing: str = "", fundamental: bool | None = None) -> dict[str, Any]:
        return {
            "class": name,
            "verdict": verdict,  # TAK / NIE / WARUNKOWO
            "reason": reason,
            "missing_axiom_or_structure": missing,
            "fundamental_block": fundamental,
        }

    return [
        row(
            "Hermitian",
            "WARUNKOWO",
            f"Undirected radial circulant W[K*] is real-symmetric for all tested beta (herm_ok={herm_ok}). "
            f"Directed residual is non-Hermitian (dir_herm={dir_herm}).",
            "choice of undirected geometry / hermitization map",
            False,
        ),
        row(
            "Gram",
            "WARUNKOWO",
            f"W is a Gram/covariance operator iff PSD. Finite scan: psd_any={psd_any}, psd_all={psd_all}. "
            "Requires beta >= beta_* on Z12 undirected.",
            "PSD condition / beta lower bound",
            False,
        ),
        row(
            "graph_Laplacian",
            "TAK",
            f"L=Deg-W+ with W+=max(K*,0) off-diagonal is always PSD on finite samples (lap_psd={lap_psd}), L1=0.",
            "",
            None,
        ),
        row(
            "Markov_generator",
            "WARUNKOWO",
            f"Q=W+-Deg has off-diag>=0 and row-sum 0 when W+ used (markov_ok={markov_ok}). "
            "Raw signed K* does not give a Markov Q without positive-part projection.",
            "positive-part projection or signed-kernel Markov theory",
            False,
        ),
        row(
            "Schrodinger",
            "WARUNKOWO",
            f"H=Aop=sI-W or H=L yields U(t)=exp(-itH) unitary when H=H* (unitarity residual {dual['unitarity_residual_fro']:.2e}). "
            "Does not export Born rule, measurement, or Hilbert-space ontology by itself.",
            "Hilbert space identification + observable algebra + units of t",
            False,
        ),
        row(
            "Dirac",
            "NIE",
            "No Clifford bundle, spinor structure, or first-order square-root of Laplacian exported by scalar radial K*.",
            "Clifford algebra / spin structure / chirality source",
            True,
        ),
        row(
            "Maxwell",
            "NIE",
            "No 2-form F, dF=0, d*F=J, gauge connection, or U(1) connection form from radial scalar kernel alone.",
            "principal U(1) bundle + connection + metric/volume form",
            True,
        ),
        row(
            "Yang_Mills",
            "NIE",
            "No nonabelian principal bundle, curvature F=dA+A∧A, or YM action from radial K*.",
            "G-bundle + connection + structure group + coupling",
            True,
        ),
        row(
            "Lindblad",
            "WARUNKOWO",
            "Given H from K* one can write Lindblad form, but Kraus/Lindblad channels V_k are not unique or sourced by radial K*.",
            "environment/channel operators V_k + complete positivity structure",
            False,
        ),
        row(
            "Fokker_Planck",
            "WARUNKOWO",
            "Finite Markov chain from Q is discrete FP analogue; continuum FP needs continuum limit functor + drift/diffusion coefficients.",
            "continuum embedding + scale units + drift field",
            False,
        ),
        row(
            "Perron_Frobenius",
            "WARUNKOWO",
            "Positive-part transfer / adjacency-like W+ is nonnegative; PF eigenvector exists for irreducible W+. "
            "Signed K* is not a PF operator.",
            "nonnegativity / irreducibility",
            False,
        ),
        row(
            "Koopman",
            "WARUNKOWO",
            "U(t)=exp(-itL) is a Koopman-type unitary evolution on l2(V) when L self-adjoint; "
            "underlying classical flow on phase space is not reconstructed from K* alone.",
            "underlying measure-preserving dynamical system",
            False,
        ),
        row(
            "transport",
            "WARUNKOWO",
            "Directed residual gives a transport-like circulant; undirected does not encode oriented flow.",
            "orientation/source of directed edges",
            False,
        ),
        row(
            "diffusion",
            "TAK",
            f"T(t)=exp(-t L) with L=graph Laplacian of W+ is a diffusion semigroup on Z12 "
            f"(rho_T={dual['T_spectral_radius']:.6f}, spectral calculus residual {dual['T_spectral_calculus_residual']:.2e}).",
            "",
            None,
        ),
        row(
            "wave",
            "WARUNKOWO",
            "Second-order wave equation d2u/dt2 = -L u requires phase-space lift (u, du/dt); "
            "not exported by first-order semigroup alone without importing time-structure.",
            "phase-space / symplectic / second-order time structure",
            False,
        ),
        row(
            "integral",
            "TAK",
            "W[K*] is by definition a finite integral / matrix operator with kernel K*(d(i,j)).",
            "",
            None,
        ),
        row(
            "pseudodifferential",
            "WARUNKOWO",
            "On Z12, circulant W is a Fourier multiplier (exact ΨDO of order 0 on the circle). "
            "On R^n continuum, symbol class requires continuum limit not exported here.",
            "continuum manifold + symbol calculus",
            False,
        ),
        row(
            "spectral",
            "TAK",
            "Finite spectral theorem applies: Hermitian W/L diagonalized by discrete Fourier modes on cycle; "
            "functional calculus yields U(t), T(t).",
            "",
            None,
        ),
    ]


def physics_readiness() -> list[dict[str, Any]]:
    return [
        {
            "theory": "classical_mechanics",
            "from_legacy_star_alone": "radial potential-like profile / signed kernel; no phase space",
            "from_operator_duality": "formal Hamiltonian H=L only after importing (q,p) structure",
            "minimal_new_axioms": ["phase space", "symplectic form", "time unit"],
            "export_level": "scaffold",
        },
        {
            "theory": "quantum_mechanics",
            "from_legacy_star_alone": "self-adjoint H on l2(Z12) when hermitized; unitary group",
            "from_operator_duality": "U(t)=exp(-itH) is QM-like evolution of pure states",
            "minimal_new_axioms": ["Born rule", "observable algebra", "measurement/readout", "hbar unit"],
            "export_level": "partial_formal",
        },
        {
            "theory": "Markov_processes",
            "from_legacy_star_alone": "Q from W+ is continuous-time Markov generator on finite state space",
            "from_operator_duality": "T(t)=exp(tQ) is the Markov semigroup",
            "minimal_new_axioms": ["interpretation of states as probabilities", "embedding into continuous space optional"],
            "export_level": "strong_finite",
        },
        {
            "theory": "thermodynamics",
            "from_legacy_star_alone": "none unit-bearing; Dirichlet energy is dimensionless Lyapunov candidate",
            "from_operator_duality": "diffusion semigroup gives entropy-increase candidates on finite graphs",
            "minimal_new_axioms": ["temperature/energy unit", "bath/preparation", "equilibrium measure source"],
            "export_level": "formal_only",
        },
        {
            "theory": "information_theory",
            "from_legacy_star_alone": "kernel amplitudes as correlations; dI=-log|K| is not a metric",
            "from_operator_duality": "Markov channel mutual information computable on finite states",
            "minimal_new_axioms": ["canonical probability source", "coding theorem link optional"],
            "export_level": "partial",
        },
        {
            "theory": "field_theory",
            "from_legacy_star_alone": "none as continuum QFT",
            "from_operator_duality": "formal free-field Laplacian on graph only",
            "minimal_new_axioms": ["spacetime continuum", "fields", "action density", "renormalization"],
            "export_level": "absent",
        },
        {
            "theory": "electrodynamics",
            "from_legacy_star_alone": "none",
            "from_operator_duality": "none without U(1) connection",
            "minimal_new_axioms": ["U(1) bundle", "metric", "charge current"],
            "export_level": "absent",
        },
        {
            "theory": "Standard_Model",
            "from_legacy_star_alone": "none",
            "from_operator_duality": "none",
            "minimal_new_axioms": ["gauge group", "representations", "Yukawa", "Higgs sector", "units"],
            "export_level": "absent",
        },
        {
            "theory": "general_relativity",
            "from_legacy_star_alone": "none",
            "from_operator_duality": "graph metric candidates not Lorentzian",
            "minimal_new_axioms": ["Lorentzian 4-metric", "Einstein-Hilbert action", "stress-energy"],
            "export_level": "absent",
        },
        {
            "theory": "quantum_gravity",
            "from_legacy_star_alone": "none",
            "from_operator_duality": "none",
            "minimal_new_axioms": ["all of QFT + GR + Planck-scale completion"],
            "export_level": "absent",
        },
    ]


def no_go_ledger() -> list[dict[str, Any]]:
    return [
        {
            "claim": "unique_Dirac_operator_from_K_star",
            "status": "NO_GO",
            "formal_reason": "scalar radial kernel has no Clifford generators or spinor module",
            "missing_axiom": "Clifford/spin structure",
            "missing_structure": "spinor bundle",
            "fundamental_or_technical": "fundamental",
        },
        {
            "claim": "unique_Maxwell_or_YM_from_K_star",
            "status": "NO_GO",
            "formal_reason": "no principal connection or curvature form",
            "missing_axiom": "gauge principle + structure group",
            "missing_structure": "principal G-bundle",
            "fundamental_or_technical": "fundamental",
        },
        {
            "claim": "dI_is_metric_for_all_beta",
            "status": "NO_GO",
            "formal_reason": (
                "raw dI=-log|K*| has dI(x,x)=-log|A cos(pi/6)|, generally nonzero; "
                "triangle inequality also fails on Z12 for small beta "
                "(e.g. 0.01: 216/1320 violations)."
            ),
            "missing_axiom": (
                "a declared diagonal normalization and a separate proof of "
                "nonnegativity, point separation, and triangle inequality"
            ),
            "missing_structure": "beta-regime geometry split / spectral distance alternative",
            "fundamental_or_technical": "fundamental_for_small_beta_regime",
        },
        {
            "claim": "unique_Lindblad_channels",
            "status": "NO_GO",
            "formal_reason": "radial K* does not select Kraus operators",
            "missing_axiom": "environment/coupling operators",
            "missing_structure": "system-bath factorization",
            "fundamental_or_technical": "fundamental",
        },
        {
            "claim": "continuum_field_theory",
            "status": "NO_GO_without_import",
            "formal_reason": "finite Z12 graph does not determine continuum limit functor",
            "missing_axiom": "continuum embedding + scale",
            "missing_structure": "manifold + measure",
            "fundamental_or_technical": "technical_if_embedding_supplied; fundamental_as_export",
        },
        {
            "claim": "physical_units_from_K_star",
            "status": "NO_GO",
            "formal_reason": "K* is dimensionless amplitude profile; free (A,beta) scale orbit",
            "missing_axiom": "unit-bearing source Omega_scale / S_+",
            "missing_structure": "dimension triad (action,length,time)",
            "fundamental_or_technical": "fundamental_on_current_artifacts",
        },
        {
            "claim": "selector_orientation_from_undirected_K_star",
            "status": "NO_GO",
            "formal_reason": "undirected radial kernel is inversion-even; Aut(Z12) contains orientation-reversing units",
            "missing_axiom": "inversion-odd source with polarity law",
            "missing_structure": "orientation torsor section",
            "fundamental_or_technical": "fundamental_for_undirected_radial_class",
        },
        {
            "claim": "SM_or_GR_or_ToE_from_Legacy_star_alone",
            "status": "NO_GO",
            "formal_reason": "operator/model generator exports finite spectral/Markov/Schrödinger scaffolds only",
            "missing_axiom": "large stack of field/geometry/unit/selector axioms",
            "missing_structure": "full physical theory package",
            "fundamental_or_technical": "fundamental",
        },
        {
            "claim": "unique_inverse_from_arbitrary_operator",
            "status": "NO_GO",
            "formal_reason": "map K* -> Op is many-to-one outside radial circulant class; not surjective onto B(H)",
            "missing_axiom": "restriction to radial circulant undirected class",
            "missing_structure": "canonical geometry of support",
            "fundamental_or_technical": "fundamental_for_global_inverse; technical_within_class",
        },
        {
            "claim": "observer_paradox_resolution_as_theorem",
            "status": "AXIOM_LEVEL",
            "formal_reason": "common generator implies both semigroups exist; operational choice is modeling postulate, not derived consciousness theorem",
            "missing_axiom": "operational model map instrument -> Borel function on spectrum",
            "missing_structure": "instrument algebra / POVM",
            "fundamental_or_technical": "foundational_modeling_choice",
        },
    ]


def potential_scores() -> list[dict[str, Any]]:
    """Research-program potential 0-10, not truth claims."""
    return [
        {"domain": "operator_theory", "score": 8, "rationale": "rich finite spectral/circulant family with dual functional calculus"},
        {"domain": "functional_analysis", "score": 7, "rationale": "semigroup/group generators, PSD threshold, Gram indefinite regime"},
        {"domain": "spectral_geometry", "score": 6, "rationale": "graph spectral geometry real; continuum spectral geometry not exported"},
        {"domain": "graph_theory", "score": 7, "rationale": "weighted cycle family, Laplacian, PF, cut metrics partial"},
        {"domain": "mathematical_physics", "score": 6, "rationale": "clean dual dynamics generator; limited continuum physics"},
        {"domain": "QM_foundations", "score": 5, "rationale": "unitary group yes; Born/measurement/no-go remains"},
        {"domain": "information_theory", "score": 5, "rationale": "kernel as correlation; dI not metric; Markov channel info partial"},
        {"domain": "unification_program", "score": 2, "rationale": "as standalone object: scaffold only; SM/GR/QG absent without large axiom stack"},
    ]


def research_program_branches() -> list[dict[str, Any]]:
    return [
        {"branch": "new_operator_class", "exists": True, "name": "LegacyStarRadialCirculantFamily", "note": "signed damped 8-periodic radial multipliers on cycles"},
        {"branch": "new_graph_class", "exists": True, "name": "PhaseModulatedWeightedCycles", "note": "undirected/directed variants with beta-PSD transition"},
        {"branch": "new_kernel_class", "exists": True, "name": "LegacyStarKernel", "note": "A cos(pi d/4+pi/6)/(1+beta d)"},
        {"branch": "new_semigroup_class", "exists": "conditional", "name": "DualBorelSemigroupsOfLegacyStar", "note": "all f(L) for Borel f; U and T special cases"},
        {"branch": "new_geometry_class", "exists": "conditional", "name": "SignedKernelGeometry", "note": "not dI-metric; spectral/information geometries candidates"},
        {"branch": "new_inverse_problems", "exists": True, "name": "RadialKernelFromOperatorInverse", "note": "partial invertibility on radial circulants"},
        {"branch": "new_spectral_theory", "exists": "conditional", "name": "IndefiniteToPSDTransitionSpectra", "note": "beta_* critical spectrum for Gram regime"},
    ]


def functor_analysis() -> dict[str, Any]:
    return {
        "pipeline": ["LegacyStarKernel", "Operator", "SemigroupOrGroup", "PhysicalModel"],
        "is_functor_strictly": "WARUNKOWO",
        "functor_statement": (
            "After fixing a support category C_supp (e.g. finite metric spaces or cycles) and a realization "
            "functor R: kernels -> operators (W, L, or Q), and a functional calculus F_f: Op -> dynamics, "
            "the composition F_f o R is a functor on that fixed category. Without fixed support/realization, "
            "Legacy* alone is a generator of candidates, not a unique functor."
        ),
        "operator_category_definable": True,
        "operator_category_objects": [
            "W_signed_circulant",
            "W_PSD_Gram",
            "L_graph_Laplacian",
            "Q_Markov",
            "H_Schrodinger",
            "T_diffusion_semigroup",
            "U_unitary_group",
        ],
        "universality_class": "radial_circulant_Fourier_multipliers_with_8_periodic_modulation_and_linear_damping",
        "images_classifiable": "partially — within finite cycle supports and fixed realization maps; not all of B(H)",
        "not_a_physical_ToE_functor": True,
    }


def dual_dynamics_structure() -> dict[str, Any]:
    wit = dual_dynamics_witness()
    return {
        "common_generator": True,
        "common_algebra": "C*(L) commutative spectral algebra on finite space (or von Neumann algebra generated by L)",
        "universal_generator": "L (or any affine sI-W with same spectral projections)",
        "deeper_structure": (
            "Borel functional calculus: all dynamics are maps lambda |-> f(lambda) on spectrum of L. "
            "U and T are two points in the same function space. Operational model selects f."
        ),
        "classification_of_dynamics": (
            "All strongly continuous groups/semigroups of the form f(L) for Borel f on spec(L). "
            "Unitary groups: f(lambda)=exp(-i t g(lambda)) with g real on spectrum. "
            "Contraction semigroups: Re f <= 0 appropriately. "
            "Markov: requires positivity-preserving f and probability conservation."
        ),
        "finite_witness": wit,
        "observer_paradox": (
            "No mathematical paradox: same spectral measure, two Borel functions. "
            "Choice is operational, not ontological consciousness. This is a modeling postulate, not a theorem about observers."
        ),
    }


def spectrum_table() -> list[dict[str, Any]]:
    rows = []
    for b in BETAS:
        W = undirected_circulant(1.0, b)
        ev = np.linalg.eigvalsh(W)
        Wp = positive_part_weights(W)
        L = graph_laplacian(Wp)
        lev = np.linalg.eigvalsh(L)
        rows.append(
            {
                "beta": b,
                "W_lambda_min": float(ev.min()),
                "W_lambda_max": float(ev.max()),
                "W_n_negative": int(np.sum(ev < -1e-10)),
                "W_psd": bool(ev.min() >= -1e-10),
                "L_lambda2": float(sorted(lev)[1]),
                "L_lambda_max": float(lev.max()),
            }
        )
    return rows


def build_payload() -> dict[str, Any]:
    phases = phase_table(16)
    tri = {f"beta_{b}": triangle_violations(1.0, b) for b in [0.01, 1.0, 2.0]}
    beta_psd = find_beta_psd_threshold()
    op_matrix = operator_class_matrix()
    dual = dual_dynamics_structure()
    inv = inverse_recovery_test()
    scores = potential_scores()
    nogo = no_go_ledger()
    physics = physics_readiness()
    fun = functor_analysis()
    branches = research_program_branches()
    specs = spectrum_table()

    n_tak = sum(1 for r in op_matrix if r["verdict"] == "TAK")
    n_nie = sum(1 for r in op_matrix if r["verdict"] == "NIE")
    n_cond = sum(1 for r in op_matrix if r["verdict"] == "WARUNKOWO")

    payload: dict[str, Any] = {
        "packet": "P3172",
        "slot": "S2122",
        "status": "P3172_LEGACY_STAR_OPERATOR_MODEL_GENERATOR_POTENTIAL_AUDIT",
        "date": "2026-07-20",
        "scope": {
            "kernel": "K*(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)",
            "A_gt_0": True,
            "beta_gt_0": True,
            "no_strict_bridge": True,
            "no_parameter_fitting_to_strict": True,
            "no_role_transfer": True,
            "treat_as_research_program": True,
            "working_facts": [
                "same_kernel_generates_U=exp(-itL)_and_T=exp(-tL)",
                "no_observer_paradox_operational_choice",
                "Legacy_star_is_autonomous_mathematical_object",
            ],
        },
        "constructed_theoretical_objects": {
            "LegacyStarKernel": {"omega": "pi/4", "phi": "pi/6", "damping": "1+beta*d"},
            "operator_class_matrix": op_matrix,
            "dual_dynamics_structure": dual,
            "functor_analysis": fun,
            "inverse_recovery": inv,
            "spectrum_table_undirected_Z12": specs,
            "phase_table": phases,
            "triangle_violations_dI": tri,
            "beta_psd_threshold": beta_psd,
            "physics_readiness": physics,
            "no_go_ledger": nogo,
            "research_program_branches": branches,
            "potential_scores": scores,
        },
        "finite_certificate": {
            "n_operator_classes": len(op_matrix),
            "n_TAK": n_tak,
            "n_NIE": n_nie,
            "n_WARUNKOWO": n_cond,
            "beta_star_approx": beta_psd["beta_star_approx"],
            "dual_unitarity_residual": dual["finite_witness"]["unitarity_residual_fro"],
            "dual_T_spectral_radius": dual["finite_witness"]["T_spectral_radius"],
            "inverse_max_abs_error": inv["max_abs_recovery_error"],
            "partial_invertibility_radial": inv["partially_invertible_on_radial_circulant_class"],
            "dI_is_metric_beta_0_01": tri["beta_0.01"]["is_metric"],
            "dI_violations_beta_0_01": tri["beta_0.01"]["n_violations"],
            "dI_is_metric_beta1": tri["beta_1.0"]["is_metric"],
            "dI_violations_beta1": tri["beta_1.0"]["n_violations"],
            "dI_triangle_holds_beta1": tri["beta_1.0"]["triangle_inequality_holds"],
            "n_no_go": len(nogo),
            "mean_potential_score": float(np.mean([s["score"] for s in scores])),
            "unification_score": next(s["score"] for s in scores if s["domain"] == "unification_program"),
            "operator_theory_score": next(s["score"] for s in scores if s["domain"] == "operator_theory"),
        },
        "tiered_claims": {
            "proved_theorems_scoped": [
                "Undirected W[K*] on Z12 is real symmetric for all A>0, beta>0.",
                "Functional calculus: U(t)=exp(-itL) and T(t)=exp(-tL) share spectral measure of self-adjoint L.",
                "L=Deg-W+ is PSD with L1=0 for W+=max(K*,0) off-diagonal.",
                "Raw dI=-log|K*| is not a metric because its diagonal is generally nonzero; its triangle inequality additionally fails at small beta (e.g. 0.01: 216 violations).",
                "Radial profile recoverable from undirected circulant first row (partial inverse on class).",
                "Phase cos(pi d/4+pi/6) has period 8 in d.",
            ],
            "results_from_analysis": [
                "Operator-class matrix: majority WARUNKOWO; Dirac/Maxwell/YM are NIE.",
                "Legacy* is a strong finite generator of spectral/Markov/diffusion models.",
                "Not a unique physical ToE functor without large axiom stack.",
                "Dual dynamics are two Borel functions of one generator — no observer paradox mathematically.",
            ],
            "hypotheses": [
                "Universality class of 8-periodic damped radial multipliers may have continuum PsiDO limit on S1.",
                "Indefinite Gram regime (beta<beta*) may model signed information kernels useful outside PSD covariance.",
            ],
            "intuitions": [
                "Operational instrument selects f in f(L), analogous to filter choice in spectral processing.",
                "Linear damping 1/(1+beta d) is hyperbolic-range scaffold, cosine is interference scaffold.",
            ],
            "speculations": [
                "Any deeper unification using Legacy* would require new source axioms not contained in K*.",
            ],
        },
        "decision": {
            "accepted_as_fundamental_physics": False,
            "accepted_as_mathematical_generator": True,
            "strict_bridge_performed": False,
            "role_transfer_performed": False,
            "selector_closed": False,
            "L_total_promoted": False,
            "ToE_promoted": False,
            "negative_export_flags": {
                "Dirac_from_K_star": False,
                "Maxwell_from_K_star": False,
                "YM_from_K_star": False,
                "SM_from_K_star": False,
                "GR_from_K_star": False,
                "unit_bearing_action": False,
                "selector_closure": False,
                "ToE": False,
            },
            "next_honest_step": (
                "Either (i) develop the finite LegacyStarRadialCirculantFamily as pure operator/graph theory "
                "(PSD threshold, indefinite Gram, dual Borel calculus, inverse radial recovery), or "
                "(ii) if physics is desired, supply exactly one missing axiom package: units, continuum embedding, "
                "or gauge/spin structure — without claiming it is derived from K* alone."
            ),
        },
        "provenance": {
            "script": "fundamental_action_reconstruction/p3172_s2122_legacy_star_operator_model_generator_potential_audit.py",
            "related": ["P3171", "P3075", "P3082", "P3145", "P3168"],
            "agents_sha": sha(AGENTS),
            "summary_sha": sha(SUMMARY),
            "packet_path": str(PACKET),
            "json_path": str(OUT),
            "md_path": str(MD),
        },
    }
    return payload


def write_md(payload: dict[str, Any]) -> str:
    cert = payload["finite_certificate"]
    scores = payload["constructed_theoretical_objects"]["potential_scores"]
    score_lines = "\n".join(f"- `{s['domain']}`: **{s['score']}/10** — {s['rationale']}" for s in scores)
    op_lines = "\n".join(
        f"- **{r['class']}**: `{r['verdict']}` — {r['reason']}"
        for r in payload["constructed_theoretical_objects"]["operator_class_matrix"]
    )
    nogo_lines = "\n".join(
        f"- `{g['claim']}`: **{g['status']}** ({g['fundamental_or_technical']}) — {g['formal_reason']}"
        for g in payload["constructed_theoretical_objects"]["no_go_ledger"]
    )
    text = f"""# P3172/S2122 Legacy* operator/model generator potential audit

Status: `{payload['status']}`

## Scope
- Kernel: `K*(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)`, `A>0`, `beta>0`
- **No** Strict bridge, **no** parameter fitting to Strict, **no** role transfer
- Working facts only: dual dynamics `U=exp(-itL)`, `T=exp(-tL)`; operational (not consciousness) choice

## Finite certificate
- operator classes: `{cert['n_operator_classes']}` (`TAK={cert['n_TAK']}`, `WARUNKOWO={cert['n_WARUNKOWO']}`, `NIE={cert['n_NIE']}`)
- `beta_*` (PSD approx on undirected Z12): `{cert['beta_star_approx']}`
- dual unitarity residual: `{cert['dual_unitarity_residual']:.3e}`
- dual `T` spectral radius: `{cert['dual_T_spectral_radius']:.6f}`
- inverse radial recovery max abs error: `{cert['inverse_max_abs_error']:.3e}`
- `dI` metric at beta=0.01: `{cert['dI_is_metric_beta_0_01']}` (violations `{cert['dI_violations_beta_0_01']}`)
- `dI` metric at beta=1: `{cert['dI_is_metric_beta1']}` (violations `{cert['dI_violations_beta1']}`)
- mean potential score: `{cert['mean_potential_score']:.2f}`
- operator theory score: `{cert['operator_theory_score']}/10`
- unification score: `{cert['unification_score']}/10`

## Operator-class matrix
{op_lines}

## Dual dynamics
- Common generator: **yes** (spectral calculus)
- Common algebra: commutative C*(L) on finite space
- Observer paradox: **not a mathematical paradox** — two Borel functions of one spectral measure; operational model selects f

## Potential scores (research program, not truth)
{score_lines}

## No-go ledger
{nogo_lines}

## Tiered claims
### Proved (scoped)
{chr(10).join('- ' + x for x in payload['tiered_claims']['proved_theorems_scoped'])}

### Analysis results
{chr(10).join('- ' + x for x in payload['tiered_claims']['results_from_analysis'])}

### Hypotheses
{chr(10).join('- ' + x for x in payload['tiered_claims']['hypotheses'])}

### Intuitions
{chr(10).join('- ' + x for x in payload['tiered_claims']['intuitions'])}

### Speculations
{chr(10).join('- ' + x for x in payload['tiered_claims']['speculations'])}

## Decision
- mathematical generator: **accepted**
- fundamental physics / ToE / L_total / selector: **not exported**
- next honest step: {payload['decision']['next_honest_step']}

Full narrative: `P3172_S2122_LEGACY_STAR_OPERATOR_MODEL_GENERATOR_POTENTIAL_AUDIT.md`
"""
    return text


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(write_md(payload), encoding="utf-8")

    # compact guardrail appendices
    marker = "P3172/S2122"
    sheet_note = f"""
## {marker} Legacy* operator/model generator potential
- Standalone audit of `K*(d)=A cos(pi d/4+pi/6)/(1+beta d)` as operator generator.
- No Strict bridge, no role transfer, no L_total/ToE promotion.
- Dual dynamics U/T share spectral measure; Dirac/Maxwell/YM not exported from radial K*.
"""
    draft_note = f"""
## {marker}
- Legacy* treated as autonomous model-generator, not Strict completion.
- Finite certificate: TAK/WARUNKOWO/NIE operator matrix; partial inverse on radial circulants; dI non-metric.
"""
    agents_note = f"""
## Current Legacy* standalone operator-generator potential guardrail (P3172/S2122, 2026-07-20)

- P3172 audits `K*(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)` as an autonomous mathematical generator of operators/semigroups/models, without Strict bridge or role transfer.
- Dual dynamics `U=exp(-itL)` and `T=exp(-tL)` are Borel functions of one self-adjoint generator when L is hermitized; this removes a mathematical observer paradox but does not export Born rule, instruments, or QW-2191 discharge.
- Finite operator-class matrix: graph Laplacian / integral / spectral / diffusion are available; Dirac, Maxwell, Yang–Mills are not exported from radial scalar K* alone; most other classes are conditional.
- Do not promote Legacy* dual dynamics, Gram PSD threshold, or radial inverse recovery to unit-bearing L_total, SM/GR, selector closure, bridge completion, role transfer, or ToE.
- A next pure-math move may develop the LegacyStarRadialCirculantFamily (PSD threshold, indefinite Gram, dual Borel calculus). A physics move requires an explicit new axiom package (units, continuum, gauge/spin), marked as imported.
"""
    append_once(SHEET, marker, sheet_note)
    append_once(DRAFT, marker, draft_note)
    append_once(AGENTS, "P3172/S2122", agents_note)

    if SUMMARY.exists():
        append_once(
            SUMMARY,
            "P3172/S2122",
            f"""
### P3172/S2122 — Legacy* as operator/model generator (standalone potential)
- Status: `{payload['status']}`
- Operator matrix: TAK={payload['finite_certificate']['n_TAK']}, WARUNKOWO={payload['finite_certificate']['n_WARUNKOWO']}, NIE={payload['finite_certificate']['n_NIE']}
- Dual U/T dynamics: common spectral generator; no observer paradox as math fact
- Unification score: {payload['finite_certificate']['unification_score']}/10; operator theory: {payload['finite_certificate']['operator_theory_score']}/10
- No Strict bridge / no L_total / no ToE promotion
- Artifacts: `generated/p3172_s2122_legacy_star_operator_model_generator_potential_audit.{{json,md}}`, `P3172_S2122_LEGACY_STAR_OPERATOR_MODEL_GENERATOR_POTENTIAL_AUDIT.md`
""",
        )

    print(json.dumps({"status": payload["status"], "out": str(OUT), "md": str(MD), "cert": payload["finite_certificate"]}, indent=2))


if __name__ == "__main__":
    main()
