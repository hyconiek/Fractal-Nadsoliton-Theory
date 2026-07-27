#!/usr/bin/env python3
"""FIN Programs 178--190: composition, process memory, and scale obstructions.

The executable is deliberately adversarial.  It distinguishes analytic
theorems, exact finite checks, synthetic evidence, conditional constructions,
and unavailable proof/certification toolchains.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import hashlib
import importlib.util
import json
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "FIN_Programs_178_190_Composition_Process_Scale_Figures"
RESULTS = ROOT / "FIN_Programs_178_190_Composition_Process_Scale_Results.json"
OPEN_SET_PROTOCOL = ROOT / "FIN_Programs_178_190_Open_Set_Preregistration.json"
PREVIOUS = ROOT / "FIN_Programs_165_177_Axiom_Falsification_Measurement_Results.json"
LEAN_SOURCE = ROOT / "FIN_Programs_178_190_Dirichlet_Core.lean"

SEED = 20260728
RNG = np.random.default_rng(SEED)
ALPHA_GEO = 4 * math.log(2)
ETA_STRICT = 1.8
OMEGA_STRICT = 0.18575
PHI_STRICT = 0.1625
OMEGA_LEGACY = math.pi / 4
PHI_LEGACY = math.pi / 6


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canonical_digest(obj: object) -> str:
    payload = json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def stable_rvs(alpha: float, size: tuple[int, ...], rng: np.random.Generator) -> np.ndarray:
    """Symmetric alpha-stable variates (Chambers--Mallows--Stuck)."""
    u = rng.uniform(-math.pi / 2, math.pi / 2, size)
    w = rng.exponential(1.0, size)
    if abs(alpha - 1.0) < 1e-12:
        return np.tan(u)
    numerator = np.sin(alpha * u)
    denominator = np.cos(u) ** (1 / alpha)
    correction = (
        np.cos((1 - alpha) * u) / w
    ) ** ((1 - alpha) / alpha)
    return numerator / denominator * correction


def program178_tensor_intensive_state_laws() -> dict:
    alpha = ALPHA_GEO
    r = math.log(4.0)
    copies = np.arange(1, 9)

    candidates = {
        "raw_cardinality_alpha_over_h": lambda a, rr: a / math.exp(rr),
        "replication_invariant_alpha_over_log_h": lambda a, rr: a / rr,
        "inverse_ratio_log_h_over_alpha": lambda a, rr: rr / a,
        "constant_log2": lambda a, rr: math.log(2),
        "extensive_alpha": lambda a, rr: a,
    }
    rows = []
    series = {}
    for name, fun in candidates.items():
        base = fun(alpha, r)
        values = np.asarray([fun(n * alpha, n * r) for n in copies], dtype=float)
        # Coarse graining: split each sector into m labels without changing alpha.
        coarse = np.asarray([fun(alpha, r + math.log(m)) for m in (1, 2, 3, 5)])
        row = {
            "candidate": name,
            "one_copy": base,
            "max_replication_defect": float(np.max(np.abs(values - base))),
            "max_label_split_defect": float(np.max(np.abs(coarse - coarse[0]))),
            "replication_intensive": bool(np.allclose(values, base, atol=1e-12)),
            "label_split_invariant": bool(np.allclose(coarse, coarse[0], atol=1e-12)),
        }
        rows.append(row)
        series[name] = values.tolist()

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    for name, values in series.items():
        ax.plot(copies, values, marker="o", label=name.replace("_", " "))
    ax.set_yscale("log")
    ax.set_xlabel("independent copy number n")
    ax.set_ylabel("candidate inverse temperature")
    ax.set_title("Program 178: replication intensivity does not select a state law")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    fig.savefig(FIG / "program178_tensor_intensive_classification.png", dpi=190)
    plt.close(fig)

    return {
        "replication_action": "(alpha,r=log h) -> (n alpha,n r)",
        "classification_theorem": (
            "Every continuous degree-zero scalar under positive replication "
            "has the form F(alpha/r) on the positive cone. Replication "
            "intensivity therefore leaves an arbitrary function F."
        ),
        "coarse_graining_no_go": (
            "If arbitrary label splitting changes r while alpha is fixed and "
            "beta is invariant, F(alpha/r) must be constant. The constant is "
            "not selected by composition and cannot derive ln(2)."
        ),
        "candidate_rows": rows,
        "unique_target_free_beta_selected": False,
        "constructed_object": (
            "StateLawModuli = continuous functions F on alpha/log(h); "
            "coarse-graining quotient collapses this moduli space to "
            "unselected constants."
        ),
        "status": "Classification/no-go theorem in the declared composition class",
        "confidence": "Proven",
        "claim_boundary": (
            "A richer state law may use new typed data, but composition and "
            "label invariance alone do not select beta=ln(2)."
        ),
    }


def program179_positive_compression_source() -> dict:
    betas = np.logspace(-2, 2, 9)
    x = np.linspace(0, 8, 500)
    orbit_residuals = []
    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    for beta in betas:
        d = x / beta ** (1 / ETA_STRICT)
        denominator_profile = 1 / (1 + beta * d**ETA_STRICT)
        universal = 1 / (1 + x**ETA_STRICT)
        residual = float(np.max(np.abs(denominator_profile - universal)))
        orbit_residuals.append({"beta": float(beta), "collapse_residual": residual})
        ax.plot(x, denominator_profile, alpha=0.55)
    ax.plot(x, 1 / (1 + x**ETA_STRICT), color="black", lw=2, ls="--", label="universal orbit profile")
    ax.set_xlabel(r"rescaled coordinate $x=\beta^{1/\eta}d$")
    ax.set_ylabel("compression factor")
    ax.set_title("Program 179: every positive beta lies on one scale orbit")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program179_compression_scale_torsor.png", dpi=190)
    plt.close(fig)

    candidates = [
        {
            "source": "positivity_of_denominator",
            "positive": True,
            "target_independent": True,
            "selects_value_one": False,
            "verdict": "selects only beta>=0",
        },
        {
            "source": "monoid_unitality_A(1)=1",
            "positive": True,
            "target_independent": True,
            "selects_value_one": False,
            "verdict": "fixes intercept, not compression slope",
        },
        {
            "source": "tail_ratio_inversion",
            "positive": True,
            "target_independent": False,
            "selects_value_one": True,
            "verdict": "receiver inversion uses strict target",
        },
        {
            "source": "A_ME_cardinality",
            "positive": True,
            "target_independent": False,
            "selects_value_one": False,
            "verdict": "fails tensor intensivity and needs normalization",
        },
        {
            "source": "coordinate_choice_beta=1",
            "positive": True,
            "target_independent": True,
            "selects_value_one": True,
            "verdict": "gauge convention, not source law",
        },
        {
            "source": "external_length_calibration",
            "positive": True,
            "target_independent": True,
            "selects_value_one": True,
            "verdict": "conditional A4 input, not strict",
        },
    ]
    accepted = [
        row
        for row in candidates
        if row["positive"] and row["target_independent"] and row["selects_value_one"]
        and row["source"] not in {"coordinate_choice_beta=1", "external_length_calibration"}
    ]
    return {
        "constructed_object": {
            "name": "CompressionScaleTorsor",
            "carrier": "R_{>0}",
            "action": "c · beta = beta*c^(-eta) under d' = c*d",
            "orbit_statement": "the action is transitive on beta>0",
        },
        "orbit_collapse": orbit_residuals,
        "source_candidates": candidates,
        "accepted_strict_sources": len(accepted),
        "source_theorem": (
            "Without a selected length coordinate, all beta>0 are equivalent "
            "under rescaling. Invariant data can select the positive orbit but "
            "cannot select beta=1."
        ),
        "strict_beta_source_exported": False,
        "status": "Positive-orbit theorem; target-independent value source absent",
        "confidence": "Proven",
        "claim_boundary": (
            "The theorem does not reject beta=1 as a working gauge; it rejects "
            "its promotion to a strict physical constant without a scale frame."
        ),
    }


def program180_ball_certificate_assembly(previous: dict) -> dict:
    p165 = previous["programs"]["165"]
    p151_path = ROOT / "FIN_Programs_151_164_Axiomatic_Operational_Results.json"
    p151 = json.loads(p151_path.read_text(encoding="utf-8"))["programs"]["151"]
    components = [
        {
            "component": "finite FFT cells",
            "available": True,
            "directed_rounding": True,
            "closed": True,
            "value": p151["maximum_relative_remainder_upper"],
        },
        {
            "component": "average/polylogarithmic term",
            "available": True,
            "directed_rounding": False,
            "closed": False,
            "value": None,
        },
        {
            "component": "cancellation correction tail",
            "available": True,
            "directed_rounding": False,
            "closed": False,
            "value": p165["best_cancellation_bound"],
        },
        {
            "component": "denominator replacement",
            "available": True,
            "directed_rounding": False,
            "closed": False,
            "value": p165["rows"][-1]["denominator_replacement_bound"],
        },
        {
            "component": "resonant interval fallbacks",
            "available": True,
            "directed_rounding": False,
            "closed": False,
            "value": p165["rows"][-1]["resonant_mode_intervals"],
        },
    ]
    arb = bool(importlib.util.find_spec("flint")) or bool(importlib.util.find_spec("sageall"))
    colors = ["#19733A" if row["closed"] else "#D55E00" for row in components]
    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    ax.barh([r["component"] for r in components], [1] * len(components), color=colors)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("Program 180: one-engine certificate assembly status")
    for i, row in enumerate(components):
        ax.text(0.03, i, "closed" if row["closed"] else "open in directed rounding", va="center", color="white", weight="bold")
    fig.savefig(FIG / "program180_ball_certificate_ledger.png", dpi=190)
    plt.close(fig)

    return {
        "arb_or_python_flint_available": arb,
        "components": components,
        "formal_worst_relative_enclosure_from_P151": p151["maximum_relative_remainder_upper"],
        "analytic_cancellation_bound_from_P165": p165["best_cancellation_bound"],
        "preregistered_target": 0.03,
        "full_one_engine_certificate": False,
        "bottleneck_migration": (
            "Cancellation is no longer the dominant analytic obstruction. "
            "The open obligation is a common directed-rounding evaluation of "
            "the average term, correction modes, resonant fallbacks, and cells."
        ),
        "status": "Certificate assembly executed; full sub-3% certificate remains open",
        "confidence": "Proven ledger; no false ball-arithmetic claim",
    }


def program181_finite_dirichlet_core() -> dict:
    # Four-cycle with W=1/2 on each of two neighbours, hence row sum one.
    W = [[Fraction(0) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        W[i][(i - 1) % 4] = Fraction(1, 2)
        W[i][(i + 1) % 4] = Fraction(1, 2)
    A = [
        [(Fraction(1) if i == j else Fraction(0)) - W[i][j] for j in range(4)]
        for i in range(4)
    ]
    f = [Fraction(1), Fraction(2), Fraction(-1), Fraction(3)]
    lhs = sum(f[i] * A[i][j] * f[j] for i in range(4) for j in range(4))
    rhs = Fraction(1, 2) * sum(
        W[i][j] * (f[i] - f[j]) ** 2 for i in range(4) for j in range(4)
    )
    constant_residual = [
        sum(A[i][j] for j in range(4)) for i in range(4)
    ]

    Af = np.asarray([[float(x) for x in row] for row in A])
    times = np.linspace(0, 3, 31)
    u_res = []
    p_res = []
    p_min = []
    for t in times:
        U = expm(-1j * t * Af)
        P = expm(-t * Af)
        u_res.append(float(np.linalg.norm(U.conj().T @ U - np.eye(4))))
        p_res.append(float(np.max(np.abs(P.sum(axis=1) - 1))))
        p_min.append(float(P.min()))

    eig = np.linalg.eigvalsh(Af)
    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.plot(times, u_res, label="unitarity residual")
    ax.plot(times, p_res, label="row-sum residual")
    ax.plot(times, np.maximum(0, -np.asarray(p_min)), label="negative-entry defect")
    ax.set_yscale("symlog", linthresh=1e-17)
    ax.set_xlabel("t")
    ax.set_ylabel("numerical residual")
    ax.set_title("Program 181: finite Dirichlet dynamics")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program181_dirichlet_core.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": "FiniteDirichletSystem(H,W,A=I-W)",
        "exact_test_vector": [int(x) for x in f],
        "exact_quadratic_form": str(lhs),
        "exact_dirichlet_sum": str(rhs),
        "exact_identity": lhs == rhs,
        "exact_constant_kernel": all(x == 0 for x in constant_residual),
        "eigenvalues": eig.tolist(),
        "max_unitarity_residual": max(u_res),
        "max_stochastic_row_residual": max(p_res),
        "minimum_heat_entry": min(p_min),
        "lean_source": LEAN_SOURCE.name,
        "lean_source_sha256": sha256(LEAN_SOURCE),
        "lean_machine_compiled": False,
        "lean_note": (
            "The core-only source contains concrete finite theorems. No local "
            "Lean toolchain was available after the documented Release-10.16 "
            "download failures, so no machine-compilation claim is made."
        ),
        "status": "Exact algebraic core executed; proof-assistant source uncompiled",
        "confidence": "Proven exact finite identity; numerical exponential checks",
        "claim_boundary": (
            "The computation proves a finite model, not a continuum theorem "
            "or a physical interpretation."
        ),
    }


def max_covariant_transverse_factor(z_source: float, z_target: float, grid: int = 4001) -> tuple[float, float, float]:
    p = (1 + z_source) / 2
    q = (1 + z_target) / 2
    if p <= 1e-12 or p >= 1 - 1e-12:
        return 0.0, math.nan, math.nan
    a = np.linspace(0, 1, grid)
    c = (q - p * a) / (1 - p)
    valid = (c >= 0) & (c <= 1)
    factors = np.full_like(a, -np.inf)
    factors[valid] = np.sqrt(a[valid] * (1 - c[valid])) + np.sqrt(
        (1 - a[valid]) * c[valid]
    )
    idx = int(np.argmax(factors))
    return float(factors[idx]), float(a[idx]), float(c[idx])


def program182_complete_reflection_convertibility() -> dict:
    xs, zs = 0.6, 0.0
    xt, zt = 0.6, 0.8
    factor, aopt, copt = max_covariant_transverse_factor(zs, zt)
    max_xt = xs * factor

    zgrid = np.linspace(-0.98, 0.98, 121)
    xgrid = np.linspace(0, 1, 121)
    feasible = np.zeros((len(zgrid), len(xgrid)))
    valid = np.zeros_like(feasible)
    m_false_positive = 0
    valid_count = 0
    for iz, z in enumerate(zgrid):
        fac, _, _ = max_covariant_transverse_factor(zs, z, grid=801)
        for ix, x in enumerate(xgrid):
            if x * x + z * z <= 1 + 1e-12:
                valid[iz, ix] = 1
                valid_count += 1
                ok = x <= xs * fac + 1e-8
                feasible[iz, ix] = 1 if ok else 0
                if x <= xs + 1e-12 and not ok:
                    m_false_positive += 1

    fig, ax = plt.subplots(figsize=(8.8, 5.5), constrained_layout=True)
    image = ax.imshow(
        np.ma.masked_where(valid == 0, feasible),
        origin="lower",
        extent=[xgrid[0], xgrid[-1], zgrid[0], zgrid[-1]],
        aspect="auto",
        cmap="RdYlGn",
        vmin=0,
        vmax=1,
    )
    ax.scatter([xt], [zt], color="black", marker="x", s=70, label="equal-M counterexample")
    ax.set_xlabel("target transverse magnitude")
    ax.set_ylabel("target z")
    ax.set_title("Program 182: exact Choi-feasibility boundary from source (0.6,0)")
    ax.legend()
    fig.colorbar(image, ax=ax, ticks=[0, 1], label="covariant conversion feasible")
    fig.savefig(FIG / "program182_reflection_convertibility.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "Z2CovariantQubitConversionLMI",
            "choi_form": "[[a,0,0,w],[0,1-a,zeta,0],[0,zeta*,c,0],[w*,0,0,1-c]]",
            "positivity": "|w|^2<=a(1-c), |zeta|^2<=(1-a)c, 0<=a,c<=1",
        },
        "complete_real_plane_criterion": (
            "For source (r,z) and target (r',z'), conversion exists iff "
            "there are a,c in [0,1] with p'=p*a+(1-p)*c and "
            "r'<=r[sqrt(a(1-c))+sqrt((1-a)c)]. Free z rotations extend "
            "the criterion to arbitrary transverse phase."
        ),
        "counterexample": {
            "source": [xs, zs],
            "target": [xt, zt],
            "same_M": True,
            "maximum_transverse_factor": factor,
            "maximum_reachable_target_M": max_xt,
            "optimal_a": aopt,
            "optimal_c": copt,
            "conversion_feasible": xt <= max_xt + 1e-9,
        },
        "grid_valid_targets": valid_count,
        "grid_false_positives_of_M_only": m_false_positive,
        "M_complete_on_full_state_space": False,
        "status": "Complete Choi/LMI criterion constructed for qubit reflection covariance",
        "confidence": "Proven criterion; finite boundary evaluated numerically",
        "claim_boundary": (
            "This classifies use of an existing asymmetry resource and does "
            "not generate the signed preparation or discharge QW-2191."
        ),
    }


def iqr_interval(sample: np.ndarray, eps: float) -> tuple[float, float]:
    lo = np.quantile(sample, max(0.0, 0.75 - eps)) - np.quantile(
        sample, min(1.0, 0.25 + eps)
    )
    hi = np.quantile(sample, min(1.0, 0.75 + eps)) - np.quantile(
        sample, max(0.0, 0.25 - eps)
    )
    return max(float(lo), np.finfo(float).tiny), max(float(hi), np.finfo(float).tiny)


def exponent_interval(x1: np.ndarray, x2: np.ndarray, eps: float, ratio: float) -> tuple[float, float]:
    l1, u1 = iqr_interval(x1, eps)
    l2, u2 = iqr_interval(x2, eps)
    return math.log(l2 / u1) / math.log(ratio), math.log(u2 / l1) / math.log(ratio)


def program183_block_robust_quantiles() -> dict:
    alpha = 0.8
    true_T = 1 / alpha
    blocks = 40
    block_size = 100
    nominal_n = blocks * block_size
    delta = 0.05
    ratio = 4.0
    eps_nominal = math.sqrt(math.log(4 / delta) / (2 * nominal_n))
    eps_block = math.sqrt(math.log(4 / delta) / (2 * blocks))
    reps = 360
    cover_nominal = 0
    cover_block = 0
    widths_nominal = []
    widths_block = []
    for _ in range(reps):
        x1 = stable_rvs(alpha, (blocks,), RNG)
        x2 = stable_rvs(alpha, (blocks,), RNG) * ratio**true_T
        lo, hi = exponent_interval(x1, x2, eps_nominal, ratio)
        blo, bhi = exponent_interval(x1, x2, eps_block, ratio)
        cover_nominal += lo <= true_T <= hi
        cover_block += blo <= true_T <= bhi
        widths_nominal.append(hi - lo)
        widths_block.append(bhi - blo)

    # An exact contiguous-run diagnostic for the declared repeated-block class.
    diagnostic_record = np.repeat(np.arange(blocks, dtype=float), block_size)
    run_lengths = []
    start = 0
    for i in range(1, len(diagnostic_record) + 1):
        if i == len(diagnostic_record) or diagnostic_record[i] != diagnostic_record[start]:
            run_lengths.append(i - start)
            start = i
    detected_block_size = int(np.median(run_lengths))

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.hist(widths_nominal, bins=35, alpha=0.65, label=f"invalid iid n={nominal_n}")
    ax.hist(widths_block, bins=35, alpha=0.65, label=f"valid block m={blocks}")
    ax.set_xlabel("two-time exponent interval width")
    ax.set_ylabel("count")
    ax.set_title("Program 183: dependence robustness costs effective sample size")
    ax.legend()
    fig.savefig(FIG / "program183_block_DKW.png", dpi=190)
    plt.close(fig)

    return {
        "dependence_class": "m iid latent draws, each repeated in one contiguous block of size b",
        "theorem": (
            "The repeated record has exactly the empirical CDF of the m latent "
            "draws. DKW therefore holds with m=n/b, not nominal n."
        ),
        "nominal_n": nominal_n,
        "independent_blocks": blocks,
        "block_size": block_size,
        "epsilon_invalid_nominal": eps_nominal,
        "epsilon_valid_block": eps_block,
        "replicates": reps,
        "coverage_invalid_nominal": cover_nominal / reps,
        "coverage_valid_block": cover_block / reps,
        "mean_width_invalid_nominal": float(np.mean(widths_nominal)),
        "mean_width_valid_block": float(np.mean(widths_block)),
        "run_length_detected_block_size": detected_block_size,
        "status": "Exact block-class theorem and finite coverage challenge executed",
        "confidence": "Proven for repeated iid blocks; strong numerical coverage evidence",
        "claim_boundary": (
            "The theorem does not cover arbitrary mixing or hidden dependence. "
            "Block size must be known or certified by the acquisition protocol."
        ),
    }


def program184_nonlinear_multi_control() -> dict:
    xvals = np.asarray([-1.0, 0.0, 1.0])
    # Columns aT,aC1,aC2,T,g1,g2 after known control exponents are subtracted.
    rows = []
    labels = []
    for group in range(3):
        for x in xvals:
            row = np.zeros(6)
            row[group] = 1
            if group == 0:
                row[3] = x
            row[4] = x
            row[5] = x * x
            rows.append(row)
            labels.append((group, x))
    Xcell = np.asarray(rows)
    rank = int(np.linalg.matrix_rank(Xcell))
    one_control_rank = int(np.linalg.matrix_rank(Xcell[:6, [0, 1, 3, 4, 5]]))

    reps = 1500
    per_cell = 12
    sigma = 0.12
    true_T, g1, g2 = 1.25, 0.7, 0.35
    known = {0: 0.0, 1: 0.5, 2: 1.0}
    intercept = {0: 0.2, 1: -0.1, 2: 0.4}
    estimates = []
    naive = []
    violation_detected = 0
    violation_delta = 0.25
    for _ in range(reps):
        X = np.repeat(Xcell, per_cell, axis=0)
        y = []
        for group, x in labels:
            exponent = true_T if group == 0 else known[group]
            mu = intercept[group] + exponent * x + g1 * x + g2 * x * x
            y.extend(mu + RNG.normal(0, sigma, per_cell))
        y = np.asarray(y)
        adjusted = y.copy()
        for cell, (group, x) in enumerate(labels):
            if group > 0:
                sl = slice(cell * per_cell, (cell + 1) * per_cell)
                adjusted[sl] -= known[group] * x
        fit = np.linalg.lstsq(X, adjusted, rcond=None)[0]
        estimates.append(fit[3])
        cell_means = y.reshape(9, per_cell).mean(axis=1)
        naive.append((cell_means[2] - cell_means[0]) / 2)

        # Separate challenge: control 2 has an extra linear gain delta.
        c1_left = known[1] * -1 + g1 * -1 + g2
        c1_right = known[1] + g1 + g2
        c2_left = known[2] * -1 + (g1 + violation_delta) * -1 + g2
        c2_right = known[2] + (g1 + violation_delta) + g2
        slope1 = (
            (c1_right + RNG.normal(0, sigma / math.sqrt(per_cell)))
            - (c1_left + RNG.normal(0, sigma / math.sqrt(per_cell)))
        ) / 2 - known[1]
        slope2 = (
            (c2_right + RNG.normal(0, sigma / math.sqrt(per_cell)))
            - (c2_left + RNG.normal(0, sigma / math.sqrt(per_cell)))
        ) / 2 - known[2]
        se_diff = sigma / math.sqrt(per_cell)
        violation_detected += abs(slope2 - slope1) / se_diff > 1.96

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.hist(naive, bins=40, alpha=0.6, label="target only")
    ax.hist(estimates, bins=40, alpha=0.65, label="two controls + quadratic gain")
    ax.axvline(true_T, color="black", ls="--", label="true target")
    ax.set_xlabel("estimated target exponent")
    ax.set_ylabel("count")
    ax.set_title("Program 184: multi-control nonlinear calibration")
    ax.legend()
    fig.savefig(FIG / "program184_multi_control.png", dpi=190)
    plt.close(fig)

    return {
        "model": "log Q_j=a_j+T_j*x+g1*x+g2*x^2, x=log t",
        "controls": {"T_C1": known[1], "T_C2": known[2]},
        "design_rank": rank,
        "parameters": 6,
        "one_control_reduced_rank": one_control_rank,
        "times_in_log_coordinates": xvals.tolist(),
        "simulation": {
            "replicates": reps,
            "records_per_cell": per_cell,
            "true_target": true_T,
            "naive_mean": float(np.mean(naive)),
            "controlled_mean": float(np.mean(estimates)),
            "controlled_sd": float(np.std(estimates, ddof=1)),
            "shared_gain_violation_delta": violation_delta,
            "shared_gain_violation_detection_rate": violation_detected / reps,
        },
        "constructed_object": (
            "MultiControlCalibrationFrame = target + two known-exponent "
            "controls + shared quadratic gain + a cross-control falsification contrast."
        ),
        "status": "Full-rank nonlinear calibration and assumption diagnostic executed",
        "confidence": "Proven rank; strong synthetic evidence",
        "claim_boundary": (
            "Identifiability remains conditional on the declared shared-gain "
            "family; the second control tests but cannot prove universality."
        ),
    }


def quantile_features(x1: np.ndarray, x2: np.ndarray, c1: np.ndarray, c2: np.ndarray, t2: float = 4.0) -> np.ndarray:
    def iqr(x: np.ndarray) -> np.ndarray:
        return np.quantile(x, 0.75, axis=1) - np.quantile(x, 0.25, axis=1)

    target_slope = np.log(iqr(x2) / iqr(x1)) / math.log(t2)
    control_slope = np.log(iqr(c2) / iqr(c1)) / math.log(t2)
    target_adjusted = target_slope - (control_slope - 0.5)
    q05, q25, q75, q95 = np.quantile(x1, [0.05, 0.25, 0.75, 0.95], axis=1)
    shape1 = (q95 - q05) / (q75 - q25)
    q05b, q25b, q75b, q95b = np.quantile(x2, [0.05, 0.25, 0.75, 0.95], axis=1)
    shape2 = (q95b - q05b) / (q75b - q25b)
    boundary_atoms = np.mean(np.isclose(np.abs(x2), 8.0), axis=1)
    return np.column_stack([target_slope, control_slope, target_adjusted, shape1, shape2, boundary_atoms])


def classify_features(features: np.ndarray, model: dict) -> tuple[np.ndarray, np.ndarray]:
    names = model["models"]
    scale = np.asarray(model["feature_scale"])
    centroids = {k: np.asarray(v) for k, v in model["centroids"].items()}
    thresholds = model["abstention_thresholds"]
    labels = []
    abstain = []
    for feature in features:
        dist = np.asarray([np.linalg.norm((feature - centroids[m]) / scale) for m in names])
        idx = int(np.argmin(dist))
        labels.append(idx)
        abstain.append(bool(dist[idx] > thresholds[names[idx]]))
    return np.asarray(labels), np.asarray(abstain)


def program185_open_set_challenge(previous: dict) -> dict:
    p174 = previous["programs"]["174"]
    core = {
        "protocol_id": "FIN-P185-OPEN-SET-CHALLENGE-001",
        "parent_protocol_sha256": p174["protocol_core_sha256"],
        "known_models": p174["models"],
        "unknown_models": [
            "copula_ordered_gaussian",
            "tempered_stable_0.8",
            "gaussian_cauchy_mixture",
            "nonshared_random_gain",
        ],
        "acceptance_rule": "unchanged P174 nearest-centroid and classwise thresholds",
        "records": 1600,
        "replicates": 180,
        "seed": SEED,
        "frozen_before_unknown_execution": True,
    }
    record = {"core": core, "canonical_core_sha256": canonical_digest(core)}
    OPEN_SET_PROTOCOL.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    reps, n, t2 = core["replicates"], core["records"], 4.0
    c1 = RNG.normal(size=(reps, n))
    c2 = RNG.normal(size=(reps, n)) * t2**0.5
    unknown_features = {}

    # Exact feature-invariance counterexample: ordering changes, multisets do not.
    base1 = RNG.normal(size=(reps, n))
    base2 = RNG.normal(size=(reps, n)) * t2**0.5
    base_features = quantile_features(base1, base2, c1, c2)
    ordered1 = np.sort(base1, axis=1)
    ordered2 = np.sort(base2, axis=1)
    ordered_features = quantile_features(ordered1, ordered2, c1, c2)
    unknown_features["copula_ordered_gaussian"] = ordered_features
    exact_order_blindness = float(np.max(np.abs(base_features - ordered_features)))

    x1 = np.clip(stable_rvs(0.8, (reps, n), RNG), -12, 12)
    x2 = np.clip(stable_rvs(0.8, (reps, n), RNG) * t2**1.25, -12, 12)
    unknown_features["tempered_stable_0.8"] = quantile_features(x1, x2, c1, c2)

    mask1 = RNG.random((reps, n)) < 0.1
    mask2 = RNG.random((reps, n)) < 0.1
    x1 = np.where(mask1, RNG.standard_cauchy((reps, n)), RNG.normal(size=(reps, n)))
    x2 = np.where(mask2, RNG.standard_cauchy((reps, n)), RNG.normal(size=(reps, n))) * t2**0.8
    unknown_features["gaussian_cauchy_mixture"] = quantile_features(x1, x2, c1, c2)

    gain = RNG.normal(0.7, 0.35, (reps, 1))
    x1 = RNG.normal(size=(reps, n))
    x2 = RNG.normal(size=(reps, n)) * t2 ** (0.5 + gain)
    unknown_features["nonshared_random_gain"] = quantile_features(x1, x2, c1, c2)

    rows = []
    known_names = p174["models"]
    for name, feature in unknown_features.items():
        labels, abstain = classify_features(feature, p174)
        counts = {known_names[i]: int(np.sum((labels == i) & (~abstain))) for i in range(len(known_names))}
        rows.append(
            {
                "unknown_model": name,
                "rejection_rate": float(np.mean(abstain)),
                "false_accept_rate": float(np.mean(~abstain)),
                "accepted_as": counts,
            }
        )

    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    ax.bar(
        [r["unknown_model"].replace("_", "\n") for r in rows],
        [r["rejection_rate"] for r in rows],
        color="#1F5A99",
    )
    ax.axhline(0.95, color="#A61B1B", ls="--", label="illustrative 95% rejection target")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("open-set rejection rate")
    ax.set_title("Program 185: closed-set abstention does not guarantee open-set rejection")
    ax.legend()
    fig.savefig(FIG / "program185_open_set.png", dpi=190)
    plt.close(fig)

    return {
        "protocol": OPEN_SET_PROTOCOL.name,
        "protocol_sha256": record["canonical_core_sha256"],
        "rows": rows,
        "exact_order_blindness_residual": exact_order_blindness,
        "impossibility_theorem": (
            "Every P174 feature is a function of one-time empirical multisets. "
            "Any reordering changes temporal dependence while preserving all "
            "features exactly; therefore the protocol cannot identify memory."
        ),
        "open_set_complete": False,
        "external_data_used": False,
        "status": "Frozen open-set challenge executed; feature-class incompleteness proven",
        "confidence": "Proven order-blindness; strong synthetic challenge",
    }


def program186_process_tensor_double_slit() -> dict:
    sigma = 0.8
    t1 = t2 = 1.0
    gamma1 = math.exp(-0.5 * sigma**2 * t1**2)
    gamma2 = math.exp(-0.5 * sigma**2 * t2**2)
    gamma_joint_static = math.exp(-0.5 * sigma**2 * (t1 + t2) ** 2)
    gamma_product = gamma1 * gamma2
    memory_witness = gamma_joint_static - gamma_product
    gamma_echo_static = 1.0
    gamma_echo_markov = gamma_product
    blur = 0.84
    path_rate = 0.3
    path_survival_contrast = math.exp(-2 * path_rate * (t1 + t2))

    models = {
        "quasistatic_memory": [gamma_joint_static, gamma_echo_static, 0.0],
        "markov_dephasing": [gamma_product, gamma_echo_markov, 0.0],
        "detector_blur": [blur, blur, 0.0],
        "path_diffusion": [0.0, 0.0, path_survival_contrast],
    }
    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    width = 0.2
    xpos = np.arange(len(models))
    for j, label in enumerate(["raw visibility", "echo visibility", "path survival contrast"]):
        ax.bar(xpos + (j - 1) * width, [v[j] for v in models.values()], width, label=label)
    ax.set_xticks(xpos, [k.replace("_", "\n") for k in models])
    ax.set_ylim(0, 1.08)
    ax.set_title("Program 186: intervention vector separates four loss mechanisms")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program186_process_tensor.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "TwoStepProcessInstrument",
            "components": [
                "path preparation",
                "two-step random-unitary process",
                "optional echo intervention",
                "path-population control",
                "phase POVM",
                "detector calibration",
                "record",
            ],
        },
        "quasistatic_one_step_coherence": gamma1,
        "quasistatic_two_step_coherence": gamma_joint_static,
        "markov_factorized_two_step": gamma_product,
        "nonfactorization_memory_witness": memory_witness,
        "echo_static": gamma_echo_static,
        "echo_markov": gamma_echo_markov,
        "intervention_vectors": models,
        "memory_present_for_quasistatic_process": abs(memory_witness) > 1e-12,
        "status": "Conditional two-step process instrument and memory witness constructed",
        "confidence": "Proven finite random-unitary formulas",
        "claim_boundary": (
            "The intervention vector separates the declared mechanisms; it is "
            "not a device-independent theorem or an external FIN prediction."
        ),
    }


def program187_environment_nonuniqueness() -> dict:
    gs = np.linspace(0, math.pi / 2, 101)
    gammas = np.cos(gs)
    visibilities = gammas.copy()

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.plot(gs, visibilities, lw=2)
    ax.set_xlabel("environment coupling g")
    ax.set_ylabel("output visibility gamma=cos(g)")
    ax.set_title("Program 187: the same system generator admits a continuum of environments")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program187_environment_nonuniqueness.png", dpi=190)
    plt.close(fig)

    sample_rows = []
    for g in (0.0, 0.3, 0.7, 1.1, math.pi / 2):
        gamma = math.cos(g)
        # Normalized Choi eigenvalues of the dephasing channel.
        choi_eigs = [(1 + gamma) / 2, (1 - gamma) / 2, 0.0, 0.0]
        sample_rows.append({"g": g, "gamma": gamma, "choi_eigenvalues": choi_eigs})

    return {
        "fixed_system_generator": "A=[[1,-1],[-1,1]]",
        "dilation_family": (
            "U_g=|0><0| tensor I + |1><1| tensor exp(-i g sigma_z), "
            "environment initially |+>."
        ),
        "reduced_channel": "off-diagonal rho_01 -> cos(g)*rho_01",
        "sample_rows": sample_rows,
        "same_A_distinct_channels": True,
        "uniqueness_no_go": (
            "A contains no value of g. Holding A fixed while varying g gives "
            "a continuum of CPTP dephasing channels and records. Therefore no "
            "unique environment or decoherence law is a function of A alone."
        ),
        "minimal_operational_replacement": (
            "EnvironmentProcessClass = reduced process tensor modulo "
            "operationally equivalent dilations."
        ),
        "status": "Environment-from-generator uniqueness refuted",
        "confidence": "Proven by explicit dilation family",
    }


def program188_phase_cocycle_provenance() -> dict:
    delta_omega = OMEGA_STRICT - OMEGA_LEGACY
    delta_phi = PHI_STRICT - PHI_LEGACY
    d = np.arange(24)
    legacy = np.exp(1j * (OMEGA_LEGACY * d + PHI_LEGACY))
    strict = np.exp(1j * (OMEGA_STRICT * d + PHI_STRICT))
    correction = strict / legacy

    fig, ax = plt.subplots(figsize=(7.2, 7.0), constrained_layout=True)
    ax.scatter(legacy.real, legacy.imag, label="legacy cyclotomic phase", s=35)
    ax.scatter(strict.real, strict.imag, label="strict infinite-order phase", s=35)
    ax.scatter(correction.real, correction.imag, label="correction cocycle", s=25)
    circle = np.linspace(0, 2 * math.pi, 400)
    ax.plot(np.cos(circle), np.sin(circle), color="black", alpha=0.25)
    ax.set_aspect("equal")
    ax.set_xlabel("real")
    ax.set_ylabel("imaginary")
    ax.set_title("Program 188: cyclotomic source versus transcendental target character")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program188_phase_cohomology.png", dpi=190)
    plt.close(fig)

    return {
        "cohomology": "H^1(Z,U(1))=Hom(Z,U(1))=U(1) for trivial coefficients",
        "legacy_generator": "exp(i*pi/4), algebraic root of unity",
        "strict_generator": "exp(i*743/4000), transcendental by Lindemann-Weierstrass",
        "correction_generator": "exp(i*(743/4000-pi/4)), transcendental",
        "offset_statement": (
            "delta_phi=13/80-pi/6 is an affine phase-origin correction, "
            "not part of a normalized group character."
        ),
        "algebraic_source_no_go": (
            "Finite algebraic operations on legacy cyclotomic values remain "
            "algebraic. They cannot produce the transcendental strict "
            "generator or correction value."
        ),
        "delta_omega": delta_omega,
        "delta_phi": delta_phi,
        "correction_residual": float(np.max(np.abs(legacy * correction - strict))),
        "strict_phase_source_exported": False,
        "status": "Cohomology class identified; algebraic legacy-source no-go proven",
        "confidence": "Proven",
        "claim_boundary": (
            "An analytic source law with genuinely new input could produce the "
            "strict character; the exact target quotient is not such a law."
        ),
    }


def program189_intrinsic_scale_obstruction() -> dict:
    A = np.asarray(
        [[1, -0.5, 0, -0.5], [-0.5, 1, -0.5, 0], [0, -0.5, 1, -0.5], [-0.5, 0, -0.5, 1]],
        dtype=float,
    )
    t = 0.7
    z = 0.4
    scales = np.logspace(-3, 3, 25)
    u_res = []
    p_res = []
    g_res = []
    raw_gap = []
    ratio_gap = []
    base_U = expm(-1j * t * A)
    base_P = expm(-t * A)
    base_G = np.linalg.inv(A + z * np.eye(4))
    base_eig = np.linalg.eigvalsh(A)
    for c in scales:
        Ac = c * A
        Uc = expm(-1j * (t / c) * Ac)
        Pc = expm(-(t / c) * Ac)
        Gc = np.linalg.inv(Ac + c * z * np.eye(4))
        eig = np.linalg.eigvalsh(Ac)
        u_res.append(float(np.linalg.norm(Uc - base_U)))
        p_res.append(float(np.linalg.norm(Pc - base_P)))
        g_res.append(float(np.linalg.norm(c * Gc - base_G)))
        raw_gap.append(float(eig[1]))
        ratio_gap.append(float(eig[1] / eig[-1]))

    fig, ax1 = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax1.loglog(scales, raw_gap, label="raw spectral gap", color="#D55E00")
    ax1.loglog(scales, ratio_gap, label="gap/max eigenvalue", color="#1F5A99")
    ax1.set_xlabel("scale action c")
    ax1.set_ylabel("spectral quantity")
    ax1.set_title("Program 189: scale orbit changes units but not normalized predictions")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program189_scale_orbit.png", dpi=190)
    plt.close(fig)

    return {
        "scale_action": "(A,t,z) -> (cA,t/c,cz)",
        "max_unitary_orbit_residual": max(u_res),
        "max_heat_orbit_residual": max(p_res),
        "max_resolvent_covariance_residual": max(g_res),
        "normalized_gap_variation": max(ratio_gap) - min(ratio_gap),
        "raw_gap_scale_ratio": max(raw_gap) / min(raw_gap),
        "no_section_theorem": (
            "All invariant spectral data are constant on the free R_{>0} "
            "orbit. An invariant rule cannot choose one orbit representative; "
            "a scale-charged datum or explicit calibration is necessary."
        ),
        "constructed_object": {
            "name": "SpectralClockFrame",
            "definition": "a chosen positive tau_* such that tau_* A is dimensionless",
            "transformation": "A->cA requires tau_*->tau_*/c",
            "status": "conditional calibration object, not strict source",
        },
        "full_physical_units": (
            "A single clock frame calibrates this generator. Independent "
            "length/action/mass dimensions still require the rank-three CA "
            "package or equivalent conversion data."
        ),
        "status": "Intrinsic one-generator scale obstruction proven",
        "confidence": "Proven",
    }


def program190_external_intake_gate() -> dict:
    required = [
        "public_source_identifier",
        "license_identifier",
        "immutable_hashes",
        "preparation_provenance",
        "raw_time_ordered_records",
        "physical_units_and_timestamps",
        "detector_calibration",
        "apparatus_memory_calibration",
        "reference_control",
        "exclusion_audit",
        "independent_analysis_boundary",
    ]

    raw_paths = list((ROOT / "raw_strain_unfiltered").glob("*.h5"))
    valid_hdf = 0
    for path in raw_paths:
        with path.open("rb") as handle:
            valid_hdf += handle.read(8) == b"\x89HDF\r\n\x1a\n"

    beta_manifest = ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/manifest_beta_channel.json"
    rebuild_manifest = ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_rebuild/manifest.json"
    beta = json.loads(beta_manifest.read_text(encoding="utf-8"))
    rebuild = json.loads(rebuild_manifest.read_text(encoding="utf-8"))
    beta_root = beta_manifest.parent
    beta_listed_present = sum((beta_root / row["path"]).exists() for row in beta["files"])
    rebuild_root = rebuild_manifest.parent
    rebuild_listed_present = sum((rebuild_root / row["path"]).exists() for row in rebuild["files"])

    bundles = [
        {
            "bundle": "raw_strain_unfiltered",
            "evidence": f"{len(raw_paths)} .h5 paths, {valid_hdf} valid HDF5 signatures, datasets lack embedded attributes",
            "passes": {
                "public_source_identifier": False,
                "license_identifier": False,
                "immutable_hashes": False,
                "preparation_provenance": False,
                "raw_time_ordered_records": valid_hdf > 0,
                "physical_units_and_timestamps": False,
                "detector_calibration": False,
                "apparatus_memory_calibration": False,
                "reference_control": False,
                "exclusion_audit": False,
                "independent_analysis_boundary": False,
            },
        },
        {
            "bundle": "beta_channel_true_external_v2",
            "evidence": f"manifest with {len(beta['files'])} hashes; {beta_listed_present} listed files present",
            "passes": {
                "public_source_identifier": True,
                "license_identifier": False,
                "immutable_hashes": True,
                "preparation_provenance": False,
                "raw_time_ordered_records": False,
                "physical_units_and_timestamps": False,
                "detector_calibration": False,
                "apparatus_memory_calibration": False,
                "reference_control": False,
                "exclusion_audit": False,
                "independent_analysis_boundary": False,
            },
        },
        {
            "bundle": "confirmatory_dataset_external_source_rebuild",
            "evidence": f"manifest with generic license text; {rebuild_listed_present} summary files present",
            "passes": {
                "public_source_identifier": False,
                "license_identifier": False,
                "immutable_hashes": True,
                "preparation_provenance": False,
                "raw_time_ordered_records": False,
                "physical_units_and_timestamps": False,
                "detector_calibration": False,
                "apparatus_memory_calibration": False,
                "reference_control": False,
                "exclusion_audit": False,
                "independent_analysis_boundary": False,
            },
        },
    ]
    matrix = np.asarray([[row["passes"][key] for key in required] for row in bundles], dtype=int)
    admitted = [row["bundle"] for row in bundles if all(row["passes"].values())]

    fig, ax = plt.subplots(figsize=(11.0, 4.6), constrained_layout=True)
    image = ax.imshow(matrix, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(required)), [x.replace("_", "\n") for x in required], fontsize=7)
    ax.set_yticks(range(len(bundles)), [x["bundle"].replace("_", "\n") for x in bundles], fontsize=8)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center", fontsize=8)
    ax.set_title("Program 190: external-data intake gate")
    fig.colorbar(image, ax=ax, ticks=[0, 1], shrink=0.75)
    fig.savefig(FIG / "program190_external_intake.png", dpi=190)
    plt.close(fig)

    return {
        "required_fields": required,
        "audited_bundles": bundles,
        "admitted_bundles": admitted,
        "external_protocol_executed": False,
        "intake_failure_record": len(admitted) == 0,
        "status": "Local external-candidate intake audited; no bundle admitted",
        "confidence": "Proven for the declared local schema audit",
        "claim_boundary": (
            "Some files have genuine external lineage, but none supplies the "
            "specific preparation/control/calibration record required by the "
            "Programs 172--185 operational protocol."
        ),
    }


def main() -> None:
    FIG.mkdir(exist_ok=True)
    previous = json.loads(PREVIOUS.read_text(encoding="utf-8"))
    results = {
        "metadata": {
            "title": "FIN Programs 178-190: Composition Laws, Process Memory, and Scale Obstructions",
            "release": "10.17",
            "version": "1.0.0",
            "date": "2026-07-27",
            "creator": "Żuchowski, Krzysztof",
            "affiliation": "Independent Researcher — Fractal Information Theory Project",
            "orcid": "0009-0002-0909-3613",
            "seed": SEED,
            "previous_results_sha256": sha256(PREVIOUS),
            "firecrawl_used": False,
            "external_web_used": False,
        },
        "programs": {
            "178": program178_tensor_intensive_state_laws(),
            "179": program179_positive_compression_source(),
            "180": program180_ball_certificate_assembly(previous),
            "181": program181_finite_dirichlet_core(),
            "182": program182_complete_reflection_convertibility(),
            "183": program183_block_robust_quantiles(),
            "184": program184_nonlinear_multi_control(),
            "185": program185_open_set_challenge(previous),
            "186": program186_process_tensor_double_slit(),
            "187": program187_environment_nonuniqueness(),
            "188": program188_phase_cocycle_provenance(),
            "189": program189_intrinsic_scale_obstruction(),
            "190": program190_external_intake_gate(),
        },
        "global_verdict": {
            "spectral_core_survives": True,
            "tensor_intensive_state_law_derived": False,
            "strict_positive_beta_value_derived": False,
            "full_ball_certificate_closed": False,
            "complete_reflection_conversion_object_constructed": True,
            "dependence_robust_block_theorem_constructed": True,
            "nonlinear_calibration_object_constructed": True,
            "process_memory_object_constructed": True,
            "environment_unique_from_A": False,
            "strict_phase_cocycle_sourced": False,
            "intrinsic_scale_selected": False,
            "external_dataset_admitted": False,
            "QW_2191_discharged": False,
            "legacy_strict_bridge_completed": False,
            "legacy_role_transfer_started": False,
            "L_total_or_ToE_claimed": False,
        },
        "recommended_next_programs": list(range(191, 204)),
    }
    RESULTS.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(RESULTS)
    print(OPEN_SET_PROTOCOL)


if __name__ == "__main__":
    main()
