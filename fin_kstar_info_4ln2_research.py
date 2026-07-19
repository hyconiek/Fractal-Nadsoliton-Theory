#!/usr/bin/env python3
"""Research suite for the information-theoretic legacy kernel freeze:

    K*_info(d; β) = 4 ln 2 * cos(π d / 4 + π / 6) / (1 + β d)

Compares dual dynamics, Shannon diagnostics, and leaf-cut status against
historical freezes and the strict profile. All quantities are dimensionless.
No SI units, selector closure, role transfer, or ToE claims are exported.

Outputs
-------
FIN_KStar_Info_4ln2_Results.json
FIN_KStar_Info_4ln2_Figures/*.png
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
FIG = OUT / "FIN_KStar_Info_4ln2_Figures"
FIG.mkdir(exist_ok=True)

N = 12
PI = math.pi
LN2 = math.log(2.0)
ALPHA_GEO = 4.0 * LN2  # 4 ln 2 = ln 16
OMEGA_L = PI / 4.0
PHI_L = PI / 6.0
OMEGA_S, PHI_S, BETA_S, ETA_S = 0.18575, 0.16250, 1.0, 1.8
T_STAR = 1.5525
SEED = 20260719


def k_info(d, beta: float, A: float = ALPHA_GEO) -> np.ndarray:
    dd = np.asarray(d, dtype=float)
    out = np.zeros_like(dd, dtype=float)
    m = dd != 0
    out[m] = A * np.cos(OMEGA_L * dd[m] + PHI_L) / (1.0 + beta * dd[m])
    return out


def k_hist(d) -> np.ndarray:
    return k_info(d, beta=0.05, A=2.9)


def k_z12(d) -> np.ndarray:
    return k_info(d, beta=1.0, A=1.0)


def k_class(d) -> np.ndarray:
    return k_info(d, beta=0.01, A=ALPHA_GEO)


def k_strict(d) -> np.ndarray:
    dd = np.asarray(d, dtype=float)
    out = np.zeros_like(dd, dtype=float)
    m = dd != 0
    out[m] = np.cos(OMEGA_S * dd[m] + PHI_S) / (1.0 + BETA_S * dd[m] ** ETA_S)
    return out


def cyclic_distance(i: int, j: int, n: int = N) -> int:
    return min(abs(i - j), n - abs(i - j))


def build_W(kernel_fn, use_abs: bool = True) -> np.ndarray:
    W = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            w = float(np.asarray(kernel_fn(cyclic_distance(i, j))).reshape(-1)[0])
            W[i, j] = abs(w) if use_abs else w
    return W


def operator_from_W(W: np.ndarray) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    s = float(W[0].sum())
    A = s * np.eye(N) - W
    evals, evecs = np.linalg.eigh(A)
    evals = evals.copy()
    evals[np.abs(evals) < 1e-13] = 0.0
    return A, s, evals, evecs


def shannon_nats(p: np.ndarray) -> float:
    p = np.clip(np.real(p), 0.0, None)
    tot = p.sum()
    if tot <= 0:
        return 0.0
    p = p / tot
    nz = p > 1e-300
    return float(-np.sum(p[nz] * np.log(p[nz])))


def shannon_bits(p: np.ndarray) -> float:
    return shannon_nats(p) / LN2


def relative_info_nats(p: np.ndarray) -> float:
    return math.log(N) - shannon_nats(p)


def dual_audit(name: str, kernel_fn, use_abs: bool = True) -> dict:
    W = build_W(kernel_fn, use_abs=use_abs)
    # if signed and not all nonnegative, report both
    n_neg = int(np.sum(W < -1e-15)) if not use_abs else 0
    if use_abs:
        n_neg_raw = int(
            np.sum(build_W(kernel_fn, use_abs=False) < -1e-15)
        )
    else:
        n_neg_raw = n_neg
    A, s, evals, evecs = operator_from_W(W)
    gap = float(np.sort(evals)[1])
    dt = 1e-4
    U = (evecs * np.exp(-1j * dt * evals)) @ evecs.T
    P = (evecs * np.exp(-dt * evals)) @ evecs.T
    u_esc = float(1.0 - abs(U[0, 0]) ** 2)
    m_esc = float(1.0 - np.real(P[0, 0]))
    weights = np.array([W[0, d] for d in range(1, 7)], dtype=float)
    # capacity metrics
    cap_l1 = float(np.sum(np.abs(weights)))  # half-row positive mass / distance class once each dir counted in s
    # entropy of normalized |weights| on d=1..6 (with multiplicity 2 for d=1..5, 1 for d=6)
    mult = np.array([2, 2, 2, 2, 2, 1], dtype=float)
    mass = np.abs(weights) * mult
    mass = mass / max(mass.sum(), 1e-300)
    H_w = float(-np.sum(mass[mass > 0] * np.log(mass[mass > 0])))
    # inversion even residual
    inv = np.zeros((N, N))
    for i in range(N):
        inv[i, (-i) % N] = 1.0
    inv_res = float(np.linalg.norm(inv @ A @ inv.T - A))
    return {
        "name": name,
        "use_abs": use_abs,
        "n_negative_raw_weights": n_neg_raw,
        "weights_d1_to_d6": weights.tolist(),
        "row_sum_s": s,
        "spectral_gap": gap,
        "A_eigenvalues": evals.tolist(),
        "A_min_eigenvalue": float(evals.min()),
        "unitary_escape_over_dt": u_esc / dt,
        "markov_escape_over_dt": m_esc / dt,
        "short_time_unitary_loglog_proxy": math.log(max(u_esc, 1e-300)) / math.log(dt),
        "short_time_markov_loglog_proxy": math.log(max(m_esc, 1e-300)) / math.log(dt),
        "capacity_l1_row_half": cap_l1,
        "weight_distribution_entropy_nats": H_w,
        "weight_distribution_entropy_bits": H_w / LN2,
        "inversion_even_residual": inv_res,
        "evals": evals,
        "evecs": evecs,
        "W": W,
        "A": A,
        "s": s,
    }


def double_slit(evecs, evals, t: float, phase: float = 0.0, coherence: float = 1.0) -> dict:
    slit_a, slit_b = (-2) % N, 2
    U = (evecs * np.exp(-1j * t * evals)) @ evecs.T
    ua, ub = U[:, slit_a], U[:, slit_b]
    incoh = 0.5 * (np.abs(ua) ** 2 + np.abs(ub) ** 2)
    cross = np.real(np.exp(-1j * phase) * ua * np.conj(ub))
    coh = incoh + coherence * cross
    # heat from equal mixture of slits
    entrance = np.zeros(N)
    entrance[[slit_a, slit_b]] = 0.5
    P = (evecs * np.exp(-t * evals)) @ evecs.T
    diff = np.real(P @ entrance)
    # visibility at site 0 over phase
    phases = np.linspace(0, 2 * PI, 361)
    p0 = []
    for ph in phases:
        cr = np.real(np.exp(-1j * ph) * ua * np.conj(ub))
        p0.append(float((incoh + coherence * cr)[0]))
    pmax, pmin = max(p0), min(p0)
    vis = (pmax - pmin) / (pmax + pmin) if (pmax + pmin) > 1e-15 else 0.0
    H_coh_nats = shannon_nats(coh)
    H_inc_nats = shannon_nats(incoh)
    H_diff_nats = shannon_nats(diff)
    # entropy reduction due to interference (nats): positive if coherent more peaked
    dH = H_inc_nats - H_coh_nats
    # proxy mutual information: I ≈ H(incoh) - H(coh) when coherence=1 vs mixture
    # also 0.5 * sum |cross| total interference strength
    strength = 0.5 * float(np.sum(np.abs(cross)))
    return {
        "t": t,
        "slit_a": slit_a,
        "slit_b": slit_b,
        "coherent": coh,
        "incoherent": incoh,
        "diffusive": diff,
        "cross": cross,
        "visibility_site0": vis,
        "p_max_site0": pmax,
        "p_min_site0": pmin,
        "interference_strength": strength,
        "H_coherent_nats": H_coh_nats,
        "H_coherent_bits": H_coh_nats / LN2,
        "H_incoherent_nats": H_inc_nats,
        "H_incoherent_bits": H_inc_nats / LN2,
        "H_diffusive_nats": H_diff_nats,
        "H_diffusive_bits": H_diff_nats / LN2,
        "entropy_reduction_nats": dH,
        "entropy_reduction_bits": dH / LN2,
        "proxy_mutual_info_nats": max(dH, 0.0),
        "proxy_mutual_info_bits": max(dH, 0.0) / LN2,
        "coherent_sum": float(np.sum(coh)),
        "relative_info_coherent_nats": relative_info_nats(coh),
        "relative_info_incoherent_nats": relative_info_nats(incoh),
    }


def markov_info_flow(evecs, evals, tmax: float = 5.0, npts: int = 201) -> dict:
    times = np.linspace(0.0, tmax, npts)
    p0 = np.zeros(N)
    p0[0] = 1.0
    I = []
    H = []
    for t in times:
        P = (evecs * np.exp(-t * evals)) @ evecs.T
        p = np.real(P @ p0)
        I.append(relative_info_nats(p))
        H.append(shannon_nats(p))
    I = np.asarray(I, dtype=float)
    H = np.asarray(H, dtype=float)
    dI = np.diff(I)
    return {
        "times": times.tolist(),
        "relative_info_nats": I.tolist(),
        "shannon_nats": H.tolist(),
        "monotone": bool(np.all(dI <= 1e-10)),
        "backflow_steps": int(np.sum(dI > 1e-10)),
        "I0": float(I[0]),
        "I_final": float(I[-1]),
    }


def constants_proxy_table(A: float = ALPHA_GEO, betas: list[float] | None = None) -> dict:
    """Legacy-style algebraic proxies only — not strict role transfer.

    From historical FIN notes (model-level, not promoted):
      α_EM^{-1} ≈ (1/2) * (A/β) * (1-β)   with β = β_tors
      sin^2 θ_W = ω/π = 1/4
      ħ_eff ≈ π^3  (geometry; independent of A)
    """
    if betas is None:
        betas = [0.01, 0.05, 0.1]
    rows = []
    for b in betas:
        if b <= 0 or b >= 1:
            alpha_inv = float("nan")
        else:
            alpha_inv = 0.5 * (A / b) * (1.0 - b)
        rows.append(
            {
                "beta": b,
                "A": A,
                "A_over_beta": A / b,
                "alpha_EM_inv_legacy_proxy": alpha_inv,
                "sin2_theta_W_omega_over_pi": OMEGA_L / PI,
                "hbar_eff_pi3_proxy": PI**3,
                "alpha_geo_equals_4ln2": abs(A - ALPHA_GEO) < 1e-15,
                "note": "model-level legacy identities only; not strict role transfer",
            }
        )
    return {
        "alpha_geo": A,
        "alpha_geo_bits": A / LN2,
        "exp_alpha_geo": math.exp(A),
        "A_phi_phase_area": 2 * PI / A,
        "rows": rows,
        "disclaimer": (
            "These are historical algebraic proxies conditioned on legacy role "
            "formulas. They do not discharge units leaf-cut or ToE closure."
        ),
    }


def information_motivation() -> dict:
    return {
        "alpha_geo": ALPHA_GEO,
        "alpha_geo_high_precision": f"{ALPHA_GEO:.18f}",
        "identity": "4 ln 2 = ln 16 = H(uniform on 16 outcomes)",
        "bits": 4.0,
        "nats": ALPHA_GEO,
        "interpretations": [
            "Shannon entropy of a uniform 4-bit / 16-state source: H = ln(16) = 4 ln 2",
            "Natural amplitude for an informational cell of capacity 4 bits (octave-related 2^4)",
            "Matches frozen legacy ontological amplitude α_geo in current FIN packages",
            "Phase-area section A_φ = 2π / α_geo is dimensionless S/ħ-shaped bookkeeping",
        ],
        "not_claimed": [
            "Does not by itself export SI units or ħ",
            "Does not select orientation / discharge QW-2191",
            "Does not force β or strict (ω,φ,η)",
        ],
        "confidence": "Proven as scalar Shannon identity; Speculative as unique physical amplitude law without further theorem",
    }


def profile_comparison() -> dict:
    d = np.arange(1, 13, dtype=float)
    profiles = {
        "info_b0.01": k_info(d, 0.01),
        "info_b0.05": k_info(d, 0.05),
        "info_b0.10": k_info(d, 0.10),
        "hist_A2.9_b0.05": k_hist(d),
        "Z12_A1_b1": k_z12(d),
        "classical_legacy": k_class(d),
        "strict": k_strict(d),
    }
    # l2 distances between abs-normalized profiles
    keys = list(profiles.keys())
    dist = {}
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            va = np.abs(profiles[a])
            vb = np.abs(profiles[b])
            va = va / max(np.linalg.norm(va), 1e-300)
            vb = vb / max(np.linalg.norm(vb), 1e-300)
            dist[f"{a}__vs__{b}"] = {
                "l2_abs_normalized": float(np.linalg.norm(va - vb)),
                "cosine_abs": float(np.dot(va, vb)),
                "sign_agreement_count": int(np.sum(np.sign(profiles[a]) == np.sign(profiles[b]))),
            }
    return {
        "d": d.tolist(),
        "profiles": {k: v.tolist() for k, v in profiles.items()},
        "pairwise_abs_normalized": dist,
    }


def make_figures(audits: dict, slit_pack: dict, flows: dict, profiles: dict) -> None:
    d = np.array(profiles["d"])
    # 1) kernel profiles
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    styles = {
        "info_b0.01": ("o-", "#1f5a99", r"$K^*_{\mathrm{info}},\beta=0.01$"),
        "info_b0.05": ("s-", "#19733a", r"$K^*_{\mathrm{info}},\beta=0.05$"),
        "info_b0.10": ("^-", "#6a3d9a", r"$K^*_{\mathrm{info}},\beta=0.10$"),
        "hist_A2.9_b0.05": ("d--", "#d55e00", r"$K^*_{\mathrm{hist}}$"),
        "Z12_A1_b1": ("x--", "#888888", r"$K^*_{\mathrm{Z12}}$"),
        "classical_legacy": (".-.", "#444444", r"$K_\ell$ classical"),
        "strict": ("D-", "#a61b1b", r"$K_s$ strict"),
    }
    for key, (sty, col, lab) in styles.items():
        ax.plot(d, profiles["profiles"][key], sty, color=col, lw=1.7, label=lab, ms=5)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("distance $d$")
    ax.set_ylabel("kernel value")
    ax.set_title(r"Kernel profiles: $K^*_{\mathrm{info}}$ vs freezes")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "kernel_profiles.png", dpi=200)
    plt.close(fig)

    # 2) double slit for best info beta (0.01 classical-like) and 0.05
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    sites = np.arange(N)
    for ax, beta, title in [
        (axes[0], 0.01, r"$\beta=0.01$, $t=1.5525$"),
        (axes[1], 0.05, r"$\beta=0.05$, $t=1.5525$"),
    ]:
        key = f"info_b{beta:.2f}"
        ds = slit_pack[key]
        ax.plot(sites, ds["coherent"], "o-", color="#1f5a99", label="coherent", lw=1.8)
        ax.plot(sites, ds["incoherent"], "s--", color="#d55e00", label="incoherent", lw=1.5)
        ax.plot(sites, ds["diffusive"], "^-.", color="#19733a", label="Markov", lw=1.5)
        ax.set_title(title)
        ax.set_xlabel("site")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("probability")
    axes[0].legend(fontsize=8)
    fig.suptitle(r"Double-slit protocol for $K^*_{\mathrm{info}}=4\ln2\cdot\cos/(1+\beta d)$")
    fig.tight_layout()
    fig.savefig(FIG / "double_slit.png", dpi=200)
    plt.close(fig)

    # 3) entropy / relative info flow
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    for beta, col in [(0.01, "#1f5a99"), (0.05, "#19733a"), (0.10, "#6a3d9a")]:
        key = f"info_b{beta:.2f}"
        fl = flows[key]
        ax.plot(fl["times"], fl["relative_info_nats"], color=col, lw=2.0, label=rf"$\beta={beta}$")
    # strict comparison
    fls = flows["strict"]
    ax.plot(fls["times"], fls["relative_info_nats"], color="#a61b1b", lw=2.0, ls="--", label="strict")
    ax.set_xlabel("dimensionless time $t$")
    ax.set_ylabel(r"relative information $\log 12 - H$ (nats)")
    ax.set_title(r"Markov relative-information flow for $K^*_{\mathrm{info}}$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "entropy_flow.png", dpi=200)
    plt.close(fig)

    # 4) spectral comparison bar
    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    names = []
    gaps = []
    ss = []
    for key in ["info_b0.01", "info_b0.05", "info_b0.10", "hist", "Z12", "classical", "strict"]:
        a = audits[key]
        names.append(a["name"])
        gaps.append(a["spectral_gap"])
        ss.append(a["row_sum_s"])
    x = np.arange(len(names))
    w = 0.38
    ax.bar(x - w / 2, ss, w, label="row sum $s$", color="#1f5a99")
    ax.bar(x + w / 2, gaps, w, label=r"gap $\gamma$", color="#d55e00")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("value")
    ax.set_title(r"Dual-dynamics scales: $s$ and spectral gap")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "spectral_scales.png", dpi=200)
    plt.close(fig)


def sanitize(obj):
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items() if k not in {"evals", "evecs", "W", "A"}}
    if isinstance(obj, list):
        return [sanitize(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return sanitize(obj.tolist())
    return obj


def main() -> None:
    motivation = information_motivation()
    constants = constants_proxy_table()
    profiles = profile_comparison()

    # dual audits
    audit_specs = [
        ("info_b0.01", lambda d: k_info(d, 0.01), "K*_info β=0.01"),
        ("info_b0.05", lambda d: k_info(d, 0.05), "K*_info β=0.05"),
        ("info_b0.10", lambda d: k_info(d, 0.10), "K*_info β=0.10"),
        ("hist", k_hist, "K*_hist A=2.9 β=0.05"),
        ("Z12", k_z12, "K*_Z12 A=1 β=1"),
        ("classical", k_class, "classical legacy 4ln2 β=0.01"),
        ("strict", k_strict, "strict gate"),
    ]
    audits = {}
    for key, fn, label in audit_specs:
        # strict already nonnegative on d=1..6; still use abs for uniformity
        use_abs = True
        a = dual_audit(label, fn, use_abs=use_abs)
        audits[key] = a

    # double slit at t*
    slit_pack = {}
    for key in ["info_b0.01", "info_b0.05", "info_b0.10", "hist", "Z12", "classical", "strict"]:
        a = audits[key]
        ds = double_slit(a["evecs"], a["evals"], T_STAR)
        # store without large arrays in nested copy later
        slit_pack[key] = {
            **{k: v for k, v in ds.items() if k not in {"coherent", "incoherent", "diffusive", "cross"}},
            "coherent": ds["coherent"].tolist(),
            "incoherent": ds["incoherent"].tolist(),
            "diffusive": ds["diffusive"].tolist(),
        }

    # markov flows
    flows = {}
    for key in ["info_b0.01", "info_b0.05", "info_b0.10", "strict", "classical"]:
        a = audits[key]
        flows[key] = markov_info_flow(a["evecs"], a["evals"])

    make_figures(audits, slit_pack, flows, profiles)

    # β ranking for dual balance: prefer gap not too large (local only) nor tiny;
    # score = cosine_to_strict_abs * markov_monotone * visibility, penalize extreme s
    strict_w = np.abs(np.array(audits["strict"]["weights_d1_to_d6"]))
    strict_w = strict_w / max(np.linalg.norm(strict_w), 1e-300)
    beta_rank = []
    for beta, key in [(0.01, "info_b0.01"), (0.05, "info_b0.05"), (0.10, "info_b0.10")]:
        w = np.abs(np.array(audits[key]["weights_d1_to_d6"]))
        w = w / max(np.linalg.norm(w), 1e-300)
        cos = float(np.dot(w, strict_w))
        vis = slit_pack[key]["visibility_site0"]
        mono = 1.0 if flows[key]["monotone"] else 0.0
        # balance: spectral gap near O(1) like strict ~0.75
        gap = audits[key]["spectral_gap"]
        gap_score = math.exp(-abs(math.log(max(gap, 1e-12) / 0.754)))
        score = 0.35 * cos + 0.25 * vis + 0.20 * mono + 0.20 * gap_score
        beta_rank.append(
            {
                "beta": beta,
                "cosine_to_strict_abs_weights": cos,
                "visibility_tstar": vis,
                "markov_info_monotone": bool(mono),
                "spectral_gap": gap,
                "row_sum_s": audits[key]["row_sum_s"],
                "gap_score": gap_score,
                "composite_score": score,
                "entropy_reduction_bits_tstar": slit_pack[key]["entropy_reduction_bits"],
                "proxy_MI_bits_tstar": slit_pack[key]["proxy_mutual_info_bits"],
            }
        )
    beta_rank.sort(key=lambda r: r["composite_score"], reverse=True)

    leaf_cuts = {
        "dimensionlessness_units": {
            "status": "still open",
            "label": "Refuted as closed by 4 ln 2 alone",
            "reason": (
                "4 ln 2 is a dimensionless Shannon cell amplitude. It improves "
                "information-theoretic naturalness of A but supplies no weight-one "
                "S_+, ħ_*, or (ℓ_*, τ_*) conversion package."
            ),
        },
        "selector_QW2191": {
            "status": "still open",
            "label": "Refuted as closed",
            "reason": "Radial K*_info Laplacians remain inversion-even (residual 0).",
        },
        "symmetry_breaking": {
            "status": "still open as strict source",
            "label": "Refuted as non-premise source",
            "reason": "cos(πd/4+π/6) is phase structure, not an origin/polarity theorem.",
        },
    }

    verdict = {
        "is_better_info_theoretic_legacy_candidate": {
            "answer": "Yes, as amplitude choice within the path-transformed class",
            "confidence": "Strong evidence",
            "detail": (
                "A=4 ln 2 is the unique Shannon entropy of a uniform 4-bit cell "
                "and matches the frozen legacy α_geo. Relative to A=2.9 or A=1 it "
                "has clearer information motivation; relative to the rejected "
                "double-damped product it is algebraically superior."
            ),
        },
        "best_beta_dual_stability_balance": {
            "recommended_beta": beta_rank[0]["beta"],
            "ranking": beta_rank,
            "confidence": "Moderate evidence",
            "detail": (
                "Composite score balances abs-weight cosine to strict, visibility, "
                "Markov info monotonicity, and gap proximity to the strict scale. "
                "β=0.01 recovers classical legacy and strongest long-range mass; "
                "β=0.05 is diagram-historical; β=0.10 is more local."
            ),
        },
        "helps_close_leaf_cuts": {
            "units": False,
            "selector": False,
            "symmetry_breaking": False,
            "confidence": "Proven (negative) under current radial/invariant audits",
            "detail": (
                "4 ln 2 closes none of the three leaf-cuts. It only naturalizes the "
                "dimensionless amplitude A inside W0 informational layer."
            ),
        },
        "overall": {
            "status": "K*_info accepted as preferred information-theoretic freeze of the K* class",
            "formula": "K*_info(d)=4 ln 2 * cos(πd/4+π/6)/(1+β d)",
            "preferred_beta_for_legacy_continuity": 0.01,
            "preferred_beta_for_diagram_history": 0.05,
            "do_not_promote": [
                "ToE",
                "strict identification",
                "role transfer without bridge theorem",
                "physical units from α_geo alone",
            ],
        },
    }

    report = {
        "metadata": {
            "title": "K*_info research: A=4 ln 2 information-theoretic legacy freeze",
            "date": "2026-07-19",
            "seed": SEED,
            "context": "Programs 41-50 + Release 10.4 supplement",
            "numpy": np.__version__,
        },
        "frozen_object": {
            "formula": "K*_info(d) = 4 ln 2 * cos(π d / 4 + π / 6) / (1 + β d)",
            "A": ALPHA_GEO,
            "A_high_precision": f"{ALPHA_GEO:.18f}",
            "betas_scanned": [0.01, 0.05, 0.1],
            "omega": "π/4",
            "phi": "π/6",
        },
        "information_motivation": motivation,
        "constants_proxy": constants,
        "profile_comparison": profiles,
        "dual_audits": {k: sanitize(v) for k, v in audits.items()},
        "double_slit_tstar": sanitize(slit_pack),
        "markov_info_flows": sanitize(flows),
        "leaf_cuts": leaf_cuts,
        "verdict": verdict,
    }

    out = OUT / "FIN_KStar_Info_4ln2_Results.json"
    with out.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    # console summary
    print("=" * 72)
    print("K*_info research  A = 4 ln 2 =", f"{ALPHA_GEO:.15f}")
    print("=" * 72)
    print("\nMotivation: H(uniform 16) = ln 16 = 4 ln 2 (Proven Shannon identity)")
    print("\nDual audits (abs-repaired):")
    print(f"{'name':28s} {'s':>12s} {'gap':>12s} {'U_esc/dt':>12s} {'M_esc/dt':>12s}")
    for key in ["info_b0.01", "info_b0.05", "info_b0.10", "hist", "Z12", "classical", "strict"]:
        a = audits[key]
        print(
            f"{a['name'][:28]:28s} {a['row_sum_s']:12.6f} {a['spectral_gap']:12.6f} "
            f"{a['unitary_escape_over_dt']:12.4e} {a['markov_escape_over_dt']:12.6f}"
        )
    print("\nDouble-slit at t*=1.5525:")
    for key in ["info_b0.01", "info_b0.05", "info_b0.10", "strict"]:
        ds = slit_pack[key]
        print(
            f"  {key:12s} vis={ds['visibility_site0']:.6f}  "
            f"H_coh={ds['H_coherent_bits']:.6f} bit  "
            f"H_inc={ds['H_incoherent_bits']:.6f} bit  "
            f"ΔH={ds['entropy_reduction_bits']:.6f} bit  "
            f"I~={ds['proxy_mutual_info_bits']:.6f} bit"
        )
    print("\nβ ranking:")
    for r in beta_rank:
        print(
            f"  β={r['beta']:.2f} score={r['composite_score']:.4f} "
            f"cos={r['cosine_to_strict_abs_weights']:.4f} gap={r['spectral_gap']:.4f}"
        )
    print("\nLeaf-cuts: units/selector/symmetry all still open")
    print("Verdict:", verdict["overall"]["status"])
    print("Wrote", out)
    print("Figures in", FIG)


if __name__ == "__main__":
    main()
