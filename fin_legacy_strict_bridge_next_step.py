#!/usr/bin/env python3
"""Next proof-computational legacy-to-strict bridge audit.

This script constructs the smallest honest bridge-obstruction/witness matrix for
newly reconstructed K*_legacy versus K_strict_gate. It avoids claiming closure:
positive scalar/envelope maps are separated from the missing phase/frequency and
source objects.
"""

from __future__ import annotations

import json
import math
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer
from reportlab.lib.units import cm


ROOT = Path(__file__).resolve().parent
OUT_MD = ROOT / "FIN_Legacy_Strict_Bridge_Next_Proof_Audit.md"
OUT_JSON = ROOT / "FIN_Legacy_Strict_Bridge_Next_Proof_Audit_Results.json"
OUT_PDF = ROOT / "FIN_Legacy_Strict_Bridge_Next_Proof_Audit.pdf"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA_L = math.pi / 4.0
PHI_L = math.pi / 6.0
OMEGA_S = 0.18575
PHI_S = 0.16250
BETA_S = 1.0
ETA_S = 1.8
D = np.arange(1, 13, dtype=float)


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def k_rej(d: np.ndarray) -> np.ndarray:
    return np.exp(-2.9 * d) * ((1.0 + 0.2 * d) / (1.0 + d)) * np.cos(OMEGA_L * d + PHI_L)


def k_star(d: np.ndarray, A: float, beta: float) -> np.ndarray:
    return A * np.cos(OMEGA_L * d + PHI_L) / (1.0 + beta * d)


def k_strict(d: np.ndarray) -> np.ndarray:
    return np.cos(OMEGA_S * d + PHI_S) / (1.0 + BETA_S * (d ** ETA_S))


def normalize(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v.copy()


def best_scalar_fit(source: np.ndarray, target: np.ndarray, nonnegative: bool = False) -> dict[str, Any]:
    denom = float(source @ source)
    c = float((source @ target) / denom) if denom else 0.0
    if nonnegative and c < 0.0:
        c = 0.0
    pred = c * source
    rel = float(np.linalg.norm(pred - target) / np.linalg.norm(target))
    cos = float((source @ target) / (np.linalg.norm(source) * np.linalg.norm(target)))
    return {"scalar": c, "relative_l2": rel, "cosine": cos, "nonnegative_constrained": nonnegative}


def grid_fit_kstar_to_strict() -> dict[str, Any]:
    rows = []
    beta_grid = np.concatenate([
        np.linspace(0.001, 0.12, 240),
        np.linspace(0.12, 2.0, 260),
    ])
    target = k_strict(D)
    for beta in beta_grid:
        base = k_star(D, 1.0, float(beta))
        fit = best_scalar_fit(base, target)
        rows.append({"beta": float(beta), **fit})
    best = min(rows, key=lambda r: r["relative_l2"])
    fixed = {}
    for beta in [0.01, 0.05, 0.1, 1.0]:
        base = k_star(D, 1.0, beta)
        fixed[str(beta)] = best_scalar_fit(base, target) | {"beta": beta, "positive_fit": best_scalar_fit(base, target, nonnegative=True)}
    return {"grid_size": len(rows), "best": best, "fixed_beta_fits": fixed}


def sign_obstruction() -> dict[str, Any]:
    target = k_strict(D)
    legacy = k_star(D, 1.0, 0.05)
    product = k_rej(D)
    s_target = np.sign(target).astype(int).tolist()
    s_legacy = np.sign(legacy).astype(int).tolist()
    s_product = np.sign(product).astype(int).tolist()
    mism_legacy = [int(d) for d, a, b in zip(D, s_legacy, s_target) if a != b]
    mism_product = [int(d) for d, a, b in zip(D, s_product, s_target) if a != b]
    theorem = (
        "For every A>0 and every positive radial envelope E(d)>0, "
        "sign(A*E(d)*cos(pi*d/4+pi/6)) = sign(cos(pi*d/4+pi/6)); "
        "therefore positive amplitude/envelope completion cannot map the K* sign row to the strict sign row."
    )
    return {
        "d": D.astype(int).tolist(),
        "strict_signs": s_target,
        "kstar_signs": s_legacy,
        "product_signs": s_product,
        "kstar_vs_strict_mismatch_distances": mism_legacy,
        "product_vs_strict_mismatch_distances": mism_product,
        "positive_scalar_envelope_no_go": theorem,
    }


def damping_witness(beta_l: float) -> dict[str, Any]:
    # After amplitude and phase/frequency have been factored away, this is the
    # exact positive envelope map from linear legacy damping to strict power damping.
    Dfac = (1.0 + beta_l * D) / (1.0 + D ** ETA_S)
    ratios = Dfac.tolist()
    diffs = np.diff(Dfac)
    return {
        "beta_legacy": beta_l,
        "formula": "D_beta(d)=(1+beta_legacy*d)/(1+d^1.8)",
        "values_d1_12": ratios,
        "positive_all_d1_12": bool(np.all(Dfac > 0)),
        "monotone_decreasing_d1_12": bool(np.all(diffs < 0)),
        "d1": float(Dfac[0]),
        "d12": float(Dfac[-1]),
        "tail_exponent": "~ beta_legacy*d / d^1.8 = beta_legacy*d^-0.8 for beta_legacy>0",
        "meaning": "exact envelope witness only after amplitude and phase/frequency maps are supplied; not a source theorem for eta=1.8",
    }


def assembly_quotient(beta_l: float, A: float) -> dict[str, Any]:
    theta_l = OMEGA_L * D + PHI_L
    theta_s = OMEGA_S * D + PHI_S
    amp = 1.0 / A
    phase = np.cos(theta_s) / np.cos(theta_l)
    damp = (1.0 + beta_l * D) / (1.0 + D ** ETA_S)
    q = amp * phase * damp
    lhs = k_star(D, A, beta_l) * q
    err = lhs - k_strict(D)
    return {
        "A": A,
        "beta_legacy": beta_l,
        "amplitude_factor": amp,
        "phase_factor_values_d1_12": phase.tolist(),
        "damping_factor_values_d1_12": damp.tolist(),
        "q_values_d1_12": q.tolist(),
        "max_abs_identity_error_d1_12": float(np.max(np.abs(err))),
        "phase_factor_sign_changes": int(np.sum(np.sign(phase[1:]) != np.sign(phase[:-1]))),
        "phase_factor_min": float(np.min(phase)),
        "phase_factor_max": float(np.max(phase)),
        "status": "finite algebraic assembly witness; not a strict source/selector/role-transfer theorem",
    }


def amplitude_audit() -> dict[str, Any]:
    target = k_strict(D)
    a1 = k_star(D, 1.0, 0.05)
    a29 = k_star(D, 2.9, 0.05)
    a4ln2 = k_star(D, ALPHA_GEO, 0.05)
    return {
        "A_interpretation": "A is a scalar coupling/measure-density normalization multiplying the whole legacy profile; it is not a shape parameter.",
        "A_drops_from_normalized_shape": bool(np.allclose(normalize(a1), normalize(a29)) and np.allclose(normalize(a1), normalize(a4ln2))),
        "A_needed_for": [
            "absolute unnormalized coupling strength",
            "legacy role formulas that use alpha_geo explicitly",
            "kernel moments before quotient/row normalization",
            "any future action-density normalization theorem",
        ],
        "A_not_needed_for": [
            "sign pattern",
            "normalized profile shape",
            "correlation/cosine similarity",
            "row-stochastic Markov operator after full row normalization",
        ],
        "A_equals_1_verdict": "If one works only in quotient-normalized shape space, A=1 is a valid gauge. If one claims alpha_geo roles or absolute physical couplings, A cannot be deleted; it becomes an open source/normalization obligation.",
        "fits_beta_0_05": {
            "A_1": best_scalar_fit(a1, target),
            "A_1_positive": best_scalar_fit(a1, target, nonnegative=True),
            "A_2_9": best_scalar_fit(a29, target),
            "A_2_9_positive": best_scalar_fit(a29, target, nonnegative=True),
            "A_4ln2": best_scalar_fit(a4ln2, target),
            "A_4ln2_positive": best_scalar_fit(a4ln2, target, nonnegative=True),
        },
        "alpha_geo_identity": {
            "value": ALPHA_GEO,
            "exp_value": math.exp(ALPHA_GEO),
            "bits": ALPHA_GEO / math.log(2.0),
            "status": "Shannon 16-state/four-bit identity; not physical unit or selector source",
        },
    }


def compare_rejected_vs_kstar() -> dict[str, Any]:
    target = k_strict(D)
    product = k_rej(D)
    star = k_star(D, 1.0, 0.05)
    return {
        "product_best_scalar_to_strict": best_scalar_fit(product, target),
        "product_positive_scalar_to_strict": best_scalar_fit(product, target, nonnegative=True),
        "kstar_beta0_05_best_scalar_to_strict": best_scalar_fit(star, target),
        "kstar_beta0_05_positive_scalar_to_strict": best_scalar_fit(star, target, nonnegative=True),
        "product_norm": float(np.linalg.norm(product)),
        "kstar_beta0_05_norm": float(np.linalg.norm(star)),
        "strict_norm": float(np.linalg.norm(target)),
        "product_abs_lt_1e_minus_3_for_d_ge_2": bool(np.all(np.abs(product[D >= 2]) < 1e-3)),
        "bridge_easier_verdict": "K* is easier than K_rej for bridge work because it preserves a non-dead long-range hyperbolic profile and separates damping from phase; it is not close to strict without a nontrivial phase/frequency map.",
    }


def theoretical_object_matrix() -> list[dict[str, Any]]:
    return [
        {
            "object": "O_A scalar amplitude/measure normalization",
            "constructed_here": "quotient factor A^{-1}; proof that A drops from normalized shape",
            "status": "witness/gauge decomposition only",
            "missing_for_closure": "source theorem fixing A as alpha_geo or proving amplitude absorption is role-safe",
        },
        {
            "object": "O_D damping/compression map",
            "constructed_here": "D_beta(d)=(1+beta_l*d)/(1+d^1.8), positive and decreasing on d=1..12",
            "status": "finite/exact envelope witness after phase separation",
            "missing_for_closure": "dynamical source of eta=1.8 and target-independent beta/Z_beta theorem",
        },
        {
            "object": "O_P phase/frequency transport",
            "constructed_here": "P(d)=cos(theta_s)/cos(theta_l) finite quotient on d=1..12",
            "status": "algebraic quotient with sign changes, not source law",
            "missing_for_closure": "phase/frequency/topological source localizer and selector-safe transport theorem",
        },
        {
            "object": "O_S selector/orientation source",
            "constructed_here": "positive scalar/envelope sign no-go matrix",
            "status": "obstruction: QW-2191 remains open",
            "missing_for_closure": "non-premise orientation/polarity source plus coupling-polarity theorem",
        },
        {
            "object": "O_U unit/action-density bridge",
            "constructed_here": "none; imported no-go status from guardrails/P2760-P2764/P3167-P3168",
            "status": "no-go on current artifacts",
            "missing_for_closure": "canonical reference-cell/action-density theorem with physical units",
        },
    ]


def build_results() -> dict[str, Any]:
    beta_l = 0.05
    return {
        "metadata": {
            "title": "FIN Legacy-Strict Bridge Next Proof Audit",
            "date_utc": datetime.now(timezone.utc).isoformat(),
            "git_head_at_generation": git_head(),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "prior_artifact_anchors": {
            "P2760": "seven open gaps: ontology-to-kernel measure, amplitude normalization, phase/frequency/topological source, damping/compression bridge, moments-to-couplings provenance, Lagrangian reverse closure, stale flags",
            "P2691": "alpha_geo=4 ln 2 and scalar-shape normalization real; no role-safe amplitude absorption source",
            "P2692": "positive beta orbit/gauge representability real; no target-independent positive beta/Z_beta source theorem",
            "P3167_P3168": "no strict unit source S+ on current monomial/state map artifacts",
            "scratch_bridge_assembly": "existing finite assembly quotient maps legacy to strict on audited nodes but explicitly does not export source, selector, role-transfer, or ToE closure",
        },
        "sign_obstruction": sign_obstruction(),
        "damping_witness_beta0_05": damping_witness(beta_l),
        "assembly_quotient_beta0_05_A_2_9": assembly_quotient(beta_l, 2.9),
        "amplitude_audit": amplitude_audit(),
        "grid_fit_kstar_to_strict": grid_fit_kstar_to_strict(),
        "rejected_vs_kstar": compare_rejected_vs_kstar(),
        "theoretical_object_matrix": theoretical_object_matrix(),
        "next_honest_step": {
            "recommendation": "Do not run another generic bridge replay. Attack exactly O_P: construct or refute a phase/frequency transport source localizer for P(d)=cos(theta_s)/cos(theta_l), while conditioning O_A and O_D as quotient witnesses only.",
            "why": "The new K* already repairs the historical envelope, so the remaining hard obstruction is not scalar A and not positive damping alone; it is the sign-changing phase/frequency/source layer plus the still-open unit/selector obligations.",
            "success_criterion": "A machine-checkable theorem or finite exhaustive certificate that P(d) comes from a strict exported source object and is selector-safe; otherwise emit a bounded no-go for this exact P object.",
        },
    }


def md_table(rows: list[list[Any]]) -> str:
    header = "| " + " | ".join(map(str, rows[0])) + " |"
    sep = "| " + " | ".join(["---"] * len(rows[0])) + " |"
    body = ["| " + " | ".join(str(x).replace("\n", " ") for x in row) + " |" for row in rows[1:]]
    return "\n".join([header, sep] + body)


def build_markdown(r: dict[str, Any]) -> str:
    lines: list[str] = []
    m = r["metadata"]
    sign = r["sign_obstruction"]
    amp = r["amplitude_audit"]
    grid = r["grid_fit_kstar_to_strict"]
    rej = r["rejected_vs_kstar"]
    damp = r["damping_witness_beta0_05"]
    asm = r["assembly_quotient_beta0_05_A_2_9"]
    lines += [
        "# FIN Legacy–Strict Bridge Next Proof Audit",
        "",
        "## Metadata",
        "",
        f"- Date UTC: `{m['date_utc']}`",
        f"- Git HEAD at generation: `{m['git_head_at_generation']}`",
        f"- Python: `{m['python']}`",
        f"- NumPy: `{m['numpy']}`",
        "",
        "## Executive verdict",
        "",
        "This is a narrower proof-computational follow-up to the legacy-star audit. It constructs the explicit scalar/envelope/phase quotient objects between the reconstructed `K*` legacy class and `K_strict_gate`, then separates finite algebraic witnesses from missing source theorems.",
        "",
        "- **[Proven finite obstruction]** Positive scalar and positive damping/envelope maps cannot fix the `K*` to strict sign mismatch; a phase/frequency transport object is mathematically necessary.",
        "- **[Proven quotient fact]** `A` is a scalar multiplier: it drops out of normalized shape/sign comparisons, but it is still needed for unnormalized coupling strength, legacy `alpha_geo` roles, moment magnitudes, and any future action-density normalization.",
        "- **[Strong evidence]** `K*` is easier to bridge than the rejected double-damped product because it is not exponentially dead and admits a clean positive damping quotient. It is still not close to strict without a sign-changing phase quotient.",
        "- **[No-go on current artifacts]** The constructed quotient `Q(d)` is an algebraic finite assembly witness, not a source theorem, selector theorem, role-transfer theorem, unit theorem, or ToE closure.",
        "",
        "## Prior repo anchors used",
        "",
    ]
    lines.append(md_table([["Anchor", "Imported status"]] + [[k, v] for k, v in r["prior_artifact_anchors"].items()]))
    lines += [
        "",
        "## Kernel definitions compared",
        "",
        "```text",
        "K*_legacy(d;A,beta_l) = A*cos(pi*d/4 + pi/6)/(1+beta_l*d)",
        "K_strict_gate(d)      = cos(0.18575*d + 0.16250)/(1+d^1.8)",
        "K_rej(d)              = exp(-2.9*d)*(1+0.2*d)/(1+d)*cos(pi*d/4+pi/6)",
        "```",
        "",
        "## Sign obstruction theorem",
        "",
        sign["positive_scalar_envelope_no_go"],
        "",
    ]
    lines.append(md_table([
        ["Row", "Signs d=1..12"],
        ["strict", sign["strict_signs"]],
        ["K*", sign["kstar_signs"]],
        ["K_rej", sign["product_signs"]],
        ["K* mismatch distances", sign["kstar_vs_strict_mismatch_distances"]],
    ]))
    lines += [
        "",
        "**Consequence:** the bridge cannot be only `A` plus a positive damping/compression map. The missing theoretical object is a phase/frequency transport/source, not another scalar amplitude fit.",
        "",
        "## Damping/compression quotient object",
        "",
        f"For `beta_l=0.05`, define `{damp['formula']}`. On `d=1..12` it is positive: `{damp['positive_all_d1_12']}`, monotone decreasing: `{damp['monotone_decreasing_d1_12']}`, with `D(1)={damp['d1']:.6g}` and `D(12)={damp['d12']:.6g}`.",
        "",
        "This is useful because it cleanly separates linear legacy damping from strict nonlinear compression after phase and amplitude have been factored away. It does **not** derive `eta=1.8`; it only says what positive envelope quotient would be needed if the phase/source layer were supplied.",
        "",
        "## Finite assembly quotient",
        "",
        "After separating amplitude, phase, and damping, the exact finite quotient is:",
        "",
        "```text",
        "Q(d) = A^-1 * [cos(theta_strict(d))/cos(theta_legacy(d))] * [(1+beta_l*d)/(1+d^1.8)]",
        "K*_legacy(d;A,beta_l) * Q(d) = K_strict_gate(d)",
        "```",
        "",
        f"For `A=2.9`, `beta_l=0.05`, the maximum absolute identity error on `d=1..12` is `{asm['max_abs_identity_error_d1_12']:.3e}`. This is an exact finite algebraic witness, but the phase factor ranges from `{asm['phase_factor_min']:.6g}` to `{asm['phase_factor_max']:.6g}` and has `{asm['phase_factor_sign_changes']}` sign changes, so it is not a positive dissipative bridge.",
        "",
        "## What is `A`?",
        "",
        amp["A_interpretation"],
        "",
    ]
    lines.append(md_table([
        ["Question", "Answer"],
        ["Does A drop from normalized shape?", amp["A_drops_from_normalized_shape"]],
        ["If A=1, is it needed?", amp["A_equals_1_verdict"]],
        ["A is needed for", amp["A_needed_for"]],
        ["A is not needed for", amp["A_not_needed_for"]],
    ]))
    lines += [
        "",
        f"The value `A=4 ln 2={amp['alpha_geo_identity']['value']:.15f}` remains **[Proven]** as `ln 16`, the entropy of a uniform four-bit/16-state cell, but **[No-go]** as a physical unit, selector, or role-transfer source on current artifacts.",
        "",
        "## Is the new legacy kernel easier to bridge to strict?",
        "",
    ]
    lines.append(md_table([
        ["Comparison", "relative L2", "cosine", "scalar"],
        ["unconstrained scalar K_rej -> strict", f"{rej['product_best_scalar_to_strict']['relative_l2']:.6g}", f"{rej['product_best_scalar_to_strict']['cosine']:.6g}", f"{rej['product_best_scalar_to_strict']['scalar']:.6g}"],
        ["positive scalar K_rej -> strict", f"{rej['product_positive_scalar_to_strict']['relative_l2']:.6g}", f"{rej['product_positive_scalar_to_strict']['cosine']:.6g}", f"{rej['product_positive_scalar_to_strict']['scalar']:.6g}"],
        ["unconstrained scalar K*(beta=0.05) -> strict", f"{rej['kstar_beta0_05_best_scalar_to_strict']['relative_l2']:.6g}", f"{rej['kstar_beta0_05_best_scalar_to_strict']['cosine']:.6g}", f"{rej['kstar_beta0_05_best_scalar_to_strict']['scalar']:.6g}"],
        ["positive scalar K*(beta=0.05) -> strict", f"{rej['kstar_beta0_05_positive_scalar_to_strict']['relative_l2']:.6g}", f"{rej['kstar_beta0_05_positive_scalar_to_strict']['cosine']:.6g}", f"{rej['kstar_beta0_05_positive_scalar_to_strict']['scalar']:.6g}"],
        ["best grid beta K* -> strict", f"{grid['best']['relative_l2']:.6g}", f"{grid['best']['cosine']:.6g}", f"beta={grid['best']['beta']:.6g}, A_fit={grid['best']['scalar']:.6g}"],
    ]))
    lines += [
        "",
        "The scalar table is intentionally harsh: because the strict row is positive until d=7 while K* has negative entries at d=2..5, a positive scalar-only fit to K* collapses to zero. The rejected product can be inflated by a huge scalar to match mostly the first entry, but this is not a bridge; its norm is tiny and it is exponentially dead for d>=2. K* is therefore easier only after allowing the necessary phase/frequency object O_P plus the clean damping quotient O_D.",
        "",
        rej["bridge_easier_verdict"],
        "",
        "## Constructed theoretical-object matrix",
        "",
    ]
    lines.append(md_table([["Object", "Constructed here", "Status", "Missing for closure"]] + [[x["object"], x["constructed_here"], x["status"], x["missing_for_closure"]] for x in r["theoretical_object_matrix"]]))
    lines += [
        "",
        "## Recommended next honest step",
        "",
        f"**Recommendation:** {r['next_honest_step']['recommendation']}",
        "",
        f"**Why:** {r['next_honest_step']['why']}",
        "",
        f"**Success criterion:** {r['next_honest_step']['success_criterion']}",
        "",
        "Do not promote this audit to strict bridge closure, selector closure, role transfer, physical units, `L_total`, or ToE.",
    ]
    return "\n".join(lines) + "\n"


def write_pdf(md: str) -> None:
    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(str(OUT_PDF), pagesize=A4, leftMargin=1.5*cm, rightMargin=1.5*cm, topMargin=1.5*cm, bottomMargin=1.5*cm)
    story = []
    code = False
    for raw in md.splitlines():
        if raw.startswith("```"):
            code = not code
            continue
        if not raw.strip():
            story.append(Spacer(1, 0.14*cm))
            continue
        line = raw.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        if raw.startswith("# "):
            story.append(Paragraph(line[2:], styles["Heading1"]))
        elif raw.startswith("## "):
            story.append(Paragraph(line[3:], styles["Heading2"]))
        else:
            story.append(Paragraph(line, styles["Code"] if code or raw.startswith("| ") else styles["BodyText"]))
    doc.build(story)


def main() -> None:
    results = build_results()
    OUT_JSON.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    md = build_markdown(results)
    OUT_MD.write_text(md, encoding="utf-8")
    write_pdf(md)
    print(f"Wrote {OUT_MD.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PDF.name}")


if __name__ == "__main__":
    main()
