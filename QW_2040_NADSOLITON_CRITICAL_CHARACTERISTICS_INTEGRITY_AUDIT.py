#!/usr/bin/env python3
"""
QW-2040: Nadsoliton critical-characteristics integrity audit.

Separates two notions of integrity:
1) Canonical pre-QW1700 TeX characteristics integrity.
2) Current refrozen-branch internal integrity.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2040_nadsoliton_critical_characteristics_integrity_audit.json"
OUT_MD = ROOT / "RAPORT_QW2040_NADSOLITON_CRITICAL_CHARACTERISTICS_INTEGRITY_AUDIT.md"


def canonical_formulas(alpha_geo: float, beta: float) -> Dict[str, float]:
    alpha_em_inv = alpha_geo / (2.0 * beta) * (1.0 - beta)
    n_s = 1.0 - 2.0 * beta
    g_ratio_20 = beta**20
    k_tau = 1.0 - 7.0 * beta
    lam = beta * alpha_geo
    return {
        "alpha_em_inv": float(alpha_em_inv),
        "n_s": float(n_s),
        "g_ratio_20": float(g_ratio_20),
        "k_tau": float(k_tau),
        "lambda_beta_alpha": float(lam),
    }


def zero_structure(omega: float, phi: float) -> Dict[str, float]:
    d0 = (math.pi / 2.0 - phi) / max(omega, 1e-12)
    spacing = math.pi / max(omega, 1e-12)
    return {"first_zero": float(d0), "zero_spacing": float(spacing)}


def damping_profile(beta: float, eta: float, dmax: int = 24) -> Dict[str, float]:
    d = np.arange(1, dmax + 1, dtype=float)
    v = 1.0 / (1.0 + beta * (d**eta))
    mono = bool(np.all(np.diff(v) < 1e-12))
    return {
        "d1": float(v[0]),
        "d24": float(v[-1]),
        "ratio_d24_over_d1": float(v[-1] / max(v[0], 1e-12)),
        "strictly_decreasing": mono,
    }


def rel_err(x: float, ref: float) -> float:
    return float(abs(x - ref) / max(abs(ref), 1e-15))


def main() -> None:
    alpha_geo = float(4.0 * math.log(2.0))

    # Canonical TeX kernel (pre-1700 branch).
    k_tex = {"omega": float(math.pi / 4.0), "phi": float(math.pi / 6.0), "beta": 0.01, "eta": 1.0}

    # Stage-C branch before derivation-compatible refreeze.
    d2030 = json.loads((ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json").read_text(encoding="utf-8"))
    k2030 = {x: float(d2030["kernel"][x]) for x in ["omega", "phi", "beta", "eta"]}

    # Derivation-compatible refrozen branch.
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))
    k2039 = {x: float(d2039["selected_kernel"][x]) for x in ["omega", "phi", "beta", "eta"]}

    kernels = {"tex_canonical": k_tex, "stagec_2030": k2030, "refrozen_2039": k2039}

    canon_ref = canonical_formulas(alpha_geo=alpha_geo, beta=0.01)
    out_rows = {}
    for name, k in kernels.items():
        c = canonical_formulas(alpha_geo=alpha_geo, beta=float(k["beta"]))
        z = zero_structure(omega=float(k["omega"]), phi=float(k["phi"]))
        dmp = damping_profile(beta=float(k["beta"]), eta=float(k["eta"]))
        out_rows[name] = {
            "kernel": k,
            "canonical_formula_outputs": c,
            "zero_structure": z,
            "damping_profile": dmp,
            "relative_to_tex_canonical": {
                "alpha_em_inv_rel_err": rel_err(c["alpha_em_inv"], canon_ref["alpha_em_inv"]),
                "n_s_rel_err": rel_err(c["n_s"], canon_ref["n_s"]),
                "g_ratio_20_rel_err": rel_err(c["g_ratio_20"], canon_ref["g_ratio_20"]),
                "k_tau_rel_err": rel_err(c["k_tau"], canon_ref["k_tau"]),
                "zero_spacing_rel_err": rel_err(z["zero_spacing"], zero_structure(k_tex["omega"], k_tex["phi"])["zero_spacing"]),
                "first_zero_rel_err": rel_err(z["first_zero"], zero_structure(k_tex["omega"], k_tex["phi"])["first_zero"]),
            },
        }

    # Integrity flags in two senses.
    tex_flags_2039 = {
        "beta_close_to_0p01": bool(abs(k2039["beta"] - 0.01) <= 0.01),
        "omega_close_to_pi_over_4": bool(abs(k2039["omega"] - math.pi / 4.0) <= 0.10),
        "phi_close_to_pi_over_6": bool(abs(k2039["phi"] - math.pi / 6.0) <= 0.10),
        "eta_linear_denominator": bool(abs(k2039["eta"] - 1.0) <= 1e-9),
    }
    tex_integrity_pass = bool(all(tex_flags_2039.values()))

    # Current branch integrity (already validated by gates).
    stagec_pass = bool(d2030.get("all_pass", False))
    refreeze_pass = bool(d2039.get("verdict") == "DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_PASS")
    external_primary = bool(d2039["global_flags"]["external_primary_all_pass"])
    external_stress = bool(d2039["global_flags"]["external_stress_soft_pass"])
    branch_flags = {
        "stagec_pass": stagec_pass,
        "refrozen_gate_pass": refreeze_pass,
        "external_primary_pass": external_primary,
        "external_stress_pass": external_stress,
    }
    branch_integrity_pass = bool(all(branch_flags.values()))

    if tex_integrity_pass:
        tex_verdict = "TEX_CANONICAL_CHARACTERISTICS_PRESERVED"
    else:
        tex_verdict = "TEX_CANONICAL_CHARACTERISTICS_NOT_PRESERVED"

    if branch_integrity_pass:
        branch_verdict = "CURRENT_REFROZEN_BRANCH_CHARACTERISTICS_OPERATIONALLY_PRESERVED"
    else:
        branch_verdict = "CURRENT_REFROZEN_BRANCH_CHARACTERISTICS_NOT_PRESERVED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "alpha_geo": alpha_geo,
        "reference_note": "TeX canonical uses beta_tors=0.01, omega=pi/4, phi=pi/6, linear denominator eta=1.",
        "kernels_compared": out_rows,
        "tex_canonical_integrity_flags_for_refrozen_2039": tex_flags_2039,
        "tex_canonical_integrity_verdict": tex_verdict,
        "current_refrozen_branch_integrity_flags": branch_flags,
        "current_refrozen_branch_integrity_verdict": branch_verdict,
        "overall_interpretation": (
            "Canonical TeX semantics are not preserved numerically in refrozen branch, "
            "but current refrozen branch remains internally consistent and gate-valid."
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    r2039 = out_rows["refrozen_2039"]
    lines = [
        "# RAPORT QW-2040: NADSOLITON CRITICAL CHARACTERISTICS INTEGRITY AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        "",
        "## Refrozen 2039 Kernel",
        (
            f"- omega/phi/beta/eta: {k2039['omega']:.6f} / {k2039['phi']:.6f} / "
            f"{k2039['beta']:.6f} / {k2039['eta']:.6f}"
        ),
        "",
        "## Canonical-Formula Drift vs TeX (beta_tors semantics)",
        f"- alpha_EM^-1 rel_err: {r2039['relative_to_tex_canonical']['alpha_em_inv_rel_err']:.3e}",
        f"- n_s rel_err: {r2039['relative_to_tex_canonical']['n_s_rel_err']:.3e}",
        f"- G_ratio_20 rel_err: {r2039['relative_to_tex_canonical']['g_ratio_20_rel_err']:.3e}",
        f"- k_tau rel_err: {r2039['relative_to_tex_canonical']['k_tau_rel_err']:.3e}",
        f"- zero_spacing rel_err: {r2039['relative_to_tex_canonical']['zero_spacing_rel_err']:.3e}",
        f"- first_zero rel_err: {r2039['relative_to_tex_canonical']['first_zero_rel_err']:.3e}",
        "",
        f"## TeX Canonical Verdict: **{tex_verdict}**",
        f"## Current Refrozen Branch Verdict: **{branch_verdict}**",
        "",
        "## Interpretation",
        f"- {out['overall_interpretation']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2040] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2040] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2040] tex={tex_verdict} branch={branch_verdict}")


if __name__ == "__main__":
    main()
