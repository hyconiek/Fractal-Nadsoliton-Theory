#!/usr/bin/env python3
"""
QW-2091: Neutrino absolute-scale gate (deterministic, no-retune, no-scan).

Goal:
- derive (m_nu1, m_nu2, m_nu3) only when an absolute-scale observable is
  provided in anchor-free form,
- prevent false strict closure when only hierarchy assumptions are available.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "neutrino_absolute_scale_input_qw2091.json"
OUT_JSON = ROOT / "report_qw2091_neutrino_absolute_scale_gate.json"
OUT_MD = ROOT / "RAPORT_QW2091_NEUTRINO_ABSOLUTE_SCALE_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def masses_from_m1_normal(m1: float, dm21: float, dm31: float) -> Tuple[float, float, float]:
    m2 = math.sqrt(max(m1 * m1 + dm21, 0.0))
    m3 = math.sqrt(max(m1 * m1 + dm31, 0.0))
    return float(m1), float(m2), float(m3)


def solve_bisection(f, lo: float, hi: float, tol: float = 1e-12, nmax: int = 200) -> Tuple[bool, float]:
    flo = f(lo)
    fhi = f(hi)
    if not np.isfinite(flo) or not np.isfinite(fhi):
        return False, float("nan")
    if flo == 0.0:
        return True, lo
    if fhi == 0.0:
        return True, hi
    if flo * fhi > 0.0:
        return False, float("nan")
    a, b = lo, hi
    fa, fb = flo, fhi
    for _ in range(nmax):
        c = 0.5 * (a + b)
        fc = f(c)
        if not np.isfinite(fc):
            return False, float("nan")
        if abs(fc) <= tol or abs(b - a) <= tol:
            return True, float(c)
        if fa * fc <= 0.0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return False, float("nan")


def source_metadata_complete(in_data: Dict) -> bool:
    src = str(in_data.get("source", "")).strip()
    citation = str(in_data.get("citation", "")).strip()
    ref_url = str(in_data.get("reference_url", "")).strip()
    src_sha = str(in_data.get("source_sha256", "")).strip()
    src_ver = str(in_data.get("source_version", "")).strip()
    low = f"{src} {citation} {ref_url}".lower()
    not_placeholder = bool(src) and ("placeholder" not in low)
    has_reference = bool(citation or ref_url)
    has_integrity = bool(src_sha or src_ver)
    return bool(not_placeholder and has_reference and has_integrity)


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2091 neutrino absolute-scale gate")
    p.add_argument("--input", default=str(DEFAULT_INPUT), help="Absolute-scale input JSON.")
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    in_path = Path(args.input).resolve()

    input_present = in_path.exists()
    in_data = load_json(in_path) if input_present else {}

    # Defaults consistent with prior internal chain, but strict closure requires external absolute input.
    dm21 = float(in_data.get("dm21_ev2", 7.42e-5))
    dm31 = float(in_data.get("dm31_ev2", 2.517e-3))
    hierarchy = str(in_data.get("hierarchy", "normal")).lower().strip()
    mode = str(in_data.get("mode", "none")).lower().strip()

    # PMNS electron-row weights for m_beta mode.
    pmns_abs = None
    for row in reg.get("groups", {}).get("flavor", []):
        if row.get("id") == "pmns_matrix_abs":
            pmns_abs = np.array(row["value"], dtype=float)
            break
    if pmns_abs is None:
        raise RuntimeError("Missing pmns_matrix_abs in registry.")
    w1, w2, w3 = [float(x * x) for x in pmns_abs[0, :]]

    absolute_observable_present = False
    solver_converged = False
    fit_reason = "no_absolute_scale_observable"
    m1 = m2 = m3 = None

    if input_present and hierarchy == "normal":
        if mode == "sum_mnu" and in_data.get("sum_mnu_ev") is not None:
            target = float(in_data["sum_mnu_ev"])
            absolute_observable_present = True

            def f(x):
                a, b, c = masses_from_m1_normal(x, dm21, dm31)
                return (a + b + c) - target

            ok, x = solve_bisection(f, 0.0, 2.0, tol=1e-12, nmax=300)
            if ok:
                m1, m2, m3 = masses_from_m1_normal(float(x), dm21, dm31)
                solver_converged = True
                fit_reason = "sum_mnu_solution"
            else:
                fit_reason = "sum_mnu_solver_failed_or_out_of_range"

        elif mode == "m_beta" and in_data.get("m_beta_eff_ev") is not None:
            target = float(in_data["m_beta_eff_ev"])
            absolute_observable_present = True

            def f(x):
                a, b, c = masses_from_m1_normal(x, dm21, dm31)
                mb = math.sqrt(max(w1 * a * a + w2 * b * b + w3 * c * c, 0.0))
                return mb - target

            ok, x = solve_bisection(f, 0.0, 2.0, tol=1e-12, nmax=300)
            if ok:
                m1, m2, m3 = masses_from_m1_normal(float(x), dm21, dm31)
                solver_converged = True
                fit_reason = "m_beta_solution"
            else:
                fit_reason = "m_beta_solver_failed_or_out_of_range"
        else:
            fit_reason = "unsupported_mode_or_missing_observable"
    elif hierarchy != "normal":
        fit_reason = "unsupported_hierarchy_in_qw2091"

    # Fallback assumption path (explicitly non-closing).
    if not solver_converged:
        m1, m2, m3 = masses_from_m1_normal(0.0, dm21, dm31)

    sum_mnu = float(m1 + m2 + m3)
    m_beta_eff = float(math.sqrt(max(w1 * m1 * m1 + w2 * m2 * m2 + w3 * m3 * m3, 0.0)))

    sum_bound = in_data.get("sum_mnu_upper_bound_ev")
    if sum_bound is None:
        sum_bound_ok = True
    else:
        sum_bound_ok = bool(sum_mnu <= float(sum_bound))

    flags = {
        "absolute_observable_present": bool(absolute_observable_present),
        "source_metadata_complete": bool(source_metadata_complete(in_data) if input_present else False),
        "hierarchy_supported": bool(hierarchy == "normal"),
        "solver_converged": bool(solver_converged),
        "masses_finite": bool(np.isfinite(m1) and np.isfinite(m2) and np.isfinite(m3)),
        "ordered_m1_le_m2_le_m3": bool(m1 <= m2 <= m3),
        "sum_bound_ok_if_provided": bool(sum_bound_ok),
        "no_anchor_feedback_loop": bool(in_data.get("provenance_anchor_free", False)),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_pass = bool(all(flags.values()))

    if strict_pass:
        verdict = "NEUTRINO_ABSOLUTE_SCALE_GATE_PASS_STRICT"
        updates = [
            {
                "id": "m_nu1",
                "predicted_value": float(m1),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": f"qw2091_absolute_scale_from_{mode}",
                "notes": "Absolute neutrino scale solved with independent absolute observable in deterministic mode.",
            },
            {
                "id": "m_nu2",
                "predicted_value": float(m2),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": f"qw2091_absolute_scale_from_{mode}",
                "notes": "Derived jointly with m_nu1/m_nu3 from strict absolute-scale gate.",
            },
            {
                "id": "m_nu3",
                "predicted_value": float(m3),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": f"qw2091_absolute_scale_from_{mode}",
                "notes": "Derived jointly with m_nu1/m_nu2 from strict absolute-scale gate.",
            },
        ]
    else:
        verdict = "NEUTRINO_ABSOLUTE_SCALE_GATE_PENDING_NONCLOSING"
        updates = [
            {
                "id": "m_nu1",
                "predicted_value": float(m1),
                "status": "derived_model_assumption_nonclosing",
                "strict_level": "model_assumption_anchor",
                "method": "qw2091_normal_hierarchy_assumption_fallback",
                "notes": "Absolute neutrino scale remains underdetermined without independent absolute observable.",
            },
            {
                "id": "m_nu2",
                "predicted_value": float(m2),
                "status": "derived_model_assumption_nonclosing",
                "strict_level": "model_assumption_anchor",
                "method": "qw2091_normal_hierarchy_assumption_fallback",
                "notes": "Depends on hierarchy/oscillation assumptions in fallback mode.",
            },
            {
                "id": "m_nu3",
                "predicted_value": float(m3),
                "status": "derived_model_assumption_nonclosing",
                "strict_level": "model_assumption_anchor",
                "method": "qw2091_normal_hierarchy_assumption_fallback",
                "notes": "Depends on hierarchy/oscillation assumptions in fallback mode.",
            },
        ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "input_json": str(in_path) if input_present else None,
        },
        "inputs": {
            "input_present": input_present,
            "mode": mode,
            "hierarchy": hierarchy,
            "fit_reason": fit_reason,
        },
        "derived": {
            "m_nu1_ev": float(m1),
            "m_nu2_ev": float(m2),
            "m_nu3_ev": float(m3),
            "sum_mnu_ev": sum_mnu,
            "m_beta_eff_ev": m_beta_eff,
            "dm21_ev2": dm21,
            "dm31_ev2": dm31,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "updates": updates,
        "verdict": verdict,
        "required_next_step": (
            "INTEGRATE_QW2091_IN_QW2069_AND_RERUN_FULL_CLOSURE_CHAIN"
            if strict_pass
            else "PROVIDE_INDEPENDENT_ABSOLUTE_NEUTRINO_SCALE_OBSERVABLE_AND_RERUN_QW2091"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2091: NEUTRINO ABSOLUTE SCALE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- input_present: {input_present}",
        f"- mode: {mode}",
        f"- fit_reason: {fit_reason}",
        "",
        "## Derived Masses",
        f"- m_nu1: {m1:.12e} eV",
        f"- m_nu2: {m2:.12e} eV",
        f"- m_nu3: {m3:.12e} eV",
        f"- sum_mnu: {sum_mnu:.12e} eV",
        f"- m_beta_eff: {m_beta_eff:.12e} eV",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2091] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2091] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2091] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
