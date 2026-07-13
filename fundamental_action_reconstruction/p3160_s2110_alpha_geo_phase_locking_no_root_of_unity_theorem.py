"""P3160/S2110: alpha_geo phase-locking no-root-of-unity theorem.

This is the stricter follow-up to P3159.  P3159 showed that alpha_geo=4 ln 2
is compatible with a dimensionless 2*pi phase-area section.  P3160 tests the
stronger and more dangerous interpretation: whether alpha_geo itself can be a
canonical rational phase slot, i.e. alpha_geo/(2*pi) in Q, so that
exp(i*alpha_geo) is a root of unity.

The exact proof obligation is symbolic: if alpha_geo/(2*pi)=p/q with p!=0,
then 16=e^(alpha_geo)=e^(2*pi*p/q).  By Gelfond-Schneider in the standard form
exp(pi*r) is transcendental for every nonzero algebraic r; hence the right side
is transcendental, contradicting the algebraicity of 16.  Therefore no exact
Z_N phase locking exists for any finite N.  This complements, rather than
replaces, the P3159 dimensionless A_phi section.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.json"
MD = GEN / "p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P3159": GEN / "p3159_s2109_alpha_geo_pi_phase_compatibility_audit.json",
    "P3111": GEN / "p3111_s2061_symplectic_phase_area_section_audit.json",
    "alpha_geo": GEN / "alpha_geo_strict_derived_v1.json",
}
ALPHA_GEO = 4.0 * math.log(2.0)  # ln(16)
TWO_PI = 2.0 * math.pi
RATIO = ALPHA_GEO / TWO_PI
A_PHI = TWO_PI / ALPHA_GEO


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def rg_hits() -> dict[str, Any]:
    patterns = {
        "phase_locking": r"root.of.unity|root-of-unity|phase.*slot|Z12 slot|phase locking|alpha_geo_over_2pi|alpha_geo/\(2\*pi\)",
        "alpha_pi_section": r"alpha_geo|A_phi|2\*pi/alpha_geo|pi/\(2 ln 2\)|phase-area",
        "closed_lanes": r"selector closure|QW-2191|bridge completion|role transfer|L_total|ToE|Planck|apparatus|physical.*unit",
    }
    out: dict[str, Any] = {}
    for name, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-g", "*.py", "-g", "*.md", "-g", "*.json"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        lines = [line for line in proc.stdout.splitlines() if line]
        out[name] = {"count": len(lines), "samples": lines[:20]}
    return out


def denominator_scan(max_den: int = 144) -> list[dict[str, Any]]:
    rows = []
    best_by_den = None
    for n in range(1, max_den + 1):
        nearest = round(n * RATIO)
        residual_cycles = n * RATIO - nearest
        residual_phase = TWO_PI * residual_cycles
        row = {
            "denominator_N": n,
            "nearest_integer_k": nearest,
            "k_over_N": nearest / n,
            "cycle_residual": residual_cycles,
            "phase_residual_mod_2pi": residual_phase,
            "absolute_phase_residual": abs(residual_phase),
            "exact_locking_excluded_by_theorem": True,
        }
        rows.append(row)
        if best_by_den is None or row["absolute_phase_residual"] < best_by_den["absolute_phase_residual"]:
            best_by_den = row
    return rows


def continued_fraction_ladder() -> list[dict[str, Any]]:
    rows = []
    for max_den in [12, 24, 48, 96, 144, 288, 576, 1152, 2304, 4608, 9216]:
        q = Fraction(RATIO).limit_denominator(max_den)
        phase_error = abs(TWO_PI * (RATIO - q.numerator / q.denominator))
        rows.append({
            "max_denominator": max_den,
            "p": q.numerator,
            "q": q.denominator,
            "p_over_q": q.numerator / q.denominator,
            "ratio_error_abs": abs(RATIO - q.numerator / q.denominator),
            "phase_error_abs": phase_error,
            "exact_equality_excluded_by_theorem": True,
        })
    return rows


def theorem_steps() -> list[dict[str, Any]]:
    return [
        {"step": 1, "claim": "Assume alpha_geo/(2*pi)=p/q with integers p,q and p!=0.", "status": "assumption_for_contradiction"},
        {"step": 2, "claim": "Since alpha_geo=4 ln 2=ln(16), exponentiation gives 16=e^(2*pi*p/q).", "status": "algebraic_rewrite"},
        {"step": 3, "claim": "For nonzero rational r=2p/q, e^(pi*r) is transcendental by the Gelfond-Schneider consequence e^(pi*algebraic_nonzero).", "status": "transcendence_input"},
        {"step": 4, "claim": "Thus e^(2*pi*p/q) is transcendental, contradicting 16 being algebraic.", "status": "contradiction"},
        {"step": 5, "claim": "Therefore alpha_geo/(2*pi) is irrational and exp(i*alpha_geo) is not a root of unity.", "status": "proved"},
    ]


def object_matrix() -> list[dict[str, Any]]:
    return [
        {"object": "A_phi_phase_area_section", "formula": "A_phi=2*pi/alpha_geo", "constructed": True, "accepted_as": "dimensionless normalization", "closure_exported": False, "blocker": "does not select a unit, origin, or physical scale"},
        {"object": "R_alpha_phase_ratio", "formula": "alpha_geo/(2*pi)", "constructed": True, "accepted_as": "irrational phase-winding ratio", "closure_exported": False, "blocker": "irrationality prevents exact finite Z_N slot source"},
        {"object": "U_alpha_root_of_unity_lock", "formula": "exp(i*alpha_geo)^N=1", "constructed": True, "accepted_as": "exact no-go object", "closure_exported": False, "blocker": "excluded for every finite N by the transcendence contradiction"},
        {"object": "Lambda_origin_from_alpha_pi", "formula": "choose phase origin from alpha_geo/pi", "constructed": True, "accepted_as": "rejected source candidate", "closure_exported": False, "blocker": "irrational winding is dense/quasiperiodic, not a canonical origin/localizer"},
        {"object": "Omega_M_from_alpha_pi", "formula": "fix positive mass/action scale from alpha_geo and pi", "constructed": True, "accepted_as": "rejected unit candidate", "closure_exported": False, "blocker": "dimensionless irrationality does not break the positive scale torsor"},
    ]


def build_payload() -> dict[str, Any]:
    p3159 = load_json(INPUTS["P3159"])
    den_rows = denominator_scan()
    cf_rows = continued_fraction_ladder()
    best = min(den_rows, key=lambda row: row["absolute_phase_residual"])
    return {
        "status": "P3160_ALPHA_GEO_PHASE_LOCKING_NO_ROOT_OF_UNITY_THEOREM",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "content_grep": rg_hits(),
        "constructed_theoretical_objects": {
            "name": "R_alpha_phase_locking_no_root_of_unity_audit",
            "objects": object_matrix(),
            "exact_theorem_steps": theorem_steps(),
        },
        "numeric_context": {
            "alpha_geo": ALPHA_GEO,
            "two_pi": TWO_PI,
            "alpha_geo_over_2pi": RATIO,
            "A_phi": A_PHI,
            "p3159_status": p3159.get("status"),
            "best_denominator_scan_row_up_to_144": best,
        },
        "denominator_scan_rows_N_1_to_144": den_rows,
        "continued_fraction_ladder": cf_rows,
        "finite_theorem": {
            "name": "P3160_T1_alpha_geo_no_finite_phase_locking",
            "statement": "Because alpha_geo=ln(16), alpha_geo/(2*pi) cannot be rational: if alpha_geo/(2*pi)=p/q with p nonzero, then 16=e^(2*pi*p/q), but e^(pi*r) is transcendental for every nonzero rational r by the Gelfond-Schneider consequence, contradiction.  Hence exp(i*alpha_geo) is not a root of unity and alpha_geo does not define an exact finite Z_N phase slot for any N.  P3159's A_phi section remains valid as dimensionless normalization, but alpha_geo/pi does not source Lambda_origin, Omega_M/K_dim, a selector, bridge completion, role transfer, L_total, or ToE.",
            "denominator_rows": len(den_rows),
            "continued_fraction_rows": len(cf_rows),
            "constructed_objects": len(object_matrix()),
            "exact_locking_rows_excluded": len(den_rows),
            "accepted_closure_exports": 0,
        },
        "decision": {
            "bounded_result": "The stronger alpha_geo/pi check is now a real theorem, not just a float scan: alpha_geo is compatible with phase normalization through A_phi, but alpha_geo itself is provably not an exact rational phase slot/root of unity.  Therefore alpha_geo/pi cannot be used as a canonical Z_N phase-origin selector or unit-source law.",
            "next_honest_step": "Do not continue alpha_geo/pi phase-locking searches.  The next proof-grade move should construct exactly one new source object outside this exhausted route: either (1) a strict positive torsor source law Omega_scale selecting the scale of Omega_M/K_dim with a nonzero source value and scale-covariance-breaking proof, or (2) a nonconventional Lambda_origin_source_localizer whose origin is not obtained from alpha_geo/pi phase locking and comes with an explicit coupling theorem to Phi_Info/A_phi.",
            "negative_export_flags": {key: False for key in ["finite_phase_slot_exported", "root_of_unity_lock_exported", "Lambda_origin_exported", "Omega_M_scale_fixed", "K_dim_functor_exported", "physical_unit_source_exported", "QW_2191_discharged", "selector_closure_exported", "bridge_completion_exported", "role_transfer_exported", "unit_bearing_L_total_exported", "ToE_closure_exported"]},
            "positive_scoped_flags": {"exact_irrationality_theorem_exported": True, "no_finite_ZN_phase_locking_proved": True, "P3159_A_phi_section_preserved": True, "repo_grep_performed": True},
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["numeric_context"]
    th = payload["finite_theorem"]
    best = cert["best_denominator_scan_row_up_to_144"]
    lines = [
        "# P3160/S2110 alpha_geo phase-locking no-root-of-unity theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact theorem",
    ]
    for step in payload["constructed_theoretical_objects"]["exact_theorem_steps"]:
        lines.append(f"- `{step['step']}` `{step['status']}`: {step['claim']}")
    lines.extend([
        "",
        "## Numeric context",
        f"- `alpha_geo = {cert['alpha_geo']:.15f}`",
        f"- `alpha_geo/(2*pi) = {cert['alpha_geo_over_2pi']:.15f}`",
        f"- `A_phi = 2*pi/alpha_geo = {cert['A_phi']:.15f}`",
        f"- Best denominator row up to 144: `N={best['denominator_N']}`, `k={best['nearest_integer_k']}`, phase residual `{best['absolute_phase_residual']:.9f}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Decision",
        payload["decision"]["bounded_result"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3160/S2110 alpha_geo phase-locking no-root-of-unity theorem", "## P3160/S2110 alpha_geo phase-locking no-root-of-unity theorem\n\n`P3160/S2110` strengthens the P3159 alpha/pi audit by proving that `alpha_geo/(2*pi)` is irrational.  If `alpha_geo/(2*pi)=p/q`, then `16=e^(2*pi*p/q)`, contradicting the Gelfond-Schneider consequence that `e^(pi*r)` is transcendental for nonzero rational `r`.  Thus `exp(i*alpha_geo)` is not a root of unity and `alpha_geo` does not define an exact finite `Z_N` phase slot.  The P3159 section `A_phi=2*pi/alpha_geo` remains valid as dimensionless normalization, but no phase-origin source, `Omega_M/K_dim` scale, selector, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3160/S2110 alpha_geo irrational phase ratio is not a Lagrangian source", "## P3160/S2110 alpha_geo irrational phase ratio is not a Lagrangian source\n\n`P3160/S2110` proves the exact no-root-of-unity boundary for `alpha_geo` as a phase: it cannot be an exact finite `Z_N` phase lock.  This is a rigorous phase arithmetic result, but it supplies no dimensionful action unit, stress-energy source, selector, bridge/role-transfer, or `L_total`/ToE datum.\n")
    append_once(AGENTS, "Current alpha_geo phase-locking no-root-of-unity guardrail (P3160/S2110, 2026-07-13)", "## Current alpha_geo phase-locking no-root-of-unity guardrail (P3160/S2110, 2026-07-13)\n\n- P3160 strengthens P3159 by proving that `alpha_geo/(2*pi)` is irrational: otherwise `16=e^(2*pi*p/q)` would contradict the Gelfond-Schneider consequence that `e^(pi*r)` is transcendental for nonzero rational `r`.\n- Therefore `exp(i*alpha_geo)` is not a root of unity and `alpha_geo` cannot be an exact finite `Z_N` phase slot or phase-origin selector; the P3159 `A_phi=2*pi/alpha_geo` section remains only a dimensionless phase-area normalization.\n- Do not continue alpha_geo/pi phase-locking or `Z_N` slot searches as closure evidence; they do not export `Lambda_origin`, `Omega_M/K_dim`, physical units, selector closure, bridge completion, role transfer, `L_total`, or ToE.\n- Next honest move: introduce exactly one genuinely new source object outside alpha_geo/pi locking, preferably a strict positive torsor source law `Omega_scale` for `Omega_M/K_dim` or a nonconventional `Lambda_origin_source_localizer` with explicit `Phi_Info/A_phi` coupling.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
