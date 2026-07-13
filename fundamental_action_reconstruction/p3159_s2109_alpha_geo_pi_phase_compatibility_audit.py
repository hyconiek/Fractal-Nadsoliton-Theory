"""P3159/S2109: alpha_geo/pi phase compatibility audit.

The user asked to check the compatibility of strict alpha_geo with the phase-pi
lane and to take the next honest proof/computation step after P3158.  This
packet constructs one bounded object, Phi_alpha_pi, testing whether
alpha_geo=4 ln 2 can canonically couple to a 2*pi phase quantum without
silently promoting it to a physical unit, selector, bridge, or role transfer.
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
OUT = GEN / "p3159_s2109_alpha_geo_pi_phase_compatibility_audit.json"
MD = GEN / "p3159_s2109_alpha_geo_pi_phase_compatibility_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "alpha_geo": GEN / "alpha_geo_strict_derived_v1.json",
    "P3111": GEN / "p3111_s2061_symplectic_phase_area_section_audit.json",
    "P3112": GEN / "p3112_s2062_c_phi_calibration_functional_audit.json",
    "P3158": GEN / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json",
}
ALPHA_GEO = 4.0 * math.log(2.0)
TWO_PI = 2.0 * math.pi
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
        "alpha_geo_pi_phase": r"alpha_geo|A_phi|2\*pi/alpha_geo|2\*pi|phase-area|pi phase|π",
        "unit_chain": r"P3111|P3112|P3116|P3158|Omega_dim|K_dim|C_phi|U_action",
        "forbidden_promotion": r"role transfer|bridge/role-transfer|selector closure|QW-2191|L_total|ToE|Planck|apparatus",
    }
    out = {}
    for name, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-g", "*.py", "-g", "*.md", "-g", "*.json"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        lines = [line for line in proc.stdout.splitlines() if line]
        out[name] = {"count": len(lines), "samples": lines[:18]}
    return out


def phase_rows() -> list[dict[str, Any]]:
    rows = []
    for n in range(1, 13):
        phase = (n * ALPHA_GEO) % TWO_PI
        nearest_slot = round(12.0 * phase / TWO_PI) % 12
        slot_phase = nearest_slot * TWO_PI / 12.0
        signed_error = ((phase - slot_phase + math.pi) % TWO_PI) - math.pi
        rows.append({
            "n": n,
            "n_alpha_geo_mod_2pi": phase,
            "nearest_Z12_phase_slot": nearest_slot,
            "signed_error_to_nearest_Z12_slot": signed_error,
            "absolute_error": abs(signed_error),
            "exact_2pi_multiple_by_A_phi": n,
            "A_phi_section_value": A_PHI,
        })
    return rows


def continued_fraction_rows() -> list[dict[str, Any]]:
    x = ALPHA_GEO / TWO_PI
    rows = []
    for max_den in [12, 24, 48, 96, 192, 384, 768, 1536]:
        q = Fraction(x).limit_denominator(max_den)
        rows.append({"max_denominator": max_den, "p": q.numerator, "q": q.denominator, "p_over_q": q.numerator / q.denominator, "error_abs": abs(x - q.numerator / q.denominator), "phase_error_2pi_units": abs(TWO_PI * (x - q.numerator / q.denominator))})
    return rows


def compatibility_matrix() -> list[dict[str, Any]]:
    specs = [
        ("dimensionless_phase_area_section", "alpha_geo*A_phi=2*pi", True, True, True, True, False, False, "This is exactly the P3111-style compatibility: alpha_geo can normalize a dimensionless 2*pi phase-area section."),
        ("Z12_slot_alignment", "n*alpha_geo mod 2*pi lands on a canonical Z12 phase slot", True, False, True, True, False, False, "Finite n<=12 scan has no exact slot-locking theorem; numeric nearest slots are approximation data only."),
        ("physical_action_calibration", "C_phi(A_phi)=U_action", True, True, False, True, False, False, "P3112/P3113 require a dimensionful reference carrier; pi-phase compatibility alone is dimensionless."),
        ("phase_origin_source_localizer", "alpha_geo/pi selects Lambda_origin", True, False, True, False, False, False, "alpha_geo and pi are even scalar constants; no nonconventional phase-origin/localizer sign is selected."),
        ("mass_unit_positive_torsor_source", "alpha_geo and 2*pi fix Omega_M scale", True, True, False, True, False, False, "The formal ratio fixes a dimensionless number, not the positive mass/action scale torsor."),
        ("legacy_role_transfer", "alpha_geo/pi phase compatibility transfers legacy EW/EM roles", False, False, False, False, False, False, "Forbidden without explicit bridge completion and role-transfer theorem."),
    ]
    keys = ["candidate", "claim", "alpha_geo_used_strictly", "pi_phase_compatible", "unit_source_exported", "selector_or_origin_exported", "accepted_strict_closure", "role_transfer_allowed", "verdict"]
    return [dict(zip(keys, s)) for s in specs]


def build_payload() -> dict[str, Any]:
    p3111 = load_json(INPUTS["P3111"])
    rows = phase_rows()
    matrix = compatibility_matrix()
    exact_section_residual = ALPHA_GEO * A_PHI - TWO_PI
    min_slot = min(rows, key=lambda r: r["absolute_error"])
    accepted = [r for r in matrix if r["accepted_strict_closure"]]
    return {
        "status": "P3159_ALPHA_GEO_PI_PHASE_COMPATIBILITY_SCOPED_SECTION_NO_UNIT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "content_grep": rg_hits(),
        "constructed_object": {
            "name": "Phi_alpha_pi",
            "formula": "A_phi := 2*pi/alpha_geo = pi/(2 ln 2), so alpha_geo*A_phi=2*pi",
            "classification": "dimensionless phase-area compatibility section, not a physical unit/source theorem",
        },
        "numeric_certificate": {
            "alpha_geo": ALPHA_GEO,
            "two_pi": TWO_PI,
            "A_phi": A_PHI,
            "alpha_geo_times_A_phi_minus_2pi": exact_section_residual,
            "alpha_geo_over_2pi": ALPHA_GEO / TWO_PI,
            "p3111_status": p3111.get("status"),
            "min_Z12_slot_error_n_le_12": min_slot,
        },
        "phase_rows_n_1_to_12": rows,
        "continued_fraction_rows_alpha_over_2pi": continued_fraction_rows(),
        "compatibility_matrix": matrix,
        "finite_theorem": {
            "name": "P3159_T1_alpha_geo_pi_phase_scoped_compatibility",
            "statement": "The strict scalar alpha_geo=4 ln 2 is exactly compatible with a dimensionless phase-area section A_phi=2*pi/alpha_geo.  This verifies the pi-phase bookkeeping route and agrees with the P3111 style phase-area section.  However, finite Z12 slot scans and rational approximants do not provide a canonical phase-origin, selector, physical action/mass unit, or legacy role-transfer theorem.  The result is positive as internal phase normalization and negative as closure/source promotion.",
            "accepted_strict_closures": len(accepted),
            "phase_rows": len(rows),
            "continued_fraction_rows": 8,
            "compatibility_rows": len(matrix),
        },
        "decision": {
            "bounded_result": "alpha_geo is compatible with pi-phase normalization only as the dimensionless A_phi=2*pi/alpha_geo section.  This is useful bookkeeping for phase-area audits, but it does not construct Lambda_origin, Omega_M/K_dim, a selector, unit-bearing action, bridge completion, or role transfer.",
            "next_honest_step": "Recommended next step: construct exactly one new source object rather than replaying alpha_geo/pi numerology.  The least-replay option is a strict positive torsor source law for Omega_M/K_dim; the alternative is a nonconventional Lambda_origin_source_localizer coupled to the already-compatible A_phi section.  In either case, acceptance must require a nonzero source value and an explicit coupling theorem, not just alpha_geo=4 ln 2 or A_phi=2*pi/alpha_geo.",
            "negative_export_flags": {key: False for key in ["physical_unit_source_exported", "Omega_M_scale_fixed", "K_dim_functor_exported", "Lambda_origin_exported", "phase_origin_selected", "QW_2191_discharged", "selector_closure_exported", "bridge_completion_exported", "role_transfer_exported", "unit_bearing_L_total_exported", "ToE_closure_exported"]},
            "positive_scoped_flags": {"alpha_geo_pi_phase_compatibility_verified": True, "A_phi_dimensionless_section_constructed": True, "repo_grep_performed": True, "next_step_narrowed_to_one_source_object": True},
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["numeric_certificate"]
    th = payload["finite_theorem"]
    lines = ["# P3159/S2109 alpha_geo/pi phase compatibility audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['name']}`", f"- Formula: `{payload['constructed_object']['formula']}`", f"- Classification: {payload['constructed_object']['classification']}", "", "## Numeric certificate", f"- `alpha_geo = {cert['alpha_geo']:.15f}`", f"- `2*pi = {cert['two_pi']:.15f}`", f"- `A_phi = 2*pi/alpha_geo = {cert['A_phi']:.15f}`", f"- `alpha_geo*A_phi - 2*pi = {cert['alpha_geo_times_A_phi_minus_2pi']:.3e}`", f"- Best n<=12 Z12 slot row: `n={cert['min_Z12_slot_error_n_le_12']['n']}`, error `{cert['min_Z12_slot_error_n_le_12']['absolute_error']:.6f}`", "", "## Finite theorem", f"`{th['name']}`: {th['statement']}", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3159/S2109 alpha_geo/pi phase compatibility audit", "## P3159/S2109 alpha_geo/pi phase compatibility audit\n\n`P3159/S2109` constructs `Phi_alpha_pi`, a scoped compatibility object for `alpha_geo=4 ln 2` and the `2*pi` phase quantum.  The exact dimensionless section `A_phi=2*pi/alpha_geo=pi/(2 ln 2)` satisfies `alpha_geo*A_phi=2*pi`, confirming the pi-phase bookkeeping route already suggested by the P3111 phase-area section.  Finite `n<=12` phase-slot scans and rational approximants do not export a phase-origin source, selector, physical action/mass unit, bridge completion, role transfer, or `L_total`/ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3159/S2109 alpha_geo/pi phase section is not a unit source", "## P3159/S2109 alpha_geo/pi phase section is not a unit source\n\n`P3159/S2109` verifies that `alpha_geo` is compatible with `2*pi` as a dimensionless phase-area normalization, but it remains below a Lagrangian/action-unit source.  No `Omega_M`, `K_dim`, `Lambda_origin`, EH/SM coupling, selector, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(AGENTS, "Current alpha_geo/pi phase compatibility guardrail (P3159/S2109, 2026-07-13)", "## Current alpha_geo/pi phase compatibility guardrail (P3159/S2109, 2026-07-13)\n\n- P3159 constructs `Phi_alpha_pi`, a scoped compatibility object between strict `alpha_geo=4 ln 2` and the `2*pi` phase quantum.\n- The audit verifies the exact dimensionless phase-area section `A_phi=2*pi/alpha_geo=pi/(2 ln 2)`, so alpha_geo is compatible with pi-phase bookkeeping at the internal normalization level.\n- Finite `Z12` slot scans and rational approximants do not export a canonical phase origin, `Lambda_origin`, positive mass/action unit, `K_dim`, selector closure, bridge completion, role transfer, `L_total`, or ToE.\n- Next honest move: introduce exactly one new source object, preferably a strict positive torsor source law for `Omega_M/K_dim` or a nonconventional `Lambda_origin_source_localizer` explicitly coupled to `A_phi`; do not replay alpha_geo/pi numerology as closure.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
