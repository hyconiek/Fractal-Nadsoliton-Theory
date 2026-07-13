"""P3163/S2113: boundary-value speed-of-light matching audit.

The user asked whether dimension search should proceed by first identifying
model boundary/limiting values and trying to match them to physical values such
as the speed of light.  This audit constructs exactly that receiver, while
keeping the unit-source guardrails from P3116-P3162.

Core theorem: any nonzero dimensionless model boundary value v_hat can be made
to match the SI speed of light c by choosing U_L/U_T = c/v_hat.  Therefore a
boundary-value fit is a useful receiver/diagnostic, but it is not a source of
physical dimensions unless U_L and U_T are independently strict-sourced.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3163_s2113_boundary_value_speed_of_light_matching_audit.json"
MD = GEN / "p3163_s2113_boundary_value_speed_of_light_matching_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P3162": GEN / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json",
    "P3158": GEN / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json",
    "P3119": GEN / "p3119_s2069_xi_lt_axis_source_object_audit.json",
    "P3120": GEN / "p3120_s2070_tau_lt_ordered_flow_source_audit.json",
    "P3077": GEN / "p3077_s2027_second_order_lift_obstruction_table.json",
}
C_SI = 299_792_458.0
ALPHA_GEO = 4.0 * math.log(2.0)
A_PHI = 2.0 * math.pi / ALPHA_GEO
OMEGA_STRICT = 0.18575
PHI_STRICT = 0.16250
BETA_STRICT = 1.0
ETA_STRICT = 1.8


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
        "boundary_light_lane": r"speed of light|observed light|light speed|c-retardation|Lorentz|spacetime embedding|unit-normalized.*metric",
        "dimension_chain": r"U_length|U_time|Xi_LT|Tau_LT|R_dim|K_dim|Omega_dim|S_\+|Omega_scale",
        "boundary_values": r"boundary value|limiting value|limit|asymptotic|Dirichlet|Laplacian|A_phi|alpha_geo|K_strict_gate",
        "closed_imports": r"Planck|apparatus|observed light|selector closure|bridge/role-transfer|role transfer|L_total|ToE",
    }
    out: dict[str, Any] = {}
    for name, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-g", "*.py", "-g", "*.md", "-g", "*.json"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        lines = [line for line in proc.stdout.splitlines() if line]
        out[name] = {"count": len(lines), "samples": lines[:20]}
    return out


def strict_kernel(d: float) -> float:
    return math.cos(OMEGA_STRICT * d + PHI_STRICT) / (1.0 + BETA_STRICT * (d ** ETA_STRICT))


def boundary_candidates() -> list[dict[str, Any]]:
    lap_gap = 2.0 - 2.0 * math.cos(2.0 * math.pi / 12.0)
    lap_max = 4.0
    rows = [
        ("unit_receiver", "normalized causal receiver speed", 1.0, "chosen receiver normalization"),
        ("alpha_over_2pi", "alpha_geo/(2*pi)", ALPHA_GEO / (2.0 * math.pi), "P3160 irrational dimensionless winding ratio"),
        ("A_phi", "2*pi/alpha_geo", A_PHI, "P3159/P3111 phase-area section"),
        ("strict_kernel_d0", "K_strict_gate(0)=cos(phi)", strict_kernel(0.0), "strict kernel zero-distance boundary"),
        ("strict_kernel_d1", "K_strict_gate(1)", strict_kernel(1.0), "strict kernel first positive sample"),
        ("strict_kernel_tail_limit", "lim_{d->infty} K_strict_gate(d)", 0.0, "tail boundary is zero and cannot fit nonzero c"),
        ("z12_laplacian_gap", "2-2*cos(2*pi/12)", lap_gap, "dimensionless first Laplacian gap"),
        ("z12_laplacian_max", "max cyclic Laplacian eigenvalue", lap_max, "dimensionless spectral boundary"),
        ("dirichlet_formal_wave_speed", "formal lattice speed from P3077 second-order lift", 1.0, "formal imported wave-compatible normalization, not internal physical c"),
        ("entropy_cell_count", "exp(alpha_geo)=16", math.exp(ALPHA_GEO), "dimensionless cardinality"),
    ]
    return [{"candidate": name, "formula": formula, "v_hat": value, "source_context": context, "nonzero": value != 0.0} for name, formula, value, context in rows]


def fit_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        v = c["v_hat"]
        if v == 0.0:
            rows.append({"candidate": c["candidate"], "v_hat": v, "c_SI_m_per_s": C_SI, "required_U_length_over_U_time": None, "can_numerically_fit_c": False, "fit_residual_after_scale_choice": None, "source_closure": False, "blocker": "zero boundary value cannot fit nonzero c"})
        else:
            ratio = C_SI / v
            rows.append({"candidate": c["candidate"], "v_hat": v, "c_SI_m_per_s": C_SI, "required_U_length_over_U_time": ratio, "can_numerically_fit_c": True, "fit_residual_after_scale_choice": abs(v * ratio - C_SI), "source_closure": False, "blocker": "fit chooses U_length/U_time; no independent strict U_length and U_time sources"})
    return rows


def scale_degeneracy_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    factors = [0.1, 0.5, 1.0, 2.0, 10.0]
    rows = []
    for c in candidates:
        if c["v_hat"] == 0.0:
            continue
        base_ratio = C_SI / c["v_hat"]
        for a in factors:
            for b in factors:
                # U_L -> a U_L and U_T -> b U_T changes the physical velocity by a/b.
                fit_ratio = base_ratio * b / a
                physical_speed = c["v_hat"] * fit_ratio * a / b
                rows.append({"candidate": c["candidate"], "length_rescale_a": a, "time_rescale_b": b, "adjusted_U_length_over_U_time": fit_ratio, "physical_speed_after_adjustment": physical_speed, "residual_to_c": abs(physical_speed - C_SI), "degenerate_fit_family": True})
    return rows


def gate_rows(candidates: list[dict[str, Any]], fits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    gates = []
    fit_by = {r["candidate"]: r for r in fits}
    for c in candidates:
        f = fit_by[c["candidate"]]
        checks = {
            "boundary_value_constructed": True,
            "nonzero_receiver_value": c["nonzero"],
            "numerically_matches_c_after_scale_choice": f["can_numerically_fit_c"],
            "independent_U_length_source_exported": False,
            "independent_U_time_source_exported": False,
            "Lorentz_spacetime_embedding_exported": False,
            "observed_light_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
        }
        for gate, passed in checks.items():
            gates.append({"candidate": c["candidate"], "gate": gate, "passed": bool(passed), "detail": "passed" if passed else f["blocker"]})
    return gates


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for candidate in sorted({g["candidate"] for g in gates}):
        sub = [g for g in gates if g["candidate"] == candidate]
        out.append({"candidate": candidate, "passed_gates": sum(g["passed"] for g in sub), "failed_gates": sum(not g["passed"] for g in sub), "accepted_physical_c_source": all(g["passed"] for g in sub)})
    return out


def theorem_steps() -> list[dict[str, str]]:
    return [
        {"step": "C1", "claim": "Let v_hat be any nonzero dimensionless boundary velocity receiver from the strict model.", "status": "receiver_definition"},
        {"step": "C2", "claim": "A physical velocity has dimension U_length/U_time, so v_phys = v_hat*(U_length/U_time).", "status": "dimension_lift"},
        {"step": "C3", "claim": "For any target c>0, choosing U_length/U_time=c/v_hat gives v_phys=c.", "status": "fit_degeneracy"},
        {"step": "C4", "claim": "Therefore numerical agreement with c is not a source theorem unless U_length and U_time are independently strict-sourced and a Lorentz/light embedding is exported.", "status": "no_closure"},
    ]


def build_payload() -> dict[str, Any]:
    candidates = boundary_candidates()
    fits = fit_rows(candidates)
    degeneracy = scale_degeneracy_rows(candidates)
    gates = gate_rows(candidates, fits)
    aggs = aggregate(gates)
    accepted = [a for a in aggs if a["accepted_physical_c_source"]]
    return {
        "status": "P3163_BOUNDARY_VALUE_SPEED_OF_LIGHT_MATCHING_UNDERDETERMINATION_AUDIT",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "input_statuses": {name: load_json(path).get("status") for name, path in INPUTS.items()},
        "content_grep": rg_hits(),
        "constructed_theoretical_objects": {
            "object": "BoundaryValuePhysicalConstantReceiverAudit",
            "target_physical_constant": "speed_of_light_c_SI_exact_299792458_m_per_s",
            "theorem_steps": theorem_steps(),
            "boundary_value_candidates": candidates,
            "speed_of_light_fit_rows": fits,
            "scale_degeneracy_rows": degeneracy,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "boundary_value_candidates": len(candidates),
            "speed_of_light_fit_rows": len(fits),
            "numerically_fit_nonzero_candidates": sum(r["can_numerically_fit_c"] for r in fits),
            "scale_degeneracy_rows": len(degeneracy),
            "candidate_gate_rows": len(gates),
            "accepted_physical_c_sources": len(accepted),
            "c_SI_m_per_s": C_SI,
        },
        "finite_theorem": {
            "name": "P3163_T1_boundary_c_fit_is_scale_underdetermined",
            "statement": "Boundary/limiting model values can be used as receiver diagnostics for physical constants, but any nonzero dimensionless v_hat can match c by choosing U_length/U_time=c/v_hat.  Since current artifacts do not export independent U_length, U_time, Lorentz/spacetime embedding, or observed-light source, no boundary-value fit to c is a strict dimension source.",
        },
        "decision": {
            "bounded_result": "P3163 constructs the requested boundary-value-to-c matching audit.  Nine nonzero boundary receivers can numerically fit c after choosing a unit ratio, and the zero tail cannot; all successful fits are scale choices rather than strict source laws.",
            "next_honest_step": "If this boundary-value route is pursued, the next proof-grade object must be a two-axis unit theorem U_LT: independently source U_length and U_time (or a Lorentzian metric/light-cone embedding) before comparing to c.  Otherwise preserve P3162's S_+ unit-source frontier or pivot to Lambda_origin_source_localizer for phase-origin work.",
            "negative_export_flags": {key: False for key in ["physical_speed_of_light_source_exported", "U_length_source_exported", "U_time_source_exported", "Lorentz_metric_exported", "observed_light_exported", "S_plus_source_exported", "Omega_scale_source_exported", "K_dim_functor_exported", "physical_unit_source_exported", "QW_2191_discharged", "selector_closure_exported", "bridge_completion_exported", "role_transfer_exported", "unit_bearing_L_total_exported", "ToE_closure_exported"]},
            "positive_scoped_flags": {"boundary_value_receiver_constructed": True, "speed_of_light_fit_matrix_built": True, "scale_degeneracy_theorem_proved": True, "repo_grep_performed": True},
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3163/S2113 boundary-value speed-of-light matching audit", "", f"Status: `{payload['status']}`", "", "## Theorem steps"]
    for step in payload["constructed_theoretical_objects"]["theorem_steps"]:
        lines.append(f"- `{step['step']}` `{step['status']}`: {step['claim']}")
    lines.extend(["", "## Finite certificate"])
    for k, v in cert.items():
        lines.append(f"- `{k}`: `{v}`")
    lines.extend(["", "## Finite theorem", f"`{payload['finite_theorem']['name']}`: {payload['finite_theorem']['statement']}", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3163/S2113 boundary-value speed-of-light matching audit", "## P3163/S2113 boundary-value speed-of-light matching audit\n\n`P3163/S2113` constructs a boundary-value receiver audit for fitting model limiting values to the physical speed of light `c=299792458 m/s`.  It tests `10` boundary receivers, including `alpha_geo/(2*pi)`, `A_phi`, strict-kernel boundary samples, Z12 Laplacian values, a formal wave-speed receiver, and `exp(alpha_geo)`.  The theorem is that any nonzero dimensionless `v_hat` fits `c` after choosing `U_length/U_time=c/v_hat`; hence the `9` successful numerical fits are scale choices, not unit-source theorems.  No independent `U_length`, `U_time`, Lorentz metric/light-cone embedding, observed-light source, `S_+`, `Omega_scale`, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3163/S2113 boundary c-fit remains receiver-only", "## P3163/S2113 boundary c-fit remains receiver-only\n\n`P3163/S2113` shows that boundary-value matching to `c` is useful as a diagnostic receiver but underdetermined as a dimension source: the fit always chooses `U_length/U_time` unless those units are independently sourced.  It does not promote a Lorentzian metric, observed-light sector, action unit, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current boundary-value speed-of-light matching guardrail (P3163/S2113, 2026-07-13)", "## Current boundary-value speed-of-light matching guardrail (P3163/S2113, 2026-07-13)\n\n- P3163 constructs a boundary-value receiver audit for comparing model limiting values with the physical speed of light `c=299792458 m/s`.\n- The audit tests `10` boundary receivers; `9` nonzero receivers can numerically match `c` after choosing `U_length/U_time=c/v_hat`, while the zero tail cannot.\n- This proves scale underdetermination, not physical closure: current artifacts do not independently export `U_length`, `U_time`, a Lorentzian metric/light-cone embedding, observed-light source, `S_+`, or `Omega_scale`.\n- Do not promote boundary-value fits to physical dimensions, observed light, spacetime EOM, selector closure, bridge completion, role transfer, `L_total`, or ToE.\n- Next honest move: construct a two-axis unit theorem `U_LT` sourcing `U_length` and `U_time` or a Lorentzian metric/light-cone embedding before fitting `c`; otherwise preserve the P3162 `S_+` frontier or pivot to `Lambda_origin_source_localizer`.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
