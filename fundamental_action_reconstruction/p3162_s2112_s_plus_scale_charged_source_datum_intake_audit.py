"""P3162/S2112: S_+ scale-charged source datum intake audit.

P3161 proved that invariant strict data cannot select the free positive
Omega_M/K_dim torsor.  The exact missing object is therefore not another number,
but a scale-charged strict datum.  P3162 constructs this missing object as a
representation-theoretic source schema:

    S_+ in V_chi, with chi(c)=c under R_{>0} rescaling.

Acceptance requires a nonzero exported strict value of weight +1, an explicit
coupling S_+ -> Omega_M/K_dim, and no Planck/apparatus/selector/bridge/L_total
imports.  The audit then inventories current artifacts and candidate classes for
such a weight-one datum.
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
OUT = GEN / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json"
MD = GEN / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P3161": GEN / "p3161_s2111_omega_scale_positive_torsor_source_law_audit.json",
    "P3158": GEN / "p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json",
    "P3117": GEN / "p3117_s2067_omega_dim_dimension_character_source_audit.json",
    "P2957": GEN / "p2957_s1907_positive_beta_scale_unit_source_obstruction.json",
}
SCALE_FACTORS = [1 / 5, 1 / 2, 1, 2, 5]
GATES = (
    "typed_S_plus_object_constructed",
    "weight_plus_one_under_Rpos",
    "nonzero_value_exported",
    "strict_nadsoliton_source_exported",
    "couples_to_Omega_M_K_dim",
    "breaks_scale_orbit_nonconventionally",
    "not_dimensionless_invariant",
    "not_external_unit_import",
    "selector_bridge_ltotal_toe_free",
)


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
        "s_plus_shape": r"S_\+|scale-charged|scale charged|weight \+1|weight-one|R_\{>0\}|positive torsor",
        "unit_source_frontier": r"Omega_scale|Omega_M|K_dim|Omega_dim|positive_torsor_source_law|unit source|mass unit",
        "prior_no_go_inputs": r"P3161|P3158|P3117|P2957|scale orbit|scale-covariance|rescaling",
        "forbidden_lanes": r"Planck|apparatus|observed light|selector closure|bridge/role-transfer|role transfer|L_total|ToE",
    }
    out: dict[str, Any] = {}
    for name, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-g", "*.py", "-g", "*.md", "-g", "*.json"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        lines = [line for line in proc.stdout.splitlines() if line]
        out[name] = {"count": len(lines), "samples": lines[:20]}
    return out


def candidate_inventory() -> list[dict[str, Any]]:
    def cand(name: str, formula: str, weight: int | None, blocker: str, **flags: bool) -> dict[str, Any]:
        row = {gate: False for gate in GATES}
        row.update({
            "typed_S_plus_object_constructed": True,
            "weight_plus_one_under_Rpos": weight == 1,
            "not_dimensionless_invariant": weight not in (0, None),
            "not_external_unit_import": True,
            "selector_bridge_ltotal_toe_free": True,
        })
        row.update(flags)
        row.update({"candidate": name, "formula": formula, "scale_weight": weight, "blocker": blocker})
        return row

    return [
        cand("alpha_geo_phase_area", "S_+ := alpha_geo or A_phi", 0, "P3159/P3160: dimensionless phase normalization, not weight +1", nonzero_value_exported=True),
        cand("entropy_cardinality", "S_+ := exp(alpha_geo)=16", 0, "dimensionless count/cardinality", nonzero_value_exported=True),
        cand("z12_spectrum", "S_+ := Laplacian eigenvalue/gap", 0, "finite spectrum is dimensionless until metric/unit supplied", nonzero_value_exported=True),
        cand("damping_beta_scale", "S_+ := beta scale", 0, "P2957/P2692: positive beta remains orbit-gauge/target-dependent", nonzero_value_exported=True, couples_to_Omega_M_K_dim=True),
        cand("formal_Omega_M_symbol", "S_+ := Omega_M", 1, "formal carrier has correct weight but no strict nonzero source value", couples_to_Omega_M_K_dim=True, breaks_scale_orbit_nonconventionally=True),
        cand("formal_higgs_mu", "S_+ := sqrt(mu2)", 1, "mu2 is the unsourced Higgs parameter, so this is circular", couples_to_Omega_M_K_dim=True),
        cand("gamma_9_5_action_unit", "S_+ := Gamma_9_5", 1, "Gamma action-unit source remains unsourced in P2918 lane", couples_to_Omega_M_K_dim=True),
        cand("u_readout_placeholder", "S_+ := U_readout", 1, "readout unit is placeholder/gauge representative, not strict source", couples_to_Omega_M_K_dim=True),
        cand("planck_mass_import", "S_+ := m_Pl", 1, "external Planck import", nonzero_value_exported=True, strict_nadsoliton_source_exported=False, not_external_unit_import=False, couples_to_Omega_M_K_dim=True, breaks_scale_orbit_nonconventionally=True),
        cand("apparatus_mass_calibration", "S_+ := detector mass calibration", 1, "apparatus import", nonzero_value_exported=True, strict_nadsoliton_source_exported=False, not_external_unit_import=False, couples_to_Omega_M_K_dim=True, breaks_scale_orbit_nonconventionally=True),
        cand("selector_fixed_unit", "S_+ := selector-chosen positive representative", 1, "selector replay is blocked", nonzero_value_exported=True, selector_bridge_ltotal_toe_free=False, couples_to_Omega_M_K_dim=True),
        cand("new_strict_scale_charged_datum_schema", "S_+ := exported nonzero value in V_chi, chi(c)=c", 1, "this is the admissible schema but no current artifact exports the value", couples_to_Omega_M_K_dim=True, breaks_scale_orbit_nonconventionally=True),
    ]


def representation_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        weight = c["scale_weight"]
        for factor in SCALE_FACTORS:
            predicted = None if weight is None else factor ** weight
            rows.append({
                "candidate": c["candidate"],
                "scale_factor": factor,
                "scale_weight": weight,
                "predicted_transform_factor": predicted,
                "matches_required_weight_one": bool(weight == 1 and abs(predicted - factor) < 1e-12),
                "has_nonzero_exported_value": bool(c["nonzero_value_exported"]),
                "accepted_weight_row": bool(weight == 1 and c["nonzero_value_exported"] and c["strict_nadsoliton_source_exported"] and c["couples_to_Omega_M_K_dim"] and c["not_external_unit_import"] and c["selector_bridge_ltotal_toe_free"]),
                "blocker": c["blocker"],
            })
    return rows


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "gate": gate, "passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for candidate in sorted({g["candidate"] for g in gates}):
        sub = [g for g in gates if g["candidate"] == candidate]
        rows.append({"candidate": candidate, "passed_gates": sum(g["passed"] for g in sub), "failed_gates": sum(not g["passed"] for g in sub), "accepted_S_plus_source": all(g["passed"] for g in sub)})
    return rows


def theorem_steps() -> list[dict[str, str]]:
    return [
        {"step": "D1", "claim": "Define S_+ as a nonzero strict datum in a one-dimensional R_{>0} representation V_chi with chi(c)=c.", "status": "object_constructed"},
        {"step": "D2", "claim": "A value of weight 0 is dimensionless/invariant and is blocked by P3161's equivariant torsor obstruction.", "status": "reuses_P3161"},
        {"step": "D3", "claim": "A formal value of weight +1 has the right representation type but is not accepted unless its nonzero value is strictly exported and coupled to Omega_M/K_dim.", "status": "acceptance_gate"},
        {"step": "D4", "claim": "Planck/apparatus/selector representatives can have weight +1 but fail import-freedom or selector guardrails.", "status": "closed_lane_filter"},
        {"step": "D5", "claim": "Current artifacts export no nonzero strict source value in V_chi; therefore S_+ is specified as the next required object, not closed.", "status": "bounded_no_go"},
    ]


def build_payload() -> dict[str, Any]:
    candidates = candidate_inventory()
    rep_rows = representation_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [a for a in aggs if a["accepted_S_plus_source"]]
    return {
        "status": "P3162_S_PLUS_SCALE_CHARGED_SOURCE_DATUM_INTAKE_BOUNDED_NO_GO",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "input_statuses": {name: load_json(path).get("status") for name, path in INPUTS.items()},
        "content_grep": rg_hits(),
        "constructed_theoretical_objects": {
            "object": "S_plus_scale_charged_source_datum",
            "definition": "S_+ in V_chi, chi(c)=c for c in R_{>0}",
            "acceptance_theorem_steps": theorem_steps(),
            "candidate_S_plus_sources": candidates,
            "representation_rows": rep_rows,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "candidate_S_plus_sources": len(candidates),
            "scale_factors": len(SCALE_FACTORS),
            "representation_rows": len(rep_rows),
            "candidate_gate_rows": len(gates),
            "weight_plus_one_candidates": sum(c["scale_weight"] == 1 for c in candidates),
            "weight_zero_candidates": sum(c["scale_weight"] == 0 for c in candidates),
            "accepted_S_plus_sources": len(accepted),
            "accepted_weight_rows": sum(r["accepted_weight_row"] for r in rep_rows),
        },
        "finite_theorem": {
            "name": "P3162_T1_S_plus_acceptance_schema_no_current_export",
            "statement": "The missing scale object can be stated precisely as S_+ in the weight-one representation of R_{>0}.  This representation type is necessary to evade the P3161 invariant-data torsor obstruction, but current artifacts only provide dimensionless invariants, formal placeholders, circular Higgs/unit symbols, or external/selector imports.  No nonzero strict S_+ value with an Omega_M/K_dim coupling is exported.",
        },
        "decision": {
            "bounded_result": "P3162 upgrades the next frontier from a vague 'unit source' to a typed object: S_+ must be a nonzero strict weight-one datum.  The finite intake finds zero accepted current exports; formal Omega_M, Gamma_9_5, U_readout, and Higgs-mu symbols have the right weight shape but no strict source value, while alpha/pi, entropy, and spectra are weight-zero invariants.",
            "next_honest_step": "The next honest step is now exactly one of two options: (1) actually export a nonzero strict value of S_+ with a theorem coupling it to Omega_M/K_dim, then rerun P3162/P3161; or (2) pivot to Lambda_origin_source_localizer if the goal is phase-origin rather than unit scale.  Without one of these supplied objects, preserve the no-strict-unit/no-new-live-frontier certificate.",
            "negative_export_flags": {key: False for key in ["S_plus_source_exported", "Omega_scale_source_exported", "Omega_M_scale_fixed", "K_dim_functor_exported", "physical_unit_source_exported", "Higgs_VEV_source_exported", "EH_coupling_exported", "QW_2191_discharged", "selector_closure_exported", "bridge_completion_exported", "role_transfer_exported", "unit_bearing_L_total_exported", "ToE_closure_exported"]},
            "positive_scoped_flags": {"S_plus_object_constructed": True, "weight_one_acceptance_schema_defined": True, "candidate_intake_matrix_built": True, "repo_grep_performed": True},
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3162/S2112 S_+ scale-charged source datum intake audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_theoretical_objects']['object']}`", f"- Definition: `{payload['constructed_theoretical_objects']['definition']}`", "", "## Acceptance theorem"]
    for step in payload["constructed_theoretical_objects"]["acceptance_theorem_steps"]:
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
    append_once(STRICT_EQUATION_SHEET, "P3162/S2112 S_plus scale-charged source datum intake audit", "## P3162/S2112 S_plus scale-charged source datum intake audit\n\n`P3162/S2112` constructs the missing `S_+` object named by P3161: a nonzero strict datum in the weight-one representation `V_chi` of `R_{>0}`, with `chi(c)=c`, intended to couple to `Omega_M/K_dim`.  The audit inventories `12` candidate classes over `5` scale factors, builds `60` representation rows and `108` gate rows, and finds `0` accepted `S_+` sources.  Weight-zero alpha/pi, entropy, and spectral data remain blocked by P3161; formal `Omega_M`, `Gamma_9_5`, `U_readout`, and Higgs symbols have the right weight shape but no strict source value; Planck/apparatus/selector representatives are forbidden imports.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3162/S2112 S_plus weight-one source remains unexported", "## P3162/S2112 S_plus weight-one source remains unexported\n\n`P3162/S2112` states the missing unit-source datum precisely as `S_+ in V_chi`, `chi(c)=c`.  Current artifacts do not export a nonzero strict `S_+` value or coupling theorem to `Omega_M/K_dim`, so no action/mass unit, Higgs VEV, EH coupling, `L_total`, or ToE datum is promoted.\n")
    append_once(AGENTS, "Current S_plus scale-charged source datum guardrail (P3162/S2112, 2026-07-13)", "## Current S_plus scale-charged source datum guardrail (P3162/S2112, 2026-07-13)\n\n- P3162 constructs the missing `S_+` object required by P3161: a nonzero strict datum in the weight-one `R_{>0}` representation `V_chi`, `chi(c)=c`, coupled to `Omega_M/K_dim`.\n- The intake audit tests `12` candidate source classes over `5` scale factors, with `60` representation rows and `108` gate rows; `0` candidates export an accepted strict `S_+` source.\n- Weight-zero alpha/pi, entropy/cardinality, and Z12 spectra remain dimensionless; formal `Omega_M`, `Gamma_9_5`, `U_readout`, and Higgs symbols lack strict source values; Planck/apparatus/selector representatives are forbidden.\n- Do not promote any current candidate to physical units, Higgs VEV, EH coupling, selector closure, bridge completion, role transfer, `L_total`, or ToE.\n- Next honest move: either actually export a nonzero strict `S_+` value with an `Omega_M/K_dim` coupling theorem, or pivot to `Lambda_origin_source_localizer`; otherwise preserve the no-strict-unit/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
