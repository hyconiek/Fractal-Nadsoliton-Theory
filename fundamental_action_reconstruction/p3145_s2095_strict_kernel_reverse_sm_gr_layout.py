"""P3145/S2095: reverse SM/GR layout over the strict kernel.

The user asked for an axiomatic reverse analysis: if the observed universe and
known SM/GR physics must fit inside the nadsoliton, where would those physical
properties have to sit in the strict kernel picture?  This packet constructs a
finite layout matrix from current repo evidence.  It is deliberately a reverse
placement / obligation map, not a closure theorem and not legacy role transfer.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.json"
MD = GEN / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3144": GEN / "p3144_s2094_upsilon_sel_unit_measure_obstruction.json",
    "P3143": GEN / "p3143_s2093_vlift_weight_scale_source_audit.json",
    "P3104": GEN / "p3104_s2054_spectral_triple_geometry_interface_obstruction_audit.json",
    "P3102": GEN / "p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit.json",
    "P3094": GEN / "p3094_s2044_stress_energy_metric_response_obstruction_audit.json",
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
}

GATES = [
    "strict_carrier_present",
    "receiver_scaffold_present",
    "source_law_present",
    "unit_or_calibration_present",
    "nonproxy_eom_or_dynamics_present",
    "selector_obstruction_resolved_when_needed",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_text_or_json(path: Path) -> Any:
    if not path.exists():
        return None
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    return {"text_sha": sha(path), "status": "markdown_input_present"}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def physical_rows() -> list[dict[str, Any]]:
    # The gate values are current-repo layout judgments, not physical claims.
    return [
        {
            "physical_property": "spacetime metric / geometric distance",
            "strict_slot": "finite Z12 geometry, Laplacian/spectral-triple-like rows, alpha_geo-scaled distance witnesses",
            "repo_backing": ["P3104", "K_strict_gate operational geometry controls"],
            "legacy_analogy": "legacy kernel treated geometry/amplitude ontologically; strict side only has finite geometry receivers so far",
            "gates": [True, True, False, False, False, True],
            "missing_object": "physical length/action unit and geometry/readout interface",
        },
        {
            "physical_property": "GR stress-energy and Einstein dynamics",
            "strict_slot": "stress-energy/metric-response obstruction rows plus Lagrangian/EOM draft lane",
            "repo_backing": ["P3094", "P2685-P2687/P3144"],
            "legacy_analogy": "legacy gravity hierarchy roles cannot transfer silently; strict side needs nonproxy metric response",
            "gates": [True, True, False, False, False, True],
            "missing_object": "nonproxy EH/ELg stress-energy coupling with unit-bearing action measure",
        },
        {
            "physical_property": "SM gauge scaffold SU(3)xSU(2)xU(1)",
            "strict_slot": "strict-core partial gauge reconstruction / finite mode and character scaffolds",
            "repo_backing": ["README A6/L18/L19", "P2991 Fourier character bookkeeping"],
            "legacy_analogy": "legacy formulas suggest gauge roles; strict side requires source and uniqueness, not inheritance",
            "gates": [True, True, False, False, False, False],
            "missing_object": "unique internal gauge-mode selection and selector-safe quotient theorem",
        },
        {
            "physical_property": "electroweak and electromagnetic coupling constants",
            "strict_slot": "strict successor constants / alpha_geo and beta/Z_beta source atoms",
            "repo_backing": ["P2691", "P2692", "legacy role-transfer guardrail"],
            "legacy_analogy": "legacy sin^2(theta_W) and alpha_EM formulas are known role-transfer hazards",
            "gates": [True, True, False, False, False, False],
            "missing_object": "role-safe amplitude/coupling source theorem plus bridge/role-transfer theorem",
        },
        {
            "physical_property": "fermion masses, Yukawa-like hierarchy, flavor",
            "strict_slot": "direct-route mass/Yukawa receiver equations and P46/P2695/P2696 residual matrices",
            "repo_backing": ["F3", "P2695", "P2696"],
            "legacy_analogy": "legacy hierarchy intuition can guide targets but cannot fill zero-witness gaps",
            "gates": [True, True, False, False, False, False],
            "missing_object": "strict carrier/source for remaining pair equations and target action coefficients",
        },
        {
            "physical_property": "chirality, orientation, matter/antimatter sign, selector",
            "strict_slot": "Z12 orientation torsors, chiral-bispectrum marker, D_HL/axiom branch, Upsilon_sel obstruction",
            "repo_backing": ["P2708-P2721", "P3133-P3144"],
            "legacy_analogy": "legacy phase/torsion split can create conditional D_HL shapes but not strict origin/polarity source",
            "gates": [True, True, False, False, False, False],
            "missing_object": "non-premise symmetry-breaking source fixing orientation/polarity without A_origin/A_lambda",
        },
        {
            "physical_property": "RG flow, damping/compression, hierarchy scaling",
            "strict_slot": "K_strict_gate beta/eta nonlinear damping and micro/Stage-C operational gate",
            "repo_backing": ["K2", "P2692", "P2377/P2378 lanes"],
            "legacy_analogy": "legacy beta_tors linear damping is an intermediate bridge clue, not a strict substitute",
            "gates": [True, True, False, False, False, True],
            "missing_object": "target-independent positive beta/Z_beta source and unit-normalized transport law",
        },
        {
            "physical_property": "quantum probabilities, Born/readout, detector calibration",
            "strict_slot": "formal probability/readout rows on Z12 Dirichlet/Laplacian branch",
            "repo_backing": ["P3101", "P3102", "P3089"],
            "legacy_analogy": "legacy observer/readout must remain downstream of nadsoliton/light/matter",
            "gates": [True, True, False, False, False, True],
            "missing_object": "non-imported Born-rule map, measurement basis source, detector calibration, empirical frequency interface",
        },
        {
            "physical_property": "light-before-matter ontology and observed radiation",
            "strict_slot": "nadsoliton -> light -> matter -> emergent observer ordering with strict operational kernels as controls",
            "repo_backing": ["AX9/K1 guardrail", "P3101/P3102 readout obstructions"],
            "legacy_analogy": "legacy ontology fits the order, but strict side must still source observation interfaces",
            "gates": [True, True, False, False, False, True],
            "missing_object": "observed-light interface and detector/readout calibration without apparatus import",
        },
        {
            "physical_property": "dimensionful units: length, action, time, calibration",
            "strict_slot": "entropy/alpha_geo/unit-map lanes, spectral geometry rows, V_lift/Upsilon unit-measure audits",
            "repo_backing": ["P2689", "P3104", "P3143", "P3144"],
            "legacy_analogy": "legacy numerical roles do not provide strict physical units by themselves",
            "gates": [True, True, False, False, False, True],
            "missing_object": "intrinsic bit-to-length/action/time conversion law",
        },
    ]


def classify(row: dict[str, Any]) -> dict[str, Any]:
    gate_map = dict(zip(GATES, row["gates"], strict=True))
    positive = sum(row["gates"])
    if all(row["gates"]):
        status = "closed"
    elif gate_map["strict_carrier_present"] and gate_map["receiver_scaffold_present"]:
        status = "receiver_layout_only"
    else:
        status = "speculative_slot_only"
    return {**row, "gate_map": gate_map, "positive_gate_count": positive, "layout_status": status}


def build_payload() -> dict[str, Any]:
    inputs = {name: load_text_or_json(path) for name, path in INPUTS.items()}
    rows = [classify(row) for row in physical_rows()]
    closed = [r for r in rows if r["layout_status"] == "closed"]
    receiver = [r for r in rows if r["layout_status"] == "receiver_layout_only"]
    theorem = {
        "name": "P3145_T1_strict_kernel_reverse_sm_gr_layout_obligation_matrix",
        "statement": "If known SM/GR physics is forced to fit inside the strict nadsoliton kernel picture, current artifacts place each physical property into receiver/scaffold slots, not into closed strict sources.  All 10 audited properties have a plausible strict carrier and receiver scaffold, but 0/10 pass the full source, unit/calibration, nonproxy dynamics, and selector-safety package required for ToE closure.",
        "finite_counts": {
            "physical_properties_audited": len(rows),
            "strict_carrier_rows": sum(r["gate_map"]["strict_carrier_present"] for r in rows),
            "receiver_scaffold_rows": sum(r["gate_map"]["receiver_scaffold_present"] for r in rows),
            "source_law_rows": sum(r["gate_map"]["source_law_present"] for r in rows),
            "unit_or_calibration_rows": sum(r["gate_map"]["unit_or_calibration_present"] for r in rows),
            "nonproxy_eom_rows": sum(r["gate_map"]["nonproxy_eom_or_dynamics_present"] for r in rows),
            "selector_resolved_rows": sum(r["gate_map"]["selector_obstruction_resolved_when_needed"] for r in rows),
            "closed_rows": len(closed),
            "receiver_layout_only_rows": len(receiver),
        },
    }
    return {
        "status": "P3145_STRICT_KERNEL_REVERSE_SM_GR_LAYOUT_MATRIX_RECEIVER_ONLY_NO_CLOSURE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_context": {name: (inputs[name] or {}).get("status") if isinstance(inputs[name], dict) else None for name in inputs},
        "constructed_object": {
            "name": "Xi_SMGR^strict reverse-layout obligation matrix",
            "classification": "axiomatic_reverse_mapping_from_observed_physics_to_strict_kernel_slots",
            "kernel_reading": "K_strict_gate is treated as completed/enriched operational strict kernel, not as silent legacy role transfer",
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
        },
        "layout_rows": rows,
        "finite_theorem": theorem,
        "decision": {
            "bounded_result": "P3145 gives the requested reverse analysis: SM/GR can be laid out over strict-kernel slots as a receiver architecture, but current artifacts do not close the source/unit/EOM/selector obligations.  The strict kernel picture is therefore a structured research map, not a ToE closure.",
            "axiomatic_route_recommendation": "Continue axiomatizing only if every axiom is labelled as conditional and used to derive downstream consequences.  For strict-core ToE progress, prioritize a non-premise source theorem for one missing bridge: either selector/orientation source, physical unit/action measure, or nonproxy metric/gauge EOM coupling.",
            "next_honest_step": "Construct one reverse-layout bridge theorem candidate for the most bottlenecked row: a non-premise physical-unit/action-measure source tying strict finite receivers to a dimensionful Lagrangian density.  If that cannot be built, pivot to a new non-premise selector/orientation source; do not claim SM/GR reduction or ToE from receiver layout alone.",
            "negative_export_flags": {
                "SM_GR_reduction_theorem_exported": False,
                "strict_selector_closure_exported": False,
                "unit_bearing_L_total_exported": False,
                "Einstein_Hilbert_derivation_exported": False,
                "gauge_uniqueness_exported": False,
                "legacy_role_transfer_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3145/S2095 Strict-kernel reverse SM/GR layout matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        f"- Ontology order: `{payload['constructed_object']['ontology_order']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Reverse-layout rows"])
    for row in payload["layout_rows"]:
        lines.append(f"- `{row['physical_property']}` -> `{row['strict_slot']}`; status `{row['layout_status']}`; missing `{row['missing_object']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Axiomatic-route recommendation", payload["decision"]["axiomatic_route_recommendation"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3145/S2095 strict-kernel reverse SM/GR layout matrix", "## P3145/S2095 strict-kernel reverse SM/GR layout matrix\n\n`P3145/S2095` constructs `Xi_SMGR^strict`, a reverse-layout obligation matrix asking where known SM/GR properties would have to live if they are to fit inside the strict nadsoliton/kernel picture.  All `10/10` audited properties have plausible strict carrier/receiver slots, but `0/10` close the full source-law, unit/calibration, nonproxy dynamics, and selector-safety package.  This is a structured receiver architecture for reverse analysis, not SM/GR reduction, legacy role transfer, `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3145/S2095 reverse SM/GR layout remains receiver-only", "## P3145/S2095 reverse SM/GR layout remains receiver-only\n\n`P3145/S2095` maps metric geometry, GR stress-energy, gauge scaffold, couplings, masses, chirality, RG damping, quantum readout, light/matter order, and physical units onto strict-kernel receiver slots.  None of the rows currently supplies a unit-bearing nonproxy Lagrangian/EOM term or full SM/GR reduction theorem.\n")
    append_once(AGENTS, "Current strict-kernel reverse SM/GR layout guardrail (P3145/S2095, 2026-07-13)", "## Current strict-kernel reverse SM/GR layout guardrail (P3145/S2095, 2026-07-13)\n\n- P3145 constructs `Xi_SMGR^strict`, a reverse-layout matrix for how SM/GR properties would have to sit inside the strict nadsoliton/kernel picture.\n- The audit finds `10/10` plausible strict carrier/receiver slots but `0/10` rows with the full source-law, unit/calibration, nonproxy dynamics, and selector-safety package.\n- Do not promote reverse layout, legacy analogies, or receiver scaffolds to SM/GR reduction, strict selector closure, unit-bearing `L_total`, bridge completion, role transfer, or ToE closure.\n- Axiomatic continuation is useful only as explicitly conditional model-building; strict-core progress should target one non-premise bridge theorem for physical units/action measure, selector/orientation source, or nonproxy metric/gauge EOM coupling.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
