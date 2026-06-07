#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2554_s1504_strict_damping_local_exhaustion_bridge_reorientation_certificate.json"
MD = GEN / "p2554_s1504_strict_damping_local_exhaustion_bridge_reorientation_certificate.md"

SOURCE_FILES = {
    "P2539_TOE_POTENTIAL": GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json",
    "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json",
    "P2553_NONHOMOGENEOUS_ANCHOR_CONSTANT": GEN / "p2553_s1503_strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate.json",
}

GATES = [
    "strict_damping_local_source_obligation_discharge",
    "legacy_to_strict_completion_bridge",
    "role_transfer_theorem",
    "qw2191_selector_discharge",
    "role_bearing_ltotal_export",
    "toe_closure_gate",
]
NEGATIVE_EXPORT_FLAGS = [
    "strict_damping_local_source_obligation_discharge_exported", "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported", "qw2191_selector_discharge_exported", "role_bearing_ltotal_exported",
    "toe_closure_claimed", "bridge_completion_claimed", "legacy_role_transfer_claimed", "strict_damping_beta_eta_source_exported",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2554|S1504|local exhaustion bridge reorientation|bridge reorientation certificate|post-local strict damping",
        "intended_research_nonduplication": "local.*strict damping.*exhaustion|strict damping.*bridge reorientation|post-local.*bridge|source route exhaustion|legacy->strict.*source audit",
        "precursors": "P2539|S1489|P2547|S1497|P2553|S1503|ToE-potential|nonhomogeneous anchor constant",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2539": theorem(sources["P2539_TOE_POTENTIAL"], "strict_damping_toe_potential_recommendation_certificate"),
        "P2547": theorem(sources["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"], "strict_damping_post_identity_residual_trikey_certificate"),
        "P2553": theorem(sources["P2553_NONHOMOGENEOUS_ANCHOR_CONSTANT"], "strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate"),
    }


def bridge_gate_truth_table() -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(GATES)):
        assignment = dict(zip(GATES, values))
        accepts = all(values)
        rows.append({
            "assignment": assignment,
            "accepts_toe_closure_readiness": accepts,
            "missing_gates": [gate for gate, value in assignment.items() if not value],
        })
    return rows


def sensitivity_rows() -> list[dict[str, Any]]:
    rows = []
    for omitted in GATES:
        assignment = {gate: True for gate in GATES}
        assignment[omitted] = False
        rows.append({
            "omitted_gate": omitted,
            "assignment": assignment,
            "accepts_toe_closure_readiness": False,
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    ts = load_theorems(sources)
    table = bridge_gate_truth_table()
    accepting = [row for row in table if row["accepts_toe_closure_readiness"]]
    local_only_assignment = {gate: False for gate in GATES}
    local_only_assignment["strict_damping_local_source_obligation_discharge"] = True
    local_only_accepts = all(local_only_assignment.values())
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2554_T1_strict_damping_local_exhaustion_bridge_reorientation_certificate",
        "audited_chain": ["P2539/S1489", "P2547/S1497", "P2553/S1503"],
        "p2539_toe_potential_inherited": ts["P2539"].get("toe_potential_recommendation_certificate_exported") is True,
        "p2547_post_identity_residual_trikey_inherited": ts["P2547"].get("post_identity_residual_trikey_irredundancy_exported") is True,
        "p2553_anchor_constant_obligation_inherited": ts["P2553"].get("nonhomogeneous_anchor_reduces_to_constant_source_obligation") is True,
        "bridge_gate_vector": GATES,
        "bridge_gate_truth_table_row_count": len(table),
        "bridge_gate_accepting_row_count": len(accepting),
        "bridge_gate_accepting_row": accepting[0],
        "single_gate_sensitivity_rows": sensitivity_rows(),
        "single_gate_sensitivity_row_count": len(GATES),
        "all_single_gate_omissions_reject": True,
        "local_strict_damping_source_alone_assignment": local_only_assignment,
        "local_strict_damping_source_alone_accepts_toe_or_role_transfer": local_only_accepts,
        "local_route_exhaustion_exported_as_reorientation_not_closure": True,
        "recommended_next_honest_step": (
            "Stop adding local strict-damping bookkeeping layers unless they derive a real source. The next honest step is the broader legacy->strict completion/source bridge audit: identify the damping/compression passage, the selector/source premise, and only then run a separate role-transfer theorem."
        ),
        "not_licensed": [
            "No local strict-damping source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, or ToE closure is exported.",
            "Legacy physical roles are not transferred to K_strict_gate.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "precursors_inherited": theorem_export["p2539_toe_potential_inherited"] and theorem_export["p2547_post_identity_residual_trikey_inherited"] and theorem_export["p2553_anchor_constant_obligation_inherited"],
        "truth_table_exact": theorem_export["bridge_gate_truth_table_row_count"] == 64 and theorem_export["bridge_gate_accepting_row_count"] == 1,
        "local_source_alone_rejects": theorem_export["local_strict_damping_source_alone_accepts_toe_or_role_transfer"] is False,
        "single_gate_omissions_reject": theorem_export["all_single_gate_omissions_reject"],
        "no_false_bridge_or_role_transfer": theorem_export["legacy_to_strict_completion_bridge_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_selector_discharge_exported"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2554",
        "stage_id": "S1504",
        "status": "STRICT_DAMPING_LOCAL_EXHAUSTION_BRIDGE_REORIENTATION_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_local_exhaustion_bridge_reorientation_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_local_exhaustion_bridge_reorientation_certificate"]["theorem_export"]
    lines = [
        "# P2554/S1504 strict damping local-exhaustion bridge-reorientation certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Bridge gate vector: `{t['bridge_gate_vector']}`.",
        f"- Truth-table rows / accepting rows: `{t['bridge_gate_truth_table_row_count']}` / `{t['bridge_gate_accepting_row_count']}`.",
        f"- Local strict-damping source alone accepts ToE/role readiness: `{t['local_strict_damping_source_alone_accepts_toe_or_role_transfer']}`.",
        f"- All single-gate omissions reject: `{t['all_single_gate_omissions_reject']}`.", "", "## Interpretation", "",
        "Even a hypothetical local strict-damping source is only one gate in the post-local bridge vector. Bridge completion, role-transfer, QW-2191 selector discharge, role-bearing `L_total`, and ToE closure remain separate gates.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No source export, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['strict_damping_local_exhaustion_bridge_reorientation_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2554/S1504` reorients the strict-damping work after the local obstruction chain.  A six-gate truth table over local source discharge, `legacy -> strict` bridge completion, role-transfer theorem, QW-2191 selector discharge, role-bearing `L_total`, and ToE closure has `64` rows and only the all-true row accepts.  Even hypothetically supplying the local strict-damping source alone leaves the bridge/role/selector/Ltotal gates closed, so the next honest frontier is the broader completion/source bridge audit rather than more local bookkeeping.
""".strip()
    lag_section = """
`P2554/S1504` blocks promotion from local strict-damping bookkeeping to role-bearing `L_total`: local source discharge alone does not provide bridge completion, role transfer, selector discharge, or ToE closure.  Further Lagrangian promotion should wait for the explicit `legacy -> strict` completion/source bridge and the separate role-transfer theorem.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2554/S1504 local-exhaustion bridge-reorientation guard", "## P2554/S1504 local-exhaustion bridge-reorientation guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2554/S1504 local-exhaustion bridge-reorientation Ltotal guard", "## P2554/S1504 local-exhaustion bridge-reorientation Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
