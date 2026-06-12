#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2640_s1590_offset_stride_topology_selector_source_no_go_audit.json"
MD = GEN / "p2640_s1590_offset_stride_topology_selector_source_no_go_audit.md"

SOURCE_FILES = {
    "P2639_OFFSET_STRIDE_EXHAUSTION": GEN / "p2639_s1589_offset_stride_metric_lift_exhaustion_closure_path.json",
    "N450_Z12_CARRIER_REGULAR_ACTION": ROOT / "N450_CURRENT_FIRST_ACTUAL_STRICT_Z12_CARRIER_AND_REGULAR_ACTION_WITNESS_THEOREM.md",
    "N457_PHASE12_AUT_Z12_QUOTIENT_HOLONOMY_BOUNDARY": ROOT / "N457_CURRENT_FIRST_STRICT_PHASE_12_AUT_Z12_QUOTIENT_ORBIT_SPACE_TOPOLOGICAL_HOLONOMY_TRIVIALITY_BOUNDARY_THEOREM.md",
    "T48_PRIMORDIAL_PREORIENTATION_OMEGA_PHI_TARGET": ROOT / "T48_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_PRIMORDIAL_PREORIENTATION_OMEGA_PHI_TYPED_TRANSPORT_TARGET_SPEC.md",
    "F672_PROJECTIVE_SELECTOR_CLOSURE": ROOT / "F672_FIRST_EXPORTED_SELECTOR_CLOSURE_GLOBAL_C_V1_PROJECTIVE_STRICT_V1_PACKET.md",
    "QW2191_MODE_INDEX_UNIQUENESS_OBSTRUCTION": REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2191_MODE_INDEX_UNIQUENESS_OBSTRUCTION_THEOREM_GATE.md",
    "QW2195_AXIOM_AUGMENTED_GENERATION_MAPPING": REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2195_GENERATION_MAPPING_AXIOM_AUGMENTED_GATE.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "canonical_offset_stride_source_exported",
    "p2639_role_like_lift_promoted",
    "phase_frequency_node_gauge_certificate_exported",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
    "strict_kernel_full_kernel_claimed",
    "toe_closure_claimed",
    "positive_beta_renormalization_source_exported",
    "inverse_hierarchy_role_transfer_exported",
    "selector_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "blind_empirical_confirmation_claimed",
]

ROLE_LIKE_TARGETS = [(4, 3), (10, 6)]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:90]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "offset_stride_source_content": (
            "canonical offset|zero-lattice offset|stride|metric lift|node lift|offset/stride|"
            "integer node|node/gauge|phase-frequency|gauge separation"
        ),
        "topology_selector_carrier_content": (
            "Z_12|Z12|Phase_12|Aut_Z12|quotient carrier|parity character|regular action|"
            "selector closure|selector source|projective selector|mode index uniqueness"
        ),
        "transport_orientation_content": (
            "primordial preorientation|omega-phi|omega,phi|typed transport|orientation|oriented cycle|"
            "successor map|holonomy|Berry|directed orientation"
        ),
        "non_fit_guard_content": (
            "without retuning|no-retune|post-hoc|target-independent|source theorem|hidden selector|"
            "axiom-free|axiom augmented|candidate_only|non_strict"
        ),
        "toe_bridge_blocker_content": (
            "bridge completion|role-transfer|full kernel|ToE closure|QW-2191|role-bearing L_total|"
            "positive beta|inverse hierarchy|blind empirical"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for topology/selector offset-stride source search", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def latest_commit_audit() -> list[dict[str, Any]]:
    proc = subprocess.run(
        ["git", "log", "-n", "6", "--oneline", "--name-only"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"^[0-9a-f]{7,12} ", line):
            if current:
                rows.append(current)
            sha, subject = line.split(" ", 1)
            current = {"sha": sha, "subject": subject, "files": []}
        elif current is not None:
            current["files"].append(line)
    if current:
        rows.append(current)
    return rows


def phrase_flags(text: str, phrases: list[str]) -> dict[str, bool]:
    lower = text.lower()
    return {phrase: phrase.lower() in lower for phrase in phrases}


def topology_selector_source_ledger() -> list[dict[str, Any]]:
    definitions = [
        {
            "id": "z12_regular_action_carrier",
            "source": "N450_Z12_CARRIER_REGULAR_ACTION",
            "path": SOURCE_FILES["N450_Z12_CARRIER_REGULAR_ACTION"],
            "exported_numbers": [12],
            "possible_stride_numbers": [],
            "hard_limit_phrases": ["no generator/orientation uniqueness", "does not attempt to export"],
            "positive_phrases": ["Z_12", "regular action"],
            "verdict": "carrier_only_no_offset_stride_selector",
        },
        {
            "id": "phase12_aut_z12_quotient_parity",
            "source": "N457_PHASE12_AUT_Z12_QUOTIENT_HOLONOMY_BOUNDARY",
            "path": SOURCE_FILES["N457_PHASE12_AUT_Z12_QUOTIENT_HOLONOMY_BOUNDARY"],
            "exported_numbers": [12, 6, 2],
            "possible_stride_numbers": [6, 2],
            "hard_limit_phrases": ["holonomy", "successor map", "hidden selector slots", "does not prove"],
            "positive_phrases": ["Q_v1", "chi_parity_Z12_v1"],
            "verdict": "quotient_parity_only_no_canonical_k0_or_successor",
        },
        {
            "id": "primordial_preorientation_omega_phi_transport_target",
            "source": "T48_PRIMORDIAL_PREORIENTATION_OMEGA_PHI_TARGET",
            "path": SOURCE_FILES["T48_PRIMORDIAL_PREORIENTATION_OMEGA_PHI_TARGET"],
            "exported_numbers": [],
            "possible_stride_numbers": [],
            "hard_limit_phrases": ["future-only typed transport target", "does not decide", "actual typed transport"],
            "positive_phrases": ["primordial-preorientation", "OmegaPhi"],
            "verdict": "future_target_only_no_actual_offset_stride_source",
        },
        {
            "id": "projective_global_selector_closure",
            "source": "F672_PROJECTIVE_SELECTOR_CLOSURE",
            "path": SOURCE_FILES["F672_PROJECTIVE_SELECTOR_CLOSURE"],
            "exported_numbers": [],
            "possible_stride_numbers": [],
            "hard_limit_phrases": ["projective", "does not claim", "physical promotion", "directed orientation"],
            "positive_phrases": ["SelectorClosure", "well-definedness"],
            "verdict": "projective_ray_closure_no_sign_sensitive_metric_lift",
        },
        {
            "id": "mode_index_uniqueness_obstruction",
            "source": "QW2191_MODE_INDEX_UNIQUENESS_OBSTRUCTION",
            "path": SOURCE_FILES["QW2191_MODE_INDEX_UNIQUENESS_OBSTRUCTION"],
            "exported_numbers": [2],
            "possible_stride_numbers": [],
            "hard_limit_phrases": ["O(2) rotation freedom", "obstructed without extra symmetry breaking"],
            "positive_phrases": ["kernel eigenspaces"],
            "verdict": "strict_obstruction_not_selector_source",
        },
        {
            "id": "axiom_augmented_generation_mapping",
            "source": "QW2195_AXIOM_AUGMENTED_GENERATION_MAPPING",
            "path": SOURCE_FILES["QW2195_AXIOM_AUGMENTED_GENERATION_MAPPING"],
            "exported_numbers": [],
            "possible_stride_numbers": [],
            "hard_limit_phrases": ["axiom-augmented", "axiom-free physical uniqueness", "open"],
            "positive_phrases": ["deterministic", "generation mapping"],
            "verdict": "deterministic_only_after_axiom_not_strict_source",
        },
    ]
    rows = []
    for item in definitions:
        text = load_text(item["path"])
        hard = phrase_flags(text, item["hard_limit_phrases"])
        positive = phrase_flags(text, item["positive_phrases"])
        exports_offset = bool(re.search(r"\bk0\b|canonical offset|zero-lattice offset", text, re.IGNORECASE))
        exports_stride = bool(re.search(r"canonical stride|zero-lattice stride|\bm\s*=\s*(3|6)\b", text, re.IGNORECASE))
        strict_source_pass = exports_offset and exports_stride and not any(hard.values())
        rows.append({
            "id": item["id"],
            "source": item["source"],
            "path": rel(item["path"]),
            "exists": item["path"].exists(),
            "sha256": sha256_file(item["path"]),
            "exported_numbers": item["exported_numbers"],
            "possible_stride_numbers": item["possible_stride_numbers"],
            "positive_phrase_hits": positive,
            "hard_limit_phrase_hits": hard,
            "exports_canonical_zero_lattice_offset": exports_offset,
            "exports_canonical_zero_lattice_stride": exports_stride,
            "strict_offset_stride_source_pass": strict_source_pass,
            "verdict": item["verdict"],
        })
    return rows


def p2639_role_like_candidates() -> list[dict[str, Any]]:
    p2639 = load_json(SOURCE_FILES["P2639_OFFSET_STRIDE_EXHAUSTION"])
    return p2639.get("offset_stride_exhaustion", {}).get("role_like_candidates", [])


def compatibility_matrix(ledger: list[dict[str, Any]], role_like: list[dict[str, Any]]) -> list[dict[str, Any]]:
    matrix = []
    targets = [(row.get("k0"), row.get("stride")) for row in role_like] or ROLE_LIKE_TARGETS
    for source in ledger:
        for k0, stride in targets:
            numeric_stride_hint = stride in source["possible_stride_numbers"]
            matrix.append({
                "source_id": source["id"],
                "candidate_k0": k0,
                "candidate_stride": stride,
                "has_numeric_stride_hint": numeric_stride_hint,
                "has_canonical_offset": source["exports_canonical_zero_lattice_offset"],
                "has_canonical_stride": source["exports_canonical_zero_lattice_stride"],
                "licensed_without_hidden_selector": source["strict_offset_stride_source_pass"],
                "decision": "reject_as_source" if not source["strict_offset_stride_source_pass"] else "would_require_role_rerun",
            })
    return matrix


def closure_decision(ledger: list[dict[str, Any]], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    gates = {
        "p2639_role_like_lifts_exist": True,
        "repo_exports_topology_selector_carriers": any(row["exists"] for row in ledger),
        "some_source_exports_canonical_k0_and_stride": any(row["strict_offset_stride_source_pass"] for row in ledger),
        "some_role_like_lift_licensed_without_hidden_selector": any(row["licensed_without_hidden_selector"] for row in matrix),
        "q_w_2191_selector_obstruction_discharged": False,
        "role_transfer_rerun_passed_under_sourced_lift": False,
    }
    return {
        "gates": gates,
        "promote_offset_stride_to_bridge_completion": all(gates.values()),
        "full_kernel_now": False,
        "classification": "TOPOLOGY_SELECTOR_CARRIERS_EXIST_BUT_NO_CANONICAL_OFFSET_STRIDE_SOURCE_FOR_P2639_LIFTS",
        "professorial_verdict": (
            "The newest P2637-P2639 path makes the phase-node mismatch mathematically structured, but the current topology/selector layer still does not choose the needed zero-lattice offset and stride. "
            "Z12, quotient, parity, projective selector, primordial-preorientation, and axiom-augmented generation artifacts are real support objects; however they either remain carrier-level, projective, future-only, axiom-augmented, or explicitly obstructed by QW-2191. "
            "Thus the strict kernel remains a strong ToE-like working kernel, not a full kernel with completed legacy node/gauge role transfer."
        ),
        "next_honest_step": (
            "Do not search wider offset/stride boxes unless a new source atom appears.  The next proof-grade move is to construct a quotient-safe successor/connection/transport object on the Z12/Phase12 carrier that canonically exports a zero-lattice origin and stride, then rerun P2639 role-transfer gates; if that object cannot be exported without hidden selector input, demote the legacy integer node/gauge role and move to beta/inverse-hierarchy blockers."
        ),
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cl = payload["closure_decision"]
    lines = [
        "# P2640/S1590 offset-stride topology/selector source no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Latest-research synchronization",
        "",
        "The audit first checked the newest commits and treats P2637-P2639 as the live frontier: constructive phase-node repair, metric-pushforward viability failure, and offset/stride lift exhaustion.",
        "",
    ]
    for row in payload["latest_commit_audit"][:4]:
        lines.append(f"- `{row['sha']}` {row['subject']} ({len(row['files'])} touched files)")
    lines.extend([
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps research content about offset/stride source, topology/selector carriers, transport/orientation, non-fit guards, and ToE blockers before adding the finite source ledger.",
        "",
    ])
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Source ledger result",
        "",
        "| source atom | verdict | exports canonical offset? | exports canonical stride? | strict source pass? |",
        "| --- | --- | --- | --- | --- |",
    ])
    for row in payload["topology_selector_source_ledger"]:
        lines.append(
            f"| `{row['id']}` | `{row['verdict']}` | `{row['exports_canonical_zero_lattice_offset']}` | "
            f"`{row['exports_canonical_zero_lattice_stride']}` | `{row['strict_offset_stride_source_pass']}` |"
        )
    role_like = payload["p2639_role_like_candidates"]
    lines.extend([
        "",
        "## Compatibility with P2639 role-like lifts",
        "",
        f"P2639's role-like UV-safe lifts are `{[(row['k0'], row['stride']) for row in role_like]}`.",
        f"Licensed role-like source matches found: `{sum(1 for row in payload['compatibility_matrix'] if row['licensed_without_hidden_selector'])}`.",
        "",
        "A numeric echo is not enough: for example, quotient/parity artifacts contain a `6`, but they do not export the zero-lattice origin `k0=10` or a successor/connection law that would license stride `6` without a hidden selector.",
        "",
        "## Closure decision",
        "",
        cl["professorial_verdict"],
        "",
        f"Promote offset/stride source to bridge completion? `{cl['promote_offset_stride_to_bridge_completion']}`.",
        f"Full kernel now? `{cl['full_kernel_now']}`.",
        f"Classification: `{cl['classification']}`.",
        "",
        "## Recommended next honest step",
        "",
        cl["next_honest_step"],
        "",
        "No ToE closure, full-kernel finality, bridge completion, role-transfer, selector-source discharge, positive beta source, inverse-hierarchy transfer, blind empirical confirmation, or role-bearing `L_total` is claimed.",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    ledger = topology_selector_source_ledger()
    role_like = p2639_role_like_candidates()
    matrix = compatibility_matrix(ledger, role_like)
    cl = closure_decision(ledger, matrix)
    payload: dict[str, Any] = {
        "status": "P2640_OFFSET_STRIDE_TOPOLOGY_SELECTOR_SOURCE_NO_GO_AUDIT_NO_PROMOTION",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()},
        "p2639_role_like_candidates": role_like,
        "topology_selector_source_ledger": ledger,
        "compatibility_matrix": matrix,
        "closure_decision": cl,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2640/S1590 offset-stride topology/selector source no-go guard",
        "\n## P2640/S1590 offset-stride topology/selector source no-go guard\n\n"
        "`P2640/S1590` checks the newest P2637-P2639 frontier against topology/selector artifacts (`Z_12`, `Phase_12/Aut_Z12`, parity, primordial-preorientation, projective selector closure, `QW-2191`, axiom-augmented generation mapping).  These artifacts provide real carrier and selector-adjacent support, but none exports a canonical zero-lattice offset `k0` and stride `m` for the P2639 role-like lifts.  Therefore the phase-node repair remains a source-guarded candidate rather than bridge completion, and strict kernel full-finality, role-transfer, beta source, inverse hierarchy, `QW-2191`, role-bearing `L_total`, and ToE closure remain closed.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2640/S1590 offset-stride topology/selector source Ltotal guard",
        "\n## P2640/S1590 offset-stride topology/selector source Ltotal guard\n\n"
        "`P2640/S1590` does not re-enable `L_total`: current topology/selector carriers do not canonically choose the offset/stride metric lift required to transfer the legacy node/gauge role into strict dynamics.  A role-bearing term would first need a quotient-safe successor/connection/transport theorem that exports that origin and stride without hidden selector input.\n",
    )
    return payload


if __name__ == "__main__":
    main()
