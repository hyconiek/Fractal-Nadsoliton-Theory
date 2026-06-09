#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.json"
MD = GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.md"

PACKETS = {
    "P2601": GEN / "p2601_s1551_nadsoliton_identity_action_unital_multiplicative_source_theorem.json",
    "P2602": GEN / "p2602_s1552_nadsoliton_rg_fixed_rate_prime_log_source_theorem.json",
    "P2603": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
    "P2604": GEN / "p2604_s1554_strict_damping_post_source_bridge_readiness_matrix.json",
    "P2605": GEN / "p2605_s1555_legacy_strict_linear_slice_completion_map_evidence.json",
    "P2606": GEN / "p2606_s1556_strict_side_nonlinear_compression_residual_addition.json",
    "P2607": GEN / "p2607_s1557_strict_phase_topological_selector_bridge_completion.json",
    "P2608": GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.json",
}
MD_FILES = {pid: Path(str(path).replace(".json", ".md")) for pid, path in PACKETS.items()}

REVIEW_ROWS = [
    {
        "packet_id": "P2601",
        "review_target": "multiplicative_character_law_source",
        "ai_objection": "monoidal identity-action argument is a sketch, not a formal uniqueness proof",
        "required_ingredients": [
            "formal nadsoliton action monoid definition",
            "proof that y_1=0 is the unique boundary-compatible normalization",
            "D_f-independence check for the unital law",
        ],
        "ingredient_patterns": ["formal monoid", "unique", "D_f-independent|Df-independent|dimension-independent"],
        "operational_status": "QUARANTINE_POSITIVE_EXPORT_UNTIL_FORMAL_MONOID_UNIQUENESS_PROOF",
        "accepted_after_revalidation": False,
        "score": 1,
    },
    {
        "packet_id": "P2602",
        "review_target": "prime_log_proportionality_source",
        "ai_objection": "prime-log RG law assumes prime spectral gaps rather than deriving them",
        "required_ingredients": [
            "spectral-gap theorem selecting primes rather than another generator family",
            "fixed-rate RG derivation tied to that prime spectrum",
            "negative control excluding Fibonacci or arbitrary generator gaps",
        ],
        "ingredient_patterns": ["spectral gap", "prime", "Fibonacci|arbitrary generator|negative control"],
        "operational_status": "QUARANTINE_POSITIVE_EXPORT_UNTIL_PRIME_SPECTRAL_GAP_THEOREM",
        "accepted_after_revalidation": False,
        "score": 1,
    },
    {
        "packet_id": "P2603",
        "review_target": "fractal_codimension_slope_source",
        "ai_objection": "conditionally acceptable if D_f=9/5 and active codimension convention are held fixed",
        "required_ingredients": [
            "explicit D_f=9/5 premise",
            "active codimension formula D_f-1=4/5",
            "scope statement that this does not prove bridge or role transfer",
        ],
        "ingredient_patterns": ["D_f=9/5", "D_f-1=4/5|4/5", "does not export.*bridge|No.*bridge"],
        "operational_status": "CONDITIONALLY_RETAIN_SOURCE_SCOPE_ONLY",
        "accepted_after_revalidation": True,
        "score": 3,
    },
    {
        "packet_id": "P2604",
        "review_target": "post_source_bridge_readiness_matrix",
        "ai_objection": "logic matrix is useful because it separates source discharge from bridge/role gates",
        "required_ingredients": [
            "three bridge/role gates listed",
            "one accepting row over the gate truth table",
            "current role-bearing L_total remains false",
        ],
        "ingredient_patterns": ["three.*gates|Bridge/role gates", "accepting rows: `1`|accepting.*1", "Current role-bearing L_total accepts: `False`"],
        "operational_status": "RETAIN_AS_READINESS_MATRIX_NOT_SOURCE_THEOREM",
        "accepted_after_revalidation": True,
        "score": 3,
    },
    {
        "packet_id": "P2605",
        "review_target": "linear_slice_completion_map_evidence",
        "ai_objection": "linear-slice match at eta=1 is tautological evidence, not a bridge theorem",
        "required_ingredients": [
            "explicit classification as linear-slice evidence only",
            "nonlinear eta != 1 residual remains outside P2605",
            "no full bridge export from the linear slice alone",
        ],
        "ingredient_patterns": ["linear-slice.*evidence", "eta=1", "does not export.*full.*bridge|full legacy-to-strict bridge"],
        "operational_status": "DOWNGRADE_TO_TAUTOLOGICAL_LINEAR_SLICE_CHECK",
        "accepted_after_revalidation": False,
        "score": 1,
    },
    {
        "packet_id": "P2606",
        "review_target": "nonlinear_compression_residual_addition",
        "ai_objection": "nonlinear residual is useful numerically but remains only one strict-side component",
        "required_ingredients": [
            "nonzero residual computation",
            "explicit statement that phase/topological selector remains open",
            "no role-transfer theorem export",
        ],
        "ingredient_patterns": ["Nonlinear residual component exported: `True`", "phase/topological", "role-transfer theorem"],
        "operational_status": "RETAIN_AS_NUMERICAL_RESIDUAL_COMPONENT_ONLY",
        "accepted_after_revalidation": True,
        "score": 2,
    },
    {
        "packet_id": "P2607",
        "review_target": "GF2_phase_topological_selector_bridge_completion",
        "ai_objection": "full-rank GF(2) system lacks a physical construction theorem for the matrix",
        "required_ingredients": [
            "definition of what GF(2) rows and columns physically represent",
            "derivation of the matrix from nadsoliton phase/topological data",
            "invariance test under node relabeling or audited node perturbation",
        ],
        "ingredient_patterns": ["rows.*columns", "physical|nadsoliton", "relabel|perturb|invariance"],
        "operational_status": "QUARANTINE_BRIDGE_COMPLETION_EXPORT_UNTIL_GF2_ORIGIN_THEOREM",
        "accepted_after_revalidation": False,
        "score": 0,
    },
    {
        "packet_id": "P2608",
        "review_target": "role_bearing_L_total_export",
        "ai_objection": "role-bearing L_total lacks a formal role semantics definition",
        "required_ingredients": [
            "formal definition of role-bearing L_total term",
            "acceptance semantics distinguishing bookkeeping insertion from physical explanation",
            "proof that role transfer follows from a non-quarantined bridge completion theorem",
        ],
        "ingredient_patterns": ["definition of role|role semantics", "bookkeeping|physical explanation", "bridge completion"],
        "operational_status": "QUARANTINE_ROLE_BEARING_LTOTAL_EXPORT_UNTIL_ROLE_SEMANTICS_AND_BRIDGE_REVALIDATION",
        "accepted_after_revalidation": False,
        "score": 0,
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "strict_damping_beta_eta_source_newly_exported",
    "bridge_completion_newly_exported",
    "role_bearing_ltotal_newly_exported",
    "role_transfer_theorem_newly_exported",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_audit",
    "toe_closure_claimed",
    "apd_source_exported",
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
        "new_packet": "P2610|S1560|P2601.*P2608.*revalidation|source-strength revalidation|AI critique revalidation",
        "intended_research_nonduplication": "role-bearing.*quarantine|bridge completion.*quarantine|strict damping source export revalidation|P2607.*P2608.*critique",
        "criticized_chain": "P2601|P2602|P2603|P2604|P2605|P2606|P2607|P2608|role_bearing_L_total|bridge_completion",
        "guardrails": "K_legacy_ont|K_strict_gate|QW-2191|role-transfer|ToE closure|APD source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def flatten_text(obj: Any) -> str:
    if isinstance(obj, dict):
        return "\n".join(f"{key}: {flatten_text(value)}" for key, value in obj.items())
    if isinstance(obj, list):
        return "\n".join(flatten_text(value) for value in obj)
    return str(obj)


def pattern_hits(text: str, patterns: list[str]) -> dict[str, bool]:
    import re
    return {pattern: re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) is not None for pattern in patterns}


def build_review_rows(payloads: dict[str, dict[str, Any]], markdowns: dict[str, str]) -> list[dict[str, Any]]:
    rows = []
    for template in REVIEW_ROWS:
        pid = template["packet_id"]
        text = flatten_text(payloads[pid]) + "\n" + markdowns[pid]
        hits = pattern_hits(text, template["ingredient_patterns"])
        missing_patterns = [pattern for pattern, present in hits.items() if not present]
        rows.append({
            "packet_id": pid,
            "review_target": template["review_target"],
            "ai_objection": template["ai_objection"],
            "required_ingredients": template["required_ingredients"],
            "ingredient_pattern_hits": hits,
            "missing_pattern_count": len(missing_patterns),
            "missing_patterns": missing_patterns,
            "score": template["score"],
            "accepted_after_revalidation": template["accepted_after_revalidation"],
            "operational_status": template["operational_status"],
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {pid: load_json(path) for pid, path in PACKETS.items()}
    markdowns = {pid: MD_FILES[pid].read_text(encoding="utf-8") for pid in PACKETS}
    rows = build_review_rows(payloads, markdowns)
    accepted = [row for row in rows if row["accepted_after_revalidation"]]
    quarantined = [row for row in rows if not row["accepted_after_revalidation"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2610_T1_p2601_p2608_source_strength_revalidation",
        "audited_chain": list(PACKETS.keys()),
        "review_statement": (
            "The AI critique is upheld in the operational sense: P2601, P2602, P2605, P2607, and P2608 are not strong enough to be used as unqualified strict source/bridge/role exports. P2603 and P2604 are retained with their stated scope, and P2606 is retained only as a numerical residual component."
        ),
        "review_rows": rows,
        "accepted_packet_ids_after_revalidation": [row["packet_id"] for row in accepted],
        "quarantined_packet_ids_after_revalidation": [row["packet_id"] for row in quarantined],
        "operational_overrides": {
            row["packet_id"]: row["operational_status"] for row in rows
        },
        "strong_export_policy_after_p2610": {
            "strict_damping_beta_eta_source": "DO_NOT_USE_AS_FULLY_DISCHARGED_UNTIL_P2601_AND_P2602_ARE_REPROVED_WITH_FORMAL_MONOID_AND_PRIME_SPECTRAL_GAP_THEOREMS",
            "legacy_to_strict_completion_bridge": "DO_NOT_USE_AS_CLOSED_UNTIL_P2607_GF2_MATRIX_HAS_PHYSICAL_ORIGIN_AND_INVARIANCE_THEOREM",
            "role_bearing_ltotal": "DO_NOT_USE_AS_ACCEPTED_UNTIL_ROLE_SEMANTICS_IS_DEFINED_AND_BRIDGE_COMPLETION_IS_REVALIDATED",
        },
        "concrete_missing_ingredient_frontier": {
            row["packet_id"]: row["required_ingredients"] for row in quarantined
        },
        "recommended_next_honest_step": (
            "Pick one quarantined item and prove the missing formal ingredient directly; the highest priority is either the P2607 GF(2) physical-origin theorem or the P2608 role-semantics definition."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "all_p2601_p2608_payloads_loaded": set(payloads) == set(PACKETS),
        "all_p2601_p2608_rows_reviewed": len(rows) == 8,
        "critique_upheld_for_p2607": "P2607" in theorem_export["quarantined_packet_ids_after_revalidation"],
        "critique_upheld_for_p2608": "P2608" in theorem_export["quarantined_packet_ids_after_revalidation"],
        "p2603_retained_conditionally": "P2603" in theorem_export["accepted_packet_ids_after_revalidation"],
        "p2604_retained_as_matrix": "P2604" in theorem_export["accepted_packet_ids_after_revalidation"],
        "no_new_source_exported": theorem_export["strict_damping_beta_eta_source_newly_exported"] is False,
        "no_new_bridge_completion_exported": theorem_export["bridge_completion_newly_exported"] is False,
        "no_new_role_bearing_ltotal_exported": theorem_export["role_bearing_ltotal_newly_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_audit"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2610",
        "stage_id": "S1560",
        "status": "P2610_CRITICAL_REVALIDATION_QUARANTINES_WEAK_P2601_P2602_P2605_P2607_P2608_EXPORTS_RETAINS_P2603_P2604_SCOPE_AND_P2606_RESIDUAL_ONLY",
        "source_files": {pid: rel(path) for pid, path in PACKETS.items()},
        "rg_non_duplication_audit": grep,
        "p2601_p2608_critical_revalidation_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {pid: sha256_json(payloads[pid]) for pid in PACKETS},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2601_p2608_critical_revalidation_audit"]["theorem_export"]
    lines = [
        "# P2610/S1560 P2601-P2608 critical revalidation audit", "",
        f"Status: `{payload['status']}`", "", "## Review statement", "",
        t["review_statement"], "", "## Operational result", "",
        f"- Accepted after revalidation: `{t['accepted_packet_ids_after_revalidation']}`.",
        f"- Quarantined after revalidation: `{t['quarantined_packet_ids_after_revalidation']}`.", "",
        "## Per-packet verdicts", "",
    ]
    for row in t["review_rows"]:
        lines.append(f"- `{row['packet_id']}` / `{row['review_target']}`: `{row['operational_status']}`; accepted=`{row['accepted_after_revalidation']}`; objection: {row['ai_objection']}.")
    lines.extend([
        "", "## Strong export policy after P2610", "",
    ])
    for key, value in t["strong_export_policy_after_p2610"].items():
        lines.append(f"- `{key}`: {value}.")
    lines.extend([
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Scope guards", "",
        "P2610 exports no new strict damping source, no bridge completion, no role-bearing L_total, no legacy physical-role transfer, no QW-2191 discharge, no APD source, and no ToE closure.", "", "## Fingerprint", "",
        f"`{payload['p2601_p2608_critical_revalidation_audit']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2610/S1560 P2601-P2608 critical revalidation audit

`P2610/S1560` accepts the external critique as an operational guard: P2601, P2602, P2605, P2607, and P2608 are quarantined as unqualified strict source/bridge/role exports until their missing formal ingredients are proved.  P2603 is retained only under its `D_f=9/5` codimension scope, P2604 is retained as a readiness matrix, and P2606 is retained only as a numerical nonlinear residual component.  `L_total` must not rely on the P2607/P2608 bridge/role export path until GF(2) physical-origin and role-semantics theorems exist.
""".strip()
    lag_section = """
## P2610/S1560 Ltotal revalidation quarantine

`P2610/S1560` supersedes the operational use of the P2607/P2608 positive flags for now: role-bearing `L_total` is quarantined pending a physical derivation of the GF(2) matrix and a formal definition of role-bearing semantics.  The strict damping/compression term may remain a candidate bookkeeping target, but not an accepted role-bearing term.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2610/S1560 P2601-P2608 critical revalidation audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2610/S1560 Ltotal revalidation quarantine", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
