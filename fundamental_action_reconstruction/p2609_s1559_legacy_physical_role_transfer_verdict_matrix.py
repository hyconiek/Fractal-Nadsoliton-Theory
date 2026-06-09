#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2609_s1559_legacy_physical_role_transfer_verdict_matrix.json"
MD = GEN / "p2609_s1559_legacy_physical_role_transfer_verdict_matrix.md"

SOURCE_FILES = {
    "P2607_BRIDGE_COMPLETION": GEN / "p2607_s1557_strict_phase_topological_selector_bridge_completion.json",
    "P2608_STRICT_DAMPING_ROLE_TRANSFER": GEN / "p2608_s1558_strict_damping_role_transfer_to_ltotal_theorem.json",
}

ROLE_CLAIMS = [
    {
        "claim_id": "legacy_electroweak_sin2_theta_w",
        "legacy_formula": "sin^2(theta_W)=alpha_geo/12",
        "legacy_search_pattern": "sin\\^2\\(theta_W\\)|sin2_theta|theta_W",
        "required_strict_ingredient": "strict electroweak gauge-coupling normalization map from completed K_strict_gate parameters to theta_W",
        "strict_transfer_theorem_available": False,
        "strict_numeric_replay_available": False,
        "strict_role_map_available": False,
        "selector_obstruction_relevant": False,
    },
    {
        "claim_id": "legacy_alpha_em_inverse_formula",
        "legacy_formula": "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "legacy_search_pattern": "alpha_EM|alpha_em|fine structure|fine-structure",
        "required_strict_ingredient": "strict electromagnetic coupling map plus beta_tors-to-strict-damping parameter transport law",
        "strict_transfer_theorem_available": False,
        "strict_numeric_replay_available": False,
        "strict_role_map_available": False,
        "selector_obstruction_relevant": False,
    },
    {
        "claim_id": "legacy_beta_power_gravity_hierarchy",
        "legacy_formula": "gravity hierarchy from beta^N scaling",
        "legacy_search_pattern": "gravity hierarchy|beta\\^N|beta\\^n|hierarchy",
        "required_strict_ingredient": "strict hierarchy exponent source after nonlinear compression and strict damping role-bearing L_total insertion",
        "strict_transfer_theorem_available": False,
        "strict_numeric_replay_available": False,
        "strict_role_map_available": False,
        "selector_obstruction_relevant": False,
    },
    {
        "claim_id": "legacy_beta_tors_chi11_orientation",
        "legacy_formula": "beta_tors -> chi_11 orientation/selector role",
        "legacy_search_pattern": "beta_tors.*chi|chi11|chi_11|orientation",
        "required_strict_ingredient": "internal chi_11 selector source or QW-2191 discharge independent of beta_tors lore",
        "strict_transfer_theorem_available": False,
        "strict_numeric_replay_available": False,
        "strict_role_map_available": False,
        "selector_obstruction_relevant": True,
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "legacy_physical_role_transfer_exported",
    "legacy_electroweak_role_exported",
    "legacy_alpha_em_role_exported",
    "legacy_gravity_hierarchy_role_exported",
    "legacy_beta_tors_orientation_role_exported",
    "qw2191_discharged_by_this_packet",
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
    return {"count": len(lines), "samples": lines[:60]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2609|S1559|legacy physical-role transfer verdict|legacy physical role transfer verdict|legacy role transfer classification",
        "intended_research_nonduplication": "post-P2608 legacy role audit|legacy physical-role verdict matrix|electroweak role transfer|alpha_EM role transfer|gravity hierarchy role transfer|beta_tors orientation transfer",
        "source_chain": "P2607|S1557|P2608|S1558|legacy_to_strict_completion_bridge_exported|strict_damping_role_transfer_theorem_exported",
        "guardrails": "sin\\^2\\(theta_W\\)|alpha_EM|gravity hierarchy|beta_tors|QW-2191|role-transfer audit",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def classify_claims() -> list[dict[str, Any]]:
    rows = []
    for claim in ROLE_CLAIMS:
        legacy_hits = rg_count(claim["legacy_search_pattern"])
        necessary_conditions = {
            "legacy_corpus_mentions_claim": legacy_hits["count"] > 0,
            "strict_transfer_theorem_available_for_this_claim": claim["strict_transfer_theorem_available"],
            "strict_numeric_replay_available_for_this_claim": claim["strict_numeric_replay_available"],
            "strict_role_map_available_for_this_claim": claim["strict_role_map_available"],
            "qw2191_not_blocking_this_claim": not claim["selector_obstruction_relevant"],
        }
        transfer_accepts = all(necessary_conditions.values())
        missing = [key for key, value in necessary_conditions.items() if not value]
        rows.append({
            "claim_id": claim["claim_id"],
            "legacy_formula": claim["legacy_formula"],
            "legacy_corpus_hit_count": legacy_hits["count"],
            "legacy_corpus_samples": legacy_hits["samples"][:10],
            "required_strict_ingredient": claim["required_strict_ingredient"],
            "necessary_conditions": necessary_conditions,
            "missing_condition_count": len(missing),
            "missing_conditions": missing,
            "transfer_accepts": transfer_accepts,
            "verdict": "REJECT_STRICT_TRANSFER_NOW" if not transfer_accepts else "ACCEPT_STRICT_TRANSFER",
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2607_payload = load_json(SOURCE_FILES["P2607_BRIDGE_COMPLETION"])
    p2608_payload = load_json(SOURCE_FILES["P2608_STRICT_DAMPING_ROLE_TRANSFER"])
    p2607 = theorem(p2607_payload, "strict_phase_topological_selector_bridge_completion")
    p2608 = theorem(p2608_payload, "strict_damping_role_transfer_to_ltotal_theorem")
    claim_rows = classify_claims()
    rejected = [row for row in claim_rows if not row["transfer_accepts"]]
    accepted = [row for row in claim_rows if row["transfer_accepts"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2609_T1_legacy_physical_role_transfer_verdict_matrix",
        "audited_chain": ["P2607/S1557", "P2608/S1558"],
        "legacy_to_strict_completion_bridge_inherited": p2607.get("legacy_to_strict_completion_bridge_exported") is True,
        "strict_damping_role_transfer_inherited": p2608.get("strict_damping_role_transfer_theorem_exported") is True,
        "role_bearing_ltotal_for_strict_damping_inherited": p2608.get("role_bearing_ltotal_exported") is True,
        "theorem_statement": (
            "After P2607 bridge completion and P2608 scoped strict-damping role transfer, no listed legacy physical-role formula transfers to K_strict_gate. Each legacy role lacks at least one claim-specific strict role map, numeric replay, transfer theorem, or selector discharge; therefore strict damping/compression survives as the only current role-bearing transfer."
        ),
        "legacy_claim_verdict_rows": claim_rows,
        "legacy_claim_count": len(claim_rows),
        "accepted_legacy_claim_count": len(accepted),
        "rejected_legacy_claim_count": len(rejected),
        "survives_without_modification": [],
        "modified_survivals": [{
            "role": "strict_damping_compression_beta_eta_term",
            "source": "P2608/S1558",
            "modification": "role-bearing transfer is scoped to strict damping/compression only, not to any legacy electroweak, alpha_EM, gravity-hierarchy, or beta_tors orientation formula",
        }],
        "rejected_legacy_claim_ids": [row["claim_id"] for row in rejected],
        "concrete_missing_ingredients": {row["claim_id"]: row["required_strict_ingredient"] for row in rejected},
        "recommended_next_honest_step": (
            "Choose exactly one rejected legacy physical role and build a claim-specific strict role map plus numeric replay; do not infer it from P2608 damping role transfer."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2607_bridge_completion_inherited": theorem_export["legacy_to_strict_completion_bridge_inherited"],
        "p2608_strict_damping_role_transfer_inherited": theorem_export["strict_damping_role_transfer_inherited"],
        "strict_damping_ltotal_inherited": theorem_export["role_bearing_ltotal_for_strict_damping_inherited"],
        "all_legacy_claims_classified": theorem_export["legacy_claim_count"] == len(ROLE_CLAIMS),
        "no_legacy_claims_accepted": theorem_export["accepted_legacy_claim_count"] == 0,
        "all_legacy_claims_rejected_now": theorem_export["rejected_legacy_claim_count"] == len(ROLE_CLAIMS),
        "each_rejection_has_concrete_missing_ingredient": all(row["required_strict_ingredient"] for row in rejected),
        "no_legacy_physical_role_transfer_exported": theorem_export["legacy_physical_role_transfer_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2609",
        "stage_id": "S1559",
        "status": "P2609_LEGACY_PHYSICAL_ROLE_TRANSFER_AUDIT_REJECTS_ALL_LISTED_LEGACY_CLAIMS_STRICT_DAMPING_SCOPE_ONLY",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_physical_role_transfer_verdict_matrix": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2607_BRIDGE_COMPLETION": sha256_json(p2607_payload),
                "P2608_STRICT_DAMPING_ROLE_TRANSFER": sha256_json(p2608_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_physical_role_transfer_verdict_matrix"]["theorem_export"]
    lines = [
        "# P2609/S1559 legacy physical-role transfer verdict matrix", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        t["theorem_statement"], "", "## Verdict counts", "",
        f"- Legacy claims audited: `{t['legacy_claim_count']}`.",
        f"- Accepted legacy claims: `{t['accepted_legacy_claim_count']}`.",
        f"- Rejected legacy claims now: `{t['rejected_legacy_claim_count']}`.",
        f"- Strict damping L_total inherited: `{t['role_bearing_ltotal_for_strict_damping_inherited']}`.", "",
        "## Claim verdicts", "",
    ]
    for row in t["legacy_claim_verdict_rows"]:
        lines.append(f"- `{row['claim_id']}`: `{row['verdict']}`; missing ingredient: {row['required_strict_ingredient']}.")
    lines.extend([
        "", "## Scope guards", "",
        "P2609 does not transfer legacy electroweak, alpha_EM, gravity hierarchy, or beta_tors orientation roles; it does not discharge QW-2191, APD sourcehood, or ToE closure.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Fingerprint", "",
        f"`{payload['legacy_physical_role_transfer_verdict_matrix']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2609/S1559 legacy physical-role transfer verdict matrix

`P2609/S1559` performs the mandatory post-bridge legacy-role audit after P2607/P2608.  It classifies the listed legacy formulas (`sin^2(theta_W)`, `alpha_EM`, gravity hierarchy, `beta_tors` orientation) as not transferable to `K_strict_gate` on the current evidence: each lacks a claim-specific strict role map, numeric replay, transfer theorem, or selector discharge.  The only role-bearing transfer inherited is the scoped strict damping/compression term from P2608.
""".strip()
    lag_section = """
## P2609/S1559 legacy-role Ltotal non-transfer guard

`P2609/S1559` keeps `L_total` role-bearing status scoped to the strict damping/compression beta/eta term.  Legacy electroweak, `alpha_EM`, gravity-hierarchy, and `beta_tors` orientation roles are rejected for strict transfer until separate claim-specific role maps and numeric replays are supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2609/S1559 legacy physical-role transfer verdict matrix", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2609/S1559 legacy-role Ltotal non-transfer guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
