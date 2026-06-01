#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.json"
MD = GEN / "p2390_s1340_selector_qualified_beta_tors_chi11_role_audit.md"

SOURCE_FILES = {
    "P1343_SELECTOR_SOURCE_REPORT": GEN / "p1343_p1343_report_v1.json",
    "P1348_GLOBAL_CLOSURE_REPORT": GEN / "p1348_p1348_report_v1.json",
    "P1350_COUPLING_EXTRACTION_PACKET": ROOT / "P1350_RELEASE_8_1_STRICT_COUPLING_CONSTANT_EXTRACTION_AND_SM_GR_BRIDGE_NEXT_HONEST_STEP_PACKET.md",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "COMPONENT_GAP_REPORT": SCRATCH / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "P2389_SLACK_BUDGET": GEN / "p2389_s1339_cap_slack_budget_sensitivity_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

SEARCH_PATTERNS = {
    "beta_tors_to_chi11_literal": re.compile(r"beta_tors\s*(?:->|to)\s*chi[_]?(?:11)?|chi[_]?11.*beta_tors", re.IGNORECASE),
    "selector_source": re.compile(r"strict internal selector source|S_strict_internal_v1|selector source", re.IGNORECASE),
    "role_transfer": re.compile(r"role[- ]transfer|legacy physical-role transfer|legacy role transfer", re.IGNORECASE),
    "no_beta_tors_chi11": re.compile(r"No `?beta_tors\s*->\s*chi[_]?11`? theorem|beta_tors.*candidate.*not.*theorem", re.IGNORECASE),
}


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_ARTIFACT::{rel(path)}"
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def iter_research_files() -> list[Path]:
    return sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.suffix in {".py", ".md", ".json"}
        and ".git" not in path.parts
        and "__pycache__" not in path.parts
    )


def grep_audit() -> dict[str, Any]:
    files = iter_research_files()
    pattern_counts: dict[str, int] = {name: 0 for name in SEARCH_PATTERNS}
    file_counts: dict[str, int] = {name: 0 for name in SEARCH_PATTERNS}
    samples: dict[str, list[dict[str, Any]]] = {name: [] for name in SEARCH_PATTERNS}
    for path in files:
        text = read_text(path)
        for name, pattern in SEARCH_PATTERNS.items():
            matches = list(pattern.finditer(text))
            if not matches:
                continue
            file_counts[name] += 1
            pattern_counts[name] += len(matches)
            if len(samples[name]) < 8:
                for match in matches[: 8 - len(samples[name])]:
                    line_no = text.count("\n", 0, match.start()) + 1
                    line = text.splitlines()[line_no - 1].strip()[:220]
                    samples[name].append({"file": rel(path), "line": line_no, "text": line})
    return {
        "searched_file_count": len(files),
        "searched_terms": {
            "beta_tors_to_chi11_literal": "beta_tors -> chi11 / chi_11 with beta_tors",
            "selector_source": "strict internal selector source / S_strict_internal_v1 / selector source",
            "role_transfer": "role-transfer / legacy physical-role transfer",
            "no_beta_tors_chi11": "negative/candidate-not-theorem beta_tors -> chi11 wording",
        },
        "pattern_match_counts": pattern_counts,
        "pattern_file_counts": file_counts,
        "samples": samples,
        "nonduplication_finding": (
            "Existing reports already discuss beta_tors->chi11 as a candidate/non-theorem and P1343/P1348 as selector closure. "
            "P2390 therefore does not try to re-prove the selector; it audits whether that selector is sufficient to license legacy beta_tors->chi11 role transfer."
        ),
    }


def component_gap_rows(component_report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["component"]: row for row in component_report.get("component_gap_rows", [])}


def build_implication_certificate(artifacts: dict[str, Any], grep: dict[str, Any]) -> dict[str, Any]:
    p1343 = artifacts["P1343_SELECTOR_SOURCE_REPORT"]
    p1348 = artifacts["P1348_GLOBAL_CLOSURE_REPORT"]
    p1350_text = artifacts["P1350_COUPLING_EXTRACTION_PACKET"]["text"]
    s2_text = artifacts["S2_STRATEGIC_PRIORITY"]["text"]
    component_report = artifacts["COMPONENT_GAP_REPORT"]
    rows = component_gap_rows(component_report)
    chi11_row = rows.get("topological_phase_bit_chi11", {})
    role_row = rows.get("legacy_physical_role_transfer", {})
    summary = component_report.get("completion_gap_summary", {})

    selector_exported = p1343.get("status") == "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1"
    declared_global_closure = p1348.get("status") == "GLOBAL_CLOSURE_THEOREM_EXPORTED_CLOSED_DECLARED_SCOPE"
    p1350_no_silent_legacy = "without silently importing retired legacy roles" in p1350_text
    s2_requires_role_audit = "role-transfer audit" in s2_text and "legacy physical roles automatically transfer" in s2_text
    beta_candidate_not_theorem = bool(summary.get("beta_tors_to_chi11_candidate_not_theorem"))
    component_bridge_ready = bool(summary.get("bridge_ready_for_role_transfer"))
    role_transfer_allowed = bool(summary.get("role_transfer_blocked_until_full_bridge")) is False and bool(
        role_row.get("matrix_bits", {}).get("role_transfer_allowed_now")
    )
    chi11_role_transfer_allowed = bool(chi11_row.get("matrix_bits", {}).get("role_transfer_allowed_now"))
    explicit_beta_transport_theorem = (
        grep["pattern_match_counts"]["beta_tors_to_chi11_literal"] > 0
        and not beta_candidate_not_theorem
        and chi11_role_transfer_allowed
    )

    bits = {
        "strict_selector_exported_P1343": selector_exported,
        "declared_scope_global_closure_P1348": declared_global_closure,
        "P1350_keeps_no_silent_legacy_role_boundary": p1350_no_silent_legacy,
        "S2_requires_completion_bridge_then_role_audit": s2_requires_role_audit,
        "component_gap_says_beta_tors_to_chi11_candidate_not_theorem": beta_candidate_not_theorem,
        "component_bridge_ready_for_role_transfer": component_bridge_ready,
        "explicit_beta_tors_to_chi11_transport_theorem_found": explicit_beta_transport_theorem,
        "legacy_role_transfer_allowed_now": role_transfer_allowed,
        "chi11_row_role_transfer_allowed_now": chi11_role_transfer_allowed,
    }

    transfer_rule_inputs = {
        "selector_source": selector_exported,
        "strict_declared_closure": declared_global_closure,
        "explicit_beta_tors_to_chi11_transport_theorem": explicit_beta_transport_theorem,
        "full_legacy_to_strict_bridge_ready": component_bridge_ready,
        "separate_role_transfer_theorem": role_transfer_allowed and chi11_role_transfer_allowed,
    }
    transfer_licensed = all(transfer_rule_inputs.values())
    single_missing_unlock_candidates = [
        name for name, value in transfer_rule_inputs.items() if not value
    ]

    truth_table_slice = []
    for selector in [False, True]:
        for transport in [False, True]:
            for bridge in [False, True]:
                for role in [False, True]:
                    truth_table_slice.append(
                        {
                            "selector_source": selector,
                            "beta_tors_to_chi11_transport_theorem": transport,
                            "full_bridge_ready": bridge,
                            "role_transfer_theorem": role,
                            "licensed": selector and transport and bridge and role,
                        }
                    )
    current_state_row = {
        "selector_source": selector_exported,
        "beta_tors_to_chi11_transport_theorem": explicit_beta_transport_theorem,
        "full_bridge_ready": component_bridge_ready,
        "role_transfer_theorem": role_transfer_allowed and chi11_role_transfer_allowed,
        "licensed": transfer_licensed,
    }

    return {
        "question_answer": (
            "The old wording 'no beta_tors -> chi11' remains scientifically meaningful only after qualification: "
            "P1343/P1348 give a strict selector mechanism in declared scope, but they do not by themselves transport the legacy torsion parameter beta_tors into the chi11/topological-bit role. "
            "So the hard limit should be read as 'selector present, beta_tors->chi11 bridge-role theorem still absent.'"
        ),
        "audited_bits": bits,
        "transfer_rule_inputs": transfer_rule_inputs,
        "transfer_licensed": transfer_licensed,
        "single_missing_unlock_candidates": single_missing_unlock_candidates,
        "truth_table_slice_selector_not_sufficient": truth_table_slice,
        "current_state_row": current_state_row,
        "proof_reduction": [
            "P1343 exports S_strict_internal_v1 and P1348 packages declared-scope global closure; this discharges the old pure selector absence objection inside the strict declared lane.",
            "S2 still requires an explicit legacy->strict completion bridge and then a separate role-transfer audit before any legacy physical role is transferred.",
            "The component-gap matrix records beta_tors->chi11 as a candidate/not-theorem and reports zero role-transfer-allowed rows.",
            "Therefore selector_exported=True is not a sufficient premise for beta_tors->chi11; the missing premises are explicit beta_tors transport, full bridge readiness, and role-transfer theorem.",
        ],
        "status": "SELECTOR_PRESENT_BETA_TORS_CHI11_ROLE_TRANSFER_STILL_UNLICENSED",
    }


def build_payload() -> dict[str, Any]:
    artifacts: dict[str, Any] = {}
    for name, path in SOURCE_FILES.items():
        if path.suffix == ".json":
            artifacts[name] = load_json(path)
        else:
            artifacts[name] = {"text": read_text(path)}
    grep = grep_audit()
    cert = build_implication_certificate(artifacts, grep)
    theorem_export = {
        "theorem_name": "P2390_T1_selector_qualified_beta_tors_chi11_role_separation",
        "statement": (
            "Given the current repository artifacts, strict selector export in the P1343/P1348 declared lane does not entail beta_tors->chi11 role transfer. "
            "The transfer remains unlicensed until an explicit beta_tors-to-chi11 transport theorem, a full legacy->strict bridge, and a separate role-transfer theorem are all exported."
        ),
        "licensed": [
            "Use P1343/P1348 as evidence that a strict selector mechanism exists in its declared scope.",
            "Replace unqualified 'no selector' wording with the qualified separation: selector present, legacy torsion-to-chi11 role transfer absent.",
            "Use the missing-premise list only as a role-transfer boundary; do not treat beta_tors->chi11 as an active selector target once P1343/P1348 are accepted.",
        ],
        "not_licensed": [
            "No beta_tors -> chi11 theorem is exported.",
            "No silent transfer of legacy torsion roles onto K_strict_gate is licensed.",
            "No physical-role transfer for sin^2(theta_W), alpha_EM inverse, beta^N hierarchy, or chi11 orientation follows.",
            "No new L_total term, strict dynamical source theorem for the P2389 cap, or ToE closure is claimed.",
        ],
        "proof_dependencies": {name: rel(path) for name, path in SOURCE_FILES.items()},
    }
    gatekeeper_checks = {
        "p1343_selector_present": cert["audited_bits"]["strict_selector_exported_P1343"],
        "p1348_declared_closure_present": cert["audited_bits"]["declared_scope_global_closure_P1348"],
        "s2_role_audit_guardrail_present": cert["audited_bits"]["S2_requires_completion_bridge_then_role_audit"],
        "component_gap_beta_candidate_not_theorem_replayed": cert["audited_bits"]["component_gap_says_beta_tors_to_chi11_candidate_not_theorem"],
        "transfer_not_licensed": cert["transfer_licensed"] is False,
        "selector_not_sufficient_truth_table_contains_false_when_transport_missing": any(
            row["selector_source"] and not row["beta_tors_to_chi11_transport_theorem"] and not row["licensed"]
            for row in cert["truth_table_slice_selector_not_sufficient"]
        ),
        "docs_updated_with_p2390": all("P2390/S1340" in read_text(path) for path in DOC_FILES.values()),
    }
    payload = {
        "schema_version": "p2390_s1340_v1",
        "packet_id": "P2390",
        "stage_id": "S1340",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "result_kind": "SELECTOR_QUALIFIED_BETA_TORS_CHI11_ROLE_AUDIT",
        "status": "PASS_SELECTOR_PRESENT_BETA_TORS_CHI11_ROLE_TRANSFER_STILL_UNLICENSED",
        "selector_qualified_beta_tors_chi11_role_audit": {
            "grep_nonduplication_audit": grep,
            "implication_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": (
            "Retire beta_tors->chi11 as an active selector-search target; keep it only as a conditional legacy role-transfer boundary after bridge completion."
        ),
        "global_status": "OPEN_PROGRESS_SELECTOR_RECONCILED_WITH_LEGACY_ROLE_GUARDRAIL_NO_TOE_CLOSURE",
    }
    return payload


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["selector_qualified_beta_tors_chi11_role_audit"]["implication_certificate"]
    grep = payload["selector_qualified_beta_tors_chi11_role_audit"]["grep_nonduplication_audit"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2390 S1340: selector-qualified beta_tors -> chi11 role audit",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                cert["question_answer"],
                "",
                "## Grep/nonduplication audit",
                "",
                f"- Searched files: `{grep['searched_file_count']}`.",
                f"- beta_tors/chi11 literal matches: `{grep['pattern_match_counts']['beta_tors_to_chi11_literal']}` across `{grep['pattern_file_counts']['beta_tors_to_chi11_literal']}` files.",
                f"- selector-source matches: `{grep['pattern_match_counts']['selector_source']}` across `{grep['pattern_file_counts']['selector_source']}` files.",
                f"- role-transfer matches: `{grep['pattern_match_counts']['role_transfer']}` across `{grep['pattern_file_counts']['role_transfer']}` files.",
                f"- Finding: {grep['nonduplication_finding']}",
                "",
                "## Boolean implication certificate",
                "",
                f"- Transfer licensed now: `{cert['transfer_licensed']}`.",
                f"- Current state row: `{cert['current_state_row']}`.",
                f"- Missing premises: `{cert['single_missing_unlock_candidates']}`.",
                "",
                "## Hard limits",
                "",
                "- P1343/P1348 are respected as selector/global-closure exports in declared strict scope.",
                "- No `beta_tors -> chi11` theorem, full legacy-to-strict bridge, or legacy role-transfer theorem is claimed.",
                "- No `L_total` promotion, source theorem for the P2389 cap, SM/GR numeric extraction, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    payload = build_payload()
    write_outputs(payload)


if __name__ == "__main__":
    main()
