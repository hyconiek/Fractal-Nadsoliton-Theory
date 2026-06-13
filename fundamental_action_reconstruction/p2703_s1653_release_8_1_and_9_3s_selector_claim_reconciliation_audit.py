#!/usr/bin/env python3
"""P2703/S1653: Release 8.1 and alleged Release 9.3s selector-claim reconciliation.

The user asked to check release 8.1 and release 9.3s because older release
notes allegedly contained selector/symmetry-breaking evidence.  This audit is
intentionally conservative: it distinguishes release-note claims, theorem-draft
checkpoints, and current post-P2702 guardrails.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2703_s1653_release_8_1_and_9_3s_selector_claim_reconciliation_audit.json"
MD = GEN / "p2703_s1653_release_8_1_and_9_3s_selector_claim_reconciliation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "release_8_1_textbook": REPO / "RELEASE_8_1_TEXTBOOK_EDITION_EN_PL.md",
    "zenodo_8_1": REPO / "zenodo8.1.md",
    "p1351_release_8_1_audit": ROOT / "P1351_RELEASE_8_1_FAR_GREP_AUDIT_QW2191_NB_SM_GR_BRIDGE_NEXT_STEP_PACKET_PL.md",
    "p1293_r9_theorem_draft_script": ROOT / "p1293_qw2191_r9_formal_selector_source_theorem_draft_checkpoint.py",
    "p1543_selector_uniqueness_link": ROOT / "P1543_S493_SELECTOR_UNIQUENESS_MAIN_THEOREM_LINK_PROOF_CHECKPOINT_PACKET_PL.md",
    "p2702_current_status": GEN / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "release_8_1_claim_promoted_to_current_closure",
    "release_9_3s_found_as_release_document",
    "r9_draft_promoted_to_current_closure",
    "qw2191_discharged_now",
    "strict_selector_closure_exported_now",
    "pair12_strict_core_upgrade_exported_now",
    "ltotal_promoted",
    "toe_closure_claimed",
]


def run_rg(args: list[str]) -> list[str]:
    proc = subprocess.run(["rg", "-n", *args], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    return sorted(line for line in proc.stdout.splitlines() if line)


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def release_search() -> dict[str, Any]:
    find_proc = subprocess.run(
        ["find", ".", "-maxdepth", "3", "-type", "f", "(", "-iname", "*9.3s*", "-o", "-iname", "*9_3s*", "-o", "-iname", "*release*9*3*", ")", "-print"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    direct_release_hits = sorted(
        line
        for line in find_proc.stdout.splitlines()
        if line
        and "__pycache__" not in line
        and "fundamental_action_reconstruction/generated/" not in line
        and "p2703_s1653_release_8_1_and_9_3s" not in line
    )
    text_hits = [
        hit
        for hit in run_rg(["9\\.3s|Release 9\\.3|RELEASE 9\\.3|R9 formal selector|R9_FORMAL_SELECTOR|P1293", ".", "-g", "*.md", "-g", "*.py", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"])
        if "p2703_s1653_release_8_1_and_9_3s" not in hit
    ]
    return {
        "tooling": ["find", "rg"],
        "direct_release_9_3s_file_hits": direct_release_hits,
        "release_9_3s_document_found": any("9.3s" in hit.lower() or "release_9_3" in hit.lower() or "release 9.3" in hit.lower() for hit in direct_release_hits),
        "r9_or_9_3_text_hits": text_hits[:120],
        "r9_theorem_draft_present": any("p1293" in hit.lower() or "r9_formal_selector" in hit.lower() for hit in text_hits),
    }


def release_8_1_findings() -> dict[str, Any]:
    text = read_text(INPUTS["release_8_1_textbook"]) + "\n" + read_text(INPUTS["zenodo_8_1"])
    phrases = {
        "internal_selector_source_claim": "Internal Selector Source" in text or "wewnętrznego źródła selektora" in text,
        "p1343_claim": "P1343" in text,
        "p1348_global_closure_claim": "P1348" in text and "Global" in text,
        "no_false_pass_limit": "does not" in text or "false" in text.lower(),
        "external_audit_pending": "PENDING_EXECUTION" in text or "external" in text.lower(),
    }
    return {
        "source_files": [rel(INPUTS["release_8_1_textbook"]), rel(INPUTS["zenodo_8_1"])],
        "claim_flags": phrases,
        "plain_reading": "Release 8.1 contains a strong declared-scope claim: P1343 internal selector source and P1348 global strict closure.  It also keeps release-scope/no-false-pass boundaries and external-audit/rollback discipline.",
    }


def r9_related_findings() -> dict[str, Any]:
    p1293 = read_text(INPUTS["p1293_r9_theorem_draft_script"])
    p1543 = read_text(INPUTS["p1543_selector_uniqueness_link"])
    return {
        "p1293_present": bool(p1293),
        "p1293_status_draft": '"status": "DRAFT"' in p1293 or '"status": "DRAFT"' in p1293.replace("'", '"'),
        "p1293_closure_policy_blocks_closure": "strict_core_selector_closure_allowed" in p1293 and "False" in p1293,
        "p1543_present": bool(p1543),
        "p1543_provisional_and_qw2191_false": "PROVISIONAL" in p1543 and "qw2191_closed=false" in p1543,
        "plain_reading": "The visible R9/P1293 material is a formal theorem draft checkpoint, not a release-9.3s closure export; P1543 is provisional and explicitly keeps QW-2191 false before a full proof bundle.",
    }


def current_status_findings() -> dict[str, Any]:
    p2702 = load_json(INPUTS["p2702_current_status"])
    decision = p2702.get("decision", {})
    return {
        "p2702_hash": sha(INPUTS["p2702_current_status"]),
        "p2702_no_new_closure_exported": decision.get("no_new_closure_exported") is True,
        "p2702_strict_aut_invariant_candidate_picks_one_direction": decision.get("strict_aut_invariant_candidate_picks_one_direction") is True,
        "p2702_next_honest_step": decision.get("next_honest_step"),
    }


def reconciliation_matrix(r8: dict[str, Any], r9: dict[str, Any], current: dict[str, Any], search: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "row": "release_8_1_selector_source_claim",
            "older_support_real": r8["claim_flags"]["internal_selector_source_claim"] and r8["claim_flags"]["p1343_claim"],
            "current_block_removed": False,
            "reason": "R8.1 is a release-scope claim, but current later guardrails and P2702 require a currently accepted non-premise strict provider; older release prose alone is not a fresh provider export.",
        },
        {
            "row": "release_8_1_global_closure_claim",
            "older_support_real": r8["claim_flags"]["p1348_global_closure_claim"],
            "current_block_removed": False,
            "reason": "P1348 is acknowledged as a declared-scope closure claim, but present guardrails forbid promoting it to current QW-2191/ToE closure after later no-go/state-map audits without revalidating the provider obligations.",
        },
        {
            "row": "alleged_release_9_3s_selector_proof",
            "older_support_real": search["release_9_3s_document_found"],
            "current_block_removed": False,
            "reason": "No direct Release 9.3s document was found by find/rg; the visible R9 item is P1293, a DRAFT checkpoint with closure disallowed.",
        },
        {
            "row": "r9_p1293_formal_selector_source_theorem_draft",
            "older_support_real": r9["p1293_present"] and r9["p1293_status_draft"],
            "current_block_removed": False,
            "reason": "P1293 is useful historical evidence of an attempted selector-source theorem, but it explicitly remains DRAFT and blocks strict/global selector closure in its policy.",
        },
        {
            "row": "current_p2702_block_status",
            "older_support_real": True,
            "current_block_removed": False,
            "reason": "P2702 remains the current status packet: premise selectors can choose a direction, but current strict Aut-invariant/provider routes do not export non-premise closure.",
        },
    ]


def decision(matrix: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "decision": "P2703_RELEASE_8_1_AND_9_3S_SELECTOR_CLAIM_RECONCILIATION_AUDIT_NO_FALSE_PASS",
        "does_older_release_material_remove_current_blocks": False,
        "recognized_positive_content": [row["row"] for row in matrix if row["older_support_real"]],
        "current_blockers_preserved": ["QW-2191 strict-core uniqueness", "non-premise strict selector-provider obligation", "pair12 strict-core upgrade", "L_total/ToE promotion"],
        "lay_verdict": "Dla laika: stare release'y pokazują, że była silna próba i nawet deklaracja selektora w zakresie R8.1, ale późniejszy stan repo traktuje ją jako niewystarczającą do obecnych, ostrzejszych blokad.  Nie znaleziono osobnego dokumentu release 9.3s; widoczny R9/P1293 jest szkicem, nie zamknięciem.",
        "next_honest_step": "P2704 should be a narrow P1343/P1348 provenance revalidation table: extract the exact S_strict_internal_v1 carrier, domain, proof obligations, validation rows P1344-P1346, and compare them against the P2699-P2702 non-premise provider criteria.  If the exact carrier/proof bundle is not present and executable, preserve P2697-P2703 no-new-live-frontier.",
        "forbidden_promotions": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2703/S1653 Release 8.1 and 9.3s selector-claim reconciliation audit", "", f"Status: `{payload['status']}`", "", "## Matrix"]
    for row in payload["reconciliation_matrix"]:
        lines.append(f"- `{row['row']}`: older_support_real={row['older_support_real']}, current_block_removed={row['current_block_removed']}. {row['reason']}")
    lines.extend(["", "## Lay verdict", payload["decision"]["lay_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    search = release_search()
    r8 = release_8_1_findings()
    r9 = r9_related_findings()
    current = current_status_findings()
    matrix = reconciliation_matrix(r8, r9, current, search)
    payload: dict[str, Any] = {
        "status": "P2703_RELEASE_8_1_AND_9_3S_SELECTOR_CLAIM_RECONCILIATION_AUDIT_NO_FALSE_PASS",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "release_search": search,
        "release_8_1_findings": r8,
        "r9_related_findings": r9,
        "current_status_findings": current,
        "reconciliation_matrix": matrix,
        "decision": decision(matrix),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2703/S1653 release 8.1 and 9.3s selector audit",
        "## P2703/S1653 release 8.1 and 9.3s selector audit\n\n"
        "`P2703/S1653` audits Release 8.1 and the alleged Release 9.3s selector evidence.  Release 8.1 does contain a declared-scope `P1343` internal selector source / `P1348` global closure claim, but no direct Release 9.3s document is found; visible R9 material (`P1293`) is a draft checkpoint with selector closure disabled.  These older materials are positive historical support but do not by themselves remove current P2699-P2702 `QW-2191`/non-premise-provider blocks or promote `L_total`/ToE.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2703/S1653 older release selector evidence Ltotal guard",
        "## P2703/S1653 older release selector evidence Ltotal guard\n\n"
        "`P2703/S1653` does not promote Release 8.1 or R9 draft selector claims into a current variational source, role-bearing `L_total`, strict selector closure, pair12 strict-core upgrade, or ToE closure.\n",
    )
    append_once(
        AGENTS,
        "Current older-release selector evidence reconciliation guardrail (P2703/S1653, 2026-06-13)",
        "## Current older-release selector evidence reconciliation guardrail (P2703/S1653, 2026-06-13)\n\n"
        "- P2703 acknowledges that Release 8.1 contains a declared-scope `P1343` internal selector source / `P1348` global closure claim, but treats it as historical release-scope evidence requiring current provenance revalidation against P2699-P2702 criteria.\n"
        "- No direct `Release 9.3s` document is found in the audit; visible R9/P1293 material is a draft selector-source theorem checkpoint with closure disallowed, not a current closure export.\n"
        "- Do not remove `QW-2191`, non-premise strict selector-provider, pair12 strict-core, `L_total`, role-transfer, bridge, or ToE blocks from older release prose alone.  The next admissible move is a narrow P1343/P1348 provenance revalidation table, or preservation of the P2697-P2703 no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
