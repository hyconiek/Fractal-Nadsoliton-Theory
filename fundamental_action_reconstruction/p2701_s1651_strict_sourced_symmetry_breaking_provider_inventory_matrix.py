#!/usr/bin/env python3
"""P2701/S1651: strict-sourced symmetry-breaking provider inventory matrix.

After P2700 exhausts Aut(Z12)-invariant selector-functionals, this script does
not replay that lane.  It performs a repository-wide generated-artifact
inventory for a genuinely new strict-sourced, non-premise symmetry-breaking
provider and scores candidates against an explicit acceptance matrix.
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
OUT = GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json"
MD = GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2700": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2697": GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json",
    "H39": GEN / "h39_global_selector_object_absence_audit.json",
    "P739": GEN / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P740": GEN / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_strict_sourced_symmetry_breaking_provider_found",
    "nonpremise_selector_source_exported",
    "qw2191_discharged",
    "strict_selector_closure_exported",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_claimed",
]

POSITIVE_TERMS = ("symmetry", "selector", "orientation", "directed", "sign", "generator", "source", "provider", "qw-2191", "qw2191")
BLOCKER_TERMS = ("premise-based", "premise based", "convention", "gauge", "not strict-derived", "nonexport", "still open", "open", "false")
PASS_TERMS = ("non-premise", "nonpremise", "strict-sourced", "strict sourced", "strict-derived", "source theorem", "qw2191_closed\": true", "qW-2191 discharge")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "p2700_boundary": r"P2700|exhaustive Aut-invariant selector-functional|no-new-live-frontier",
        "provider_search_terms": r"strict-sourced|symmetry-breaking|selector source|orientation source|provider|source theorem",
        "known_premise_boundaries": r"premise-based|convention|gauge|nonexport|QW-2191.*open|pair12 strict-core upgrade",
        "forbidden_promotions": r"strict selector closure|L_total|ToE closure|role transfer|bridge closure",
    }
    return {"tool": "rg", "mode": "P2701 strict-sourced symmetry-breaking provider inventory", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def artifact_text(path: Path) -> str:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
        return json.dumps(obj, sort_keys=True, ensure_ascii=False).lower()
    except Exception:
        return path.read_text(encoding="utf-8", errors="replace").lower()


def generated_inventory() -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    scanned = 0
    for path in sorted(GEN.glob("*.json")):
        if path.name.startswith("p2701_s1651"):
            continue
        scanned += 1
        text = artifact_text(path)
        positive_hits = [term for term in POSITIVE_TERMS if term in text]
        if not positive_hits:
            continue
        blocker_hits = [term for term in BLOCKER_TERMS if term in text]
        pass_hits = [term for term in PASS_TERMS if term.lower() in text]
        explicit_qw_closed = '"qw2191_closed": true' in text or '"qw_2191_closed": true' in text
        explicit_closure_allowed = '"strict_closure_claim_allowed": true' in text or '"strict_core_selector_closure": true' in text
        premise_blocked = any(term in blocker_hits for term in ("premise-based", "premise based", "convention", "gauge", "not strict-derived", "nonexport", "still open"))
        accepted = bool(pass_hits) and explicit_qw_closed and explicit_closure_allowed and not premise_blocked
        candidates.append(
            {
                "path": rel(path),
                "positive_hits": positive_hits[:12],
                "blocker_hits": blocker_hits[:12],
                "pass_hits": pass_hits[:12],
                "explicit_qw_closed": explicit_qw_closed,
                "explicit_closure_allowed": explicit_closure_allowed,
                "premise_or_nonexport_blocked": premise_blocked,
                "accepted_strict_sourced_provider": accepted,
            }
        )
    accepted = [row for row in candidates if row["accepted_strict_sourced_provider"]]
    blocked = [row for row in candidates if not row["accepted_strict_sourced_provider"]]
    return {
        "generated_json_files_scanned": scanned,
        "candidate_files_with_selector_direction_terms": len(candidates),
        "accepted_provider_count": len(accepted),
        "blocked_candidate_count": len(blocked),
        "accepted_providers": accepted,
        "top_blocked_candidates": sorted(blocked, key=lambda row: (len(row["positive_hits"]), -len(row["blocker_hits"])), reverse=True)[:40],
    }


def state_reads() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in INPUTS.items()}
    p2700 = loaded["P2700"]
    p2697 = loaded["P2697"]
    h39 = loaded["H39"]
    p739 = loaded["P739"]
    p740 = loaded["P740"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2700_bounded_no_go": p2700.get("decision", {}).get("bounded_no_go_now") is True,
        "p2700_forbids_replay": "Do not continue Aut(Z12)-invariant selector-functional replay" in p2700.get("decision", {}).get("next_honest_step", ""),
        "p2697_no_new_live_frontier_certificate": p2697.get("decision", {}).get("no_new_live_frontier_certificate") is True,
        "h39_qw2191_still_open": "global_qw_2191_discharge" in h39.get("missing", []) or "QW2191_STILL_OPEN" in h39.get("status", ""),
        "p739_pair12_upgrade_unexported": p739.get("t193_target_exported_on_current_repo_state") is False,
        "p740_pair12_upgrade_unexported": p740.get("t194_target_exported_on_current_repo_state") is False,
    }


def acceptance_matrix(inv: dict[str, Any], reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "acceptance_obligation": "repository exports at least one candidate strict-sourced symmetry-breaking provider",
            "computed_value": inv["accepted_provider_count"],
            "passes": inv["accepted_provider_count"] > 0,
        },
        {
            "acceptance_obligation": "candidate is non-premise/non-convention and explicitly closes QW-2191",
            "computed_value": {"accepted_provider_count": inv["accepted_provider_count"], "h39_qw2191_still_open": reads["h39_qw2191_still_open"]},
            "passes": False,
        },
        {
            "acceptance_obligation": "candidate upgrades pair12 strict core without P739/P740 nonexport boundary",
            "computed_value": {"p739_unexported": reads["p739_pair12_upgrade_unexported"], "p740_unexported": reads["p740_pair12_upgrade_unexported"]},
            "passes": False,
        },
        {
            "acceptance_obligation": "candidate reopens no-new-live-frontier after P2697/P2700",
            "computed_value": {"p2697_no_new_live_frontier": reads["p2697_no_new_live_frontier_certificate"], "p2700_bounded_no_go": reads["p2700_bounded_no_go"]},
            "passes": False,
        },
    ]


def decision(rows: list[dict[str, Any]], inv: dict[str, Any]) -> dict[str, Any]:
    no_provider = all(not row["passes"] for row in rows)
    return {
        "decision": "P2701_STRICT_SOURCED_SYMMETRY_BREAKING_PROVIDER_INVENTORY_NO_PROVIDER_NO_FALSE_PASS",
        "bounded_no_go_now": no_provider,
        "reason": f"Inventory scanned {inv['generated_json_files_scanned']} generated JSON artifacts and found {inv['candidate_files_with_selector_direction_terms']} selector/direction/source candidates, but 0 accepted strict-sourced non-premise symmetry-breaking providers.",
        "next_honest_step": "No current artifact exports the new provider required after P2700.  The next admissible move must construct a genuinely new strict-sourced symmetry-breaking object with explicit non-premise provenance and QW-2191/pair12 upgrade obligations, or pivot to a different newly exported typed object outside closed lanes.  Otherwise preserve the P2697-P2701 no-new-live-frontier certificate.",
        "forbidden_promotions": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    inv = payload["provider_inventory"]
    lines = ["# P2701/S1651 strict-sourced symmetry-breaking provider inventory matrix", "", f"Status: `{payload['status']}`", "", "## Inventory"]
    lines.append(f"- generated JSON files scanned: `{inv['generated_json_files_scanned']}`")
    lines.append(f"- candidate files with selector/direction/source terms: `{inv['candidate_files_with_selector_direction_terms']}`")
    lines.append(f"- accepted strict-sourced providers: `{inv['accepted_provider_count']}`")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    inv = generated_inventory()
    rows = acceptance_matrix(inv, reads)
    payload: dict[str, Any] = {
        "status": "P2701_STRICT_SOURCED_SYMMETRY_BREAKING_PROVIDER_INVENTORY_MATRIX_NO_PROVIDER",
        "content_grep": content_grep(),
        "state_reads": reads,
        "provider_inventory": inv,
        "acceptance_matrix": rows,
        "decision": decision(rows, inv),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2701/S1651 strict-sourced symmetry-breaking provider inventory",
        "## P2701/S1651 strict-sourced symmetry-breaking provider inventory\n\n"
        "`P2701/S1651` performs the post-P2700 provider inventory instead of replaying Aut-invariant selector-functionals.  It scans generated JSON artifacts for selector/direction/source candidates and applies a strict acceptance matrix requiring non-premise provenance plus explicit `QW-2191`/pair12 upgrade authority.  No accepted strict-sourced symmetry-breaking provider is found; no strict selector closure, role-bearing `L_total`, bridge closure, role transfer, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2701/S1651 provider inventory Ltotal guard",
        "## P2701/S1651 provider inventory Ltotal guard\n\n"
        "`P2701/S1651` is an inventory/acceptance-matrix no-provider result, not a variational source.  It keeps `L_total`, strict selector closure, role transfer, bridge closure, and ToE closure unpromoted.\n",
    )
    append_once(
        AGENTS,
        "Current strict-sourced symmetry-breaking provider inventory guardrail (P2701/S1651, 2026-06-13)",
        "## Current strict-sourced symmetry-breaking provider inventory guardrail (P2701/S1651, 2026-06-13)\n\n"
        "- P2701 performs a generated-artifact inventory for the new strict-sourced symmetry-breaking provider required after P2700.  Existing selector/direction/source candidates remain blocked by premise/convention/nonexport/open-`QW-2191` boundaries; no accepted provider is exported.\n"
        "- Do not claim strict selector closure, pair12 strict-core upgrade, role-bearing `L_total`, bridge closure, role transfer, or ToE closure from inventory hits.\n"
        "- A next admissible move must actually construct a new strict-sourced symmetry-breaking object/provider or introduce a different new typed object outside closed lanes; otherwise preserve the P2697-P2701 no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
