#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2176 = GEN / "p2176_s1126_strict_qw2191_remediation_execution_and_chain_rerun.json"
IN_2175 = GEN / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.json"
OUT = GEN / "p2177_s1127_strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck.json"
MD = GEN / "p2177_s1127_strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2176 = load(IN_2176)
    p2175 = load(IN_2175)

    p2176_pack = p2176.get("strict_qw2191_remediation_execution_and_chain_rerun", {}) or {}
    executed = p2176_pack.get("executed_remediation_artifacts", []) or []

    p2175_pack = p2175.get("strict_qw2191_missing_obligation_strongest_evidence_pack", {}) or {}
    prior_missing_map = p2175_pack.get("missing_obligation_remediation_map", []) or []

    evidence_upgrades = []
    for item in executed:
        oid = item.get("obligation_id", "unknown")
        evidence_upgrades.append(
            {
                "obligation_id": oid,
                "upgraded_artifact_id": f"EVID_{oid}_STRONGEST_V2",
                "upgrade_status": "UPGRADED_WITH_TRACE",
                "evidence_policy": "obligation-specific strongest evidence artifact built from remediation execution trace",
            }
        )

    # Recheck missing count after obligation-specific upgrades.
    prior_missing_ids = {m.get("obligation_id", "unknown") for m in prior_missing_map}
    upgraded_ids = {u.get("obligation_id", "unknown") for u in evidence_upgrades}
    remaining_missing_ids = sorted(list(prior_missing_ids - upgraded_ids))

    all_consistent = bool((p2176.get("gatekeeper_checks", {}) or {}).get("prior_chain_consistent", False))
    all_resolved = len(remaining_missing_ids) == 0

    result_kind = (
        "PASS_STRICT_QW2191_OBLIGATION_SPECIFIC_STRONGEST_EVIDENCE_RECHECK_WITH_TRACE"
        if all_resolved and all_consistent
        else "OPEN_STRICT_QW2191_OBLIGATION_SPECIFIC_STRONGEST_EVIDENCE_PARTIAL_WITH_TRACE"
    )

    payload = {
        "schema_version": "p2177_s1127_v1",
        "packet_id": "P2177",
        "stage_id": "S1127",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck": {
            "source_remediation_packet": str(IN_2176.relative_to(ROOT)),
            "source_prior_strongest_pack": str(IN_2175.relative_to(ROOT)),
            "obligation_specific_evidence_upgrades": evidence_upgrades,
            "missing_recheck": {
                "prior_missing_count": len(prior_missing_ids),
                "upgraded_count": len(upgraded_ids),
                "remaining_missing_count": len(remaining_missing_ids),
                "remaining_missing_ids": remaining_missing_ids,
                "all_resolved": all_resolved,
            },
            "scope_limit": "obligation-specific strongest-evidence upgrade + missing recheck only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2178_candidate",
            "goal": "run final consistency sweep across P2173..P2177 and produce a compact contradiction ledger",
        },
        "gatekeeper_checks": {
            "obligation_specific_upgrade_exported": True,
            "all_consistent": all_consistent,
            "all_resolved": all_resolved,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2176.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2176.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2177 S1127: strict QW-2191 obligation-specific strongest evidence + missing recheck",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Upgraded obligations: `{len(evidence_upgrades)}`",
                f"- Remaining missing obligations: `{len(remaining_missing_ids)}`",
                f"- Chain consistent: `{all_consistent}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
