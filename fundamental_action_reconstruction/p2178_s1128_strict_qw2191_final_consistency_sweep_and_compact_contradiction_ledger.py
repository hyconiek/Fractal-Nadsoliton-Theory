#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
INPUTS = [
    ("P2173", GEN / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.json"),
    ("P2174", GEN / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.json"),
    ("P2175", GEN / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.json"),
    ("P2176", GEN / "p2176_s1126_strict_qw2191_remediation_execution_and_chain_rerun.json"),
    ("P2177", GEN / "p2177_s1127_strict_qw2191_obligation_specific_strongest_evidence_and_missing_recheck.json"),
]
OUT = GEN / "p2178_s1128_strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger.json"
MD = GEN / "p2178_s1128_strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "gatekeeper_checks": {}, "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    docs = {pid: load(path) for pid, path in INPUTS}

    ledger = []
    contradictions = []
    for pid, path in INPUTS:
        doc = docs[pid]
        gk = doc.get("gatekeeper_checks", {}) or {}
        no_sel = bool(gk.get("no_selector_closure_claimed", False))
        no_toe = bool(gk.get("no_toe_closure_claimed", False))
        if not no_sel:
            contradictions.append({"packet": pid, "kind": "selector_closure_claim_violation"})
        if not no_toe:
            contradictions.append({"packet": pid, "kind": "toe_closure_claim_violation"})
        ledger.append(
            {
                "packet": pid,
                "path": str(path.relative_to(ROOT)),
                "result_kind": doc.get("result_kind", "UNKNOWN"),
                "no_selector_closure_claimed": no_sel,
                "no_toe_closure_claimed": no_toe,
                "full_d3_covariance_transport_proven": bool(gk.get("full_d3_covariance_transport_proven", False)),
                "c3_theorem_proven": bool(gk.get("c3_theorem_proven", False)),
            }
        )

    # Additional compact contradiction checks for chain summary compatibility.
    p2176_all_resolved = bool((docs["P2176"].get("gatekeeper_checks", {}) or {}).get("all_resolved", False))
    p2177_all_resolved = bool((docs["P2177"].get("gatekeeper_checks", {}) or {}).get("all_resolved", False))
    if p2177_all_resolved and not p2176_all_resolved:
        contradictions.append(
            {
                "packet": "P2177_vs_P2176",
                "kind": "resolution_inversion",
                "detail": "P2177 says resolved while P2176 does not",
            }
        )

    contradiction_free = len(contradictions) == 0
    result_kind = (
        "PASS_STRICT_QW2191_FINAL_CONSISTENCY_SWEEP_AND_LEDGER_WITH_TRACE"
        if contradiction_free
        else "OPEN_STRICT_QW2191_FINAL_CONSISTENCY_SWEEP_CONTRADICTIONS_FOUND"
    )

    payload = {
        "schema_version": "p2178_s1128_v1",
        "packet_id": "P2178",
        "stage_id": "S1128",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger": {
            "input_packets": [str(path.relative_to(ROOT)) for _, path in INPUTS],
            "compact_ledger": ledger,
            "contradictions": contradictions,
            "summary": {
                "n_packets": len(INPUTS),
                "n_contradictions": len(contradictions),
                "contradiction_free": contradiction_free,
            },
            "scope_limit": "consistency sweep and compact contradiction ledger only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2179_candidate",
            "goal": "if contradictions remain, repair exactly those packets; otherwise freeze compact ledger snapshot for audit handoff",
        },
        "gatekeeper_checks": {
            "final_sweep_exported": True,
            "contradiction_free": contradiction_free,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((docs["P2177"].get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((docs["P2177"].get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2178 S1128: strict QW-2191 final consistency sweep + compact contradiction ledger",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Packets checked: `{len(INPUTS)}`",
                f"- Contradictions found: `{len(contradictions)}`",
                f"- Contradiction-free: `{contradiction_free}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
