#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2172 = GEN / "p2172_s1122_strict_qw2191_obligation_validator_o1_o3_o4_update.json"
IN_2173 = GEN / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.json"
OUT = GEN / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.json"
MD = GEN / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def classify_obligation(e: dict[str, Any]) -> dict[str, Any]:
    obligation_id = str(e.get("id", "unknown"))
    is_pass = bool(e.get("pass", False))
    status = str(e.get("status", ""))

    # Explicit, machine-readable split: what is evidence vs what is premise.
    if obligation_id.startswith("O2"):
        evidence = ["P2169 dedicated witness object", "P2170 validator promotion"]
        premise = ["selected strict branch exists"]
    elif obligation_id.startswith("O5"):
        evidence = ["P2169 transition rule witness", "P2170 validator promotion"]
        premise = ["admission policy remains no-shortcut"]
    elif obligation_id.startswith("O1"):
        evidence = ["P2171 branch identifiability witness", "P2172 validator promotion"]
        premise = ["instantiated branch remains admissible"]
    elif obligation_id.startswith("O3"):
        evidence = ["P2171 noncyclic evidence bundle", "P2172 validator promotion"]
        premise = ["noncyclic anchor provider remains valid"]
    elif obligation_id.startswith("O4"):
        evidence = ["P2171 guardrail certificate", "P2172 validator promotion"]
        premise = ["no legacy-role transfer claim policy"]
    else:
        evidence = ["unclassified"]
        premise = ["unclassified"]

    strength = "strong_evidence" if is_pass and "PASS" in status else "missing_or_partial"
    return {
        "id": obligation_id,
        "status": status,
        "pass": is_pass,
        "evidence_items": evidence,
        "premise_items": premise,
        "evidence_strength": strength,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2172 = load(IN_2172)
    p2173 = load(IN_2173)

    evals = (
        (p2172.get("strict_qw2191_obligation_validator_o1_o3_o4_update", {}) or {}).get("updated_obligation_evaluations", [])
        or []
    )
    normalized = [classify_obligation(e) for e in evals]

    missing_or_partial = [n for n in normalized if n["evidence_strength"] != "strong_evidence"]

    # Release 8.1 historical closure scope note (non-authoritative in current lane).
    release_81_scope_note = {
        "reported_history": "Release 8.1 reports QW-2191 closure within declared Release-8 scope.",
        "scope_limit": "current strict CMP2/QW-2191 lane treats that as historical scoped claim, not automatic global discharge for this lane",
        "operational_policy": "require executable per-obligation evidence in current lane regardless of historical scoped closure claims",
    }

    result_kind = (
        "PASS_STRICT_QW2191_EVIDENTIARY_TRACE_NORMALIZATION_WITH_TRACE"
        if not missing_or_partial
        else "OPEN_STRICT_QW2191_MISSING_OBLIGATION_EVIDENCE_WITH_TRACE"
    )

    payload = {
        "schema_version": "p2174_s1124_v1",
        "packet_id": "P2174",
        "stage_id": "S1124",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack": {
            "source_validator_update": str(IN_2172.relative_to(ROOT)),
            "source_independent_crosscheck": str(IN_2173.relative_to(ROOT)),
            "normalized_obligation_traces": normalized,
            "missing_or_partial_obligations": missing_or_partial,
            "release_81_scope_note": release_81_scope_note,
            "summary": {
                "n_total": len(normalized),
                "n_missing_or_partial": len(missing_or_partial),
            },
            "scope_limit": "trace normalization + missing-obligation witness pack only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2175_candidate",
            "goal": "export dedicated strongest-evidence artifacts for any still-missing obligations and rerun P2173/P2174 consistency chain",
        },
        "gatekeeper_checks": {
            "trace_normalization_exported": True,
            "missing_obligations_identified": len(missing_or_partial) > 0,
            "all_consistent": bool((p2173.get("gatekeeper_checks", {}) or {}).get("all_consistent", False)),
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2173.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2173.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2174 S1124: strict QW-2191 evidentiary-trace normalizer + missing-obligation witness pack",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Missing/partial obligations: `{len(missing_or_partial)}`",
                "- Includes explicit evidence vs premise split per obligation.",
                "- Includes Release 8.1 scoped-closure note (non-automatic for current lane).",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
