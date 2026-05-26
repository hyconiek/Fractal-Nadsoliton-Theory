#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2168 = GEN / "p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json"
IN_2172 = GEN / "p2172_s1122_strict_qw2191_obligation_validator_o1_o3_o4_update.json"
OUT = GEN / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.json"
MD = GEN / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def _index_evals(evals: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {e.get("id", f"missing_id_{i}"): e for i, e in enumerate(evals)}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2168 = load(IN_2168)
    p2172 = load(IN_2172)

    base_evals = (
        (p2168.get("strict_qw2191_theorem_obligations_executable_validator", {}) or {}).get("obligation_evaluations", [])
        or []
    )
    upd_evals = (
        (p2172.get("strict_qw2191_obligation_validator_o1_o3_o4_update", {}) or {}).get("updated_obligation_evaluations", [])
        or []
    )

    bidx = _index_evals(base_evals)
    uidx = _index_evals(upd_evals)

    obligation_ids = sorted(set(bidx) | set(uidx))
    reconciled = []
    for oid in obligation_ids:
        b = bidx.get(oid, {})
        u = uidx.get(oid, {})
        promoted = (not bool(b.get("pass", False))) and bool(u.get("pass", False))
        consistent = True
        if u:
            # basic consistency rule: pass=True must imply PASS-like status tag
            status = str(u.get("status", ""))
            if bool(u.get("pass", False)) and "PASS" not in status:
                consistent = False
        reconciled.append(
            {
                "id": oid,
                "base_pass": bool(b.get("pass", False)),
                "updated_pass": bool(u.get("pass", False)),
                "promoted": promoted,
                "consistent": consistent,
            }
        )

    n_total = len(reconciled)
    n_promoted = sum(1 for r in reconciled if r["promoted"])
    n_inconsistent = sum(1 for r in reconciled if not r["consistent"])
    all_consistent = n_inconsistent == 0
    all_pass = all(r["updated_pass"] for r in reconciled) if reconciled else False

    result_kind = (
        "PASS_STRICT_QW2191_OBLIGATION_INDEPENDENT_CROSSCHECK_WITH_TRACE"
        if all_consistent
        else "OPEN_STRICT_QW2191_OBLIGATION_INDEPENDENT_CROSSCHECK_INCONSISTENT"
    )

    payload = {
        "schema_version": "p2173_s1123_v1",
        "packet_id": "P2173",
        "stage_id": "S1123",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_obligation_independent_crosschecker": {
            "source_base_validator": str(IN_2168.relative_to(ROOT)),
            "source_updated_validator": str(IN_2172.relative_to(ROOT)),
            "reconciled_obligations": reconciled,
            "summary": {
                "n_total": n_total,
                "n_promoted": n_promoted,
                "n_inconsistent": n_inconsistent,
                "all_consistent": all_consistent,
                "all_pass": all_pass,
            },
            "scope_limit": "cross-checks obligation-evaluator consistency only; no selector-closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2174_candidate",
            "goal": "export evidentiary artifacts for any still-open obligations and iterate cross-check until all required pass",
        },
        "gatekeeper_checks": {
            "crosscheck_exported": True,
            "all_consistent": all_consistent,
            "all_required_pass": bool((p2172.get("gatekeeper_checks", {}) or {}).get("all_required_pass", False)),
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2172.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2172.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2173 S1123: strict QW-2191 obligation independent cross-checker",
                "",
                f"- Result kind: `{result_kind}`",
                f"- all_consistent: `{all_consistent}`",
                f"- n_promoted: `{n_promoted}`",
                f"- all_pass: `{all_pass}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
