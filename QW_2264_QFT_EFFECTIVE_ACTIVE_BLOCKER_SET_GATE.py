#!/usr/bin/env python3
"""QW-2264: QFT effective active blocker-set gate.

Builds conservative effective blocker set:
- declared active blockers from QW-2256,
- dangling refs from QW-2262 locality integrity.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2256 = load("report_qw2256_qft_active_path_blocker_reduction_gate.json")
    p2262 = load("proof_object_qw2262_qft_active_reference_locality_integrity.json")

    active_blockers = set(q2256.get("active_blockers", []))
    dangling_refs = set()
    for row in p2262.get("active_instances_integrity", []):
        for ref in row.get("unresolved_refs", []):
            dangling_refs.add(ref)

    effective_blockers = sorted(active_blockers | dangling_refs)

    flags = {
        "q2256_present": q2256.get("verdict")
        == "QFT_ACTIVE_PATH_BLOCKER_REDUCTION_GATE_PASS_PARTIAL_SINGLE_CORE_BLOCKER",
        "q2262_present": True,
        "declared_active_blockers_present": len(active_blockers) > 0,
        "effective_blocker_set_computed": len(effective_blockers) > 0,
        "effective_set_equals_declared_singleton": effective_blockers
        == ["PositivityToReconstruction_DerivedOrPending"],
        "effective_set_expands_declared_due_locality": len(effective_blockers) > len(active_blockers),
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2256_present"]
        and flags["q2262_present"]
        and flags["declared_active_blockers_present"]
        and flags["effective_blocker_set_computed"]
    )

    verdict = (
        "QFT_EFFECTIVE_ACTIVE_BLOCKER_SET_GATE_PASS_PARTIAL_EXPANDED_BLOCKER_SET"
        if core_ok
        else "QFT_EFFECTIVE_ACTIVE_BLOCKER_SET_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2256": "report_qw2256_qft_active_path_blocker_reduction_gate.json",
            "p2262": "proof_object_qw2262_qft_active_reference_locality_integrity.json",
        },
        "declared_active_blockers": sorted(active_blockers),
        "dangling_refs": sorted(dangling_refs),
        "effective_active_blockers": effective_blockers,
        "scope_boundary": {
            "single_core_blocker_eliminated": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2264_qft_effective_active_blocker_set.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "declared_active_blockers": sorted(active_blockers),
        "n_declared_active_blockers": len(active_blockers),
        "n_dangling_refs": len(dangling_refs),
        "effective_active_blockers": effective_blockers,
        "n_effective_active_blockers": len(effective_blockers),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DISCHARGE_ALL_EFFECTIVE_QFT_BLOCKERS_INCLUDING_LOCALITY_DANGLING_SYMBOLS",
    }

    out_json = ROOT / "report_qw2264_qft_effective_active_blocker_set_gate.json"
    out_md = ROOT / "RAPORT_QW2264_QFT_EFFECTIVE_ACTIVE_BLOCKER_SET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2264: QFT EFFECTIVE ACTIVE BLOCKER-SET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- declared_active_blockers: `{sorted(active_blockers)}`",
                f"- dangling_refs: `{sorted(dangling_refs)}`",
                f"- effective_active_blockers: `{effective_blockers}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "n_declared_active_blockers": len(active_blockers),
                "n_dangling_refs": len(dangling_refs),
                "n_effective_active_blockers": len(effective_blockers),
            }
        )
    )


if __name__ == "__main__":
    main()
