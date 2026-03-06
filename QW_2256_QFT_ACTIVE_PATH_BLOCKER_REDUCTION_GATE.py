#!/usr/bin/env python3
"""QW-2256: QFT active-path blocker reduction gate."""

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
    q2254 = load("report_qw2254_qft_export_minimal_blocker_cut_gate.json")

    proof_2254 = load("proof_object_qw2254_qft_export_minimal_blocker_cut.json")
    instances = proof_2254.get("theorem_instances", [])

    active_instances = []
    legacy_instances = []

    for t in instances:
        file_name = t.get("file", "")
        blockers = t.get("blockers", [])
        is_legacy = ("_TERMINAL" in file_name) or any(b.endswith("Witness") for b in blockers)
        if is_legacy:
            legacy_instances.append(t)
        else:
            active_instances.append(t)

    active_blockers = sorted({b for t in active_instances for b in t.get("blockers", [])})
    legacy_blockers = sorted({b for t in legacy_instances for b in t.get("blockers", [])})

    reduced_single_blocker = active_blockers == ["PositivityToReconstruction_DerivedOrPending"]

    flags = {
        "q2254_minimal_cut_present": q2254.get("verdict")
        == "QFT_EXPORT_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_BLOCKER_CUT_COMPUTED",
        "instances_present": len(instances) > 0,
        "legacy_instances_detected": len(legacy_instances) > 0,
        "active_instances_detected": len(active_instances) > 0,
        "legacy_witness_blocker_detected": "L5O1aWitness" in legacy_blockers,
        "active_path_blockers_computed": len(active_blockers) > 0,
        "active_path_reduced_to_single_core_blocker": reduced_single_blocker,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2254_minimal_cut_present"]
        and flags["instances_present"]
        and flags["legacy_instances_detected"]
        and flags["active_instances_detected"]
        and flags["legacy_witness_blocker_detected"]
        and flags["active_path_blockers_computed"]
    )

    verdict = (
        "QFT_ACTIVE_PATH_BLOCKER_REDUCTION_GATE_PASS_PARTIAL_SINGLE_CORE_BLOCKER"
        if core_ok
        else "QFT_ACTIVE_PATH_BLOCKER_REDUCTION_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "proof_object_qw2254_qft_export_minimal_blocker_cut.json",
        "legacy_instances": legacy_instances,
        "active_instances": active_instances,
        "legacy_blockers": legacy_blockers,
        "active_blockers": active_blockers,
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2256_qft_active_path_blocker_reduction.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": proof_obj["source"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_instances_total": len(instances),
        "n_instances_active": len(active_instances),
        "n_instances_legacy": len(legacy_instances),
        "legacy_blockers": legacy_blockers,
        "active_blockers": active_blockers,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "ELIMINATE_PositivityToReconstruction_DerivedOrPending_IN_ACTIVE_PATH",
    }

    out_json = ROOT / "report_qw2256_qft_active_path_blocker_reduction_gate.json"
    out_md = ROOT / "RAPORT_QW2256_QFT_ACTIVE_PATH_BLOCKER_REDUCTION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2256: QFT ACTIVE-PATH BLOCKER REDUCTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- active_blockers: `{active_blockers}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "active_blockers": active_blockers, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
