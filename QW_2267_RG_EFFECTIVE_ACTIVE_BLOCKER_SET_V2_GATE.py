#!/usr/bin/env python3
"""QW-2267: RG effective active blocker-set v2 gate.

Contracts effective blocker set by removing refs that are resolved by explicit
axiomatic bridge availability (QW-2265), while keeping non-axiomatic blockers.
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
    p2263 = load("proof_object_qw2263_rg_effective_active_blocker_set.json")
    p2265 = load("proof_object_qw2265_rg_canonical_export_bridge_availability.json")

    effective_v1 = set(p2263.get("effective_active_blockers", []))
    bridged_refs = set(p2265.get("bridged_unresolved_refs", []))
    residual = sorted(effective_v1 - bridged_refs)

    flags = {
        "q2263_effective_set_present": len(effective_v1) > 0,
        "q2265_bridge_available": "RG_CanonicalAction_to_WellPosedness_EXPORT" in set(p2265.get("bridge_symbols", [])),
        "bridged_refs_identified": len(bridged_refs) > 0,
        "residual_set_computed": len(residual) > 0,
        "residual_single_core_blocker": residual == ["RGGlobalWellPosednessAllScales_DerivedOrPending"],
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2263_effective_set_present"]
        and flags["q2265_bridge_available"]
        and flags["residual_set_computed"]
    )

    verdict = (
        "RG_EFFECTIVE_ACTIVE_BLOCKER_SET_V2_GATE_PASS_PARTIAL_SINGLE_NON_AXIOMATIC_CORE_BLOCKER"
        if core_ok
        else "RG_EFFECTIVE_ACTIVE_BLOCKER_SET_V2_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "p2263": "proof_object_qw2263_rg_effective_active_blocker_set.json",
            "p2265": "proof_object_qw2265_rg_canonical_export_bridge_availability.json",
        },
        "effective_active_blockers_v1": sorted(effective_v1),
        "bridged_refs_from_q2265": sorted(bridged_refs),
        "effective_active_blockers_v2_residual": residual,
        "scope_boundary": {
            "bridge_reduction_is_axiomatic_layer_only": True,
            "non_axiomatic_core_blocker_remaining": True,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2267_rg_effective_active_blocker_set_v2.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "effective_active_blockers_v1": sorted(effective_v1),
        "n_effective_active_blockers_v1": len(effective_v1),
        "bridged_refs": sorted(bridged_refs),
        "n_bridged_refs": len(bridged_refs),
        "effective_active_blockers_v2_residual": residual,
        "n_effective_active_blockers_v2_residual": len(residual),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DISCHARGE_RESIDUAL_NON_AXIOMATIC_RG_CORE_BLOCKER",
    }

    out_json = ROOT / "report_qw2267_rg_effective_active_blocker_set_v2_gate.json"
    out_md = ROOT / "RAPORT_QW2267_RG_EFFECTIVE_ACTIVE_BLOCKER_SET_V2_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2267: RG EFFECTIVE ACTIVE BLOCKER-SET V2 GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- effective_active_blockers_v1: `{sorted(effective_v1)}`",
                f"- bridged_refs: `{sorted(bridged_refs)}`",
                f"- effective_active_blockers_v2_residual: `{residual}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "n_effective_active_blockers_v1": len(effective_v1),
                "n_effective_active_blockers_v2_residual": len(residual),
            }
        )
    )


if __name__ == "__main__":
    main()
