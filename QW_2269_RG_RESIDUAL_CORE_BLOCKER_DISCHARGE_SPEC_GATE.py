#!/usr/bin/env python3
"""QW-2269: RG residual core-blocker discharge spec gate.

Builds explicit single-obligation discharge packet for the residual RG non-axiomatic
core blocker identified by QW-2267.
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
    p2267 = load("proof_object_qw2267_rg_effective_active_blocker_set_v2.json")
    residual = p2267.get("effective_active_blockers_v2_residual", [])

    expected = "RGGlobalWellPosednessAllScales_DerivedOrPending"

    obligations: list[dict[str, Any]] = []
    if residual:
        target = residual[0]
        obligations.append(
            {
                "id": "RG_RESIDUAL_O1",
                "target_symbol": target,
                "target_kind": "non_axiomatic_provider_replacement",
                "required_outcome": {
                    "introduce_symbol": target.replace("_DerivedOrPending", "_Derived"),
                    "remove_active_path_dependency_on": target,
                    "keep_o1c_scope_boundary": True,
                },
                "acceptance_criteria": [
                    "NON_AXIOMATIC_DERIVED_SYMBOL_EXISTS",
                    "ACTIVE_PATH_NO_LONGER_REFERENCES_DERIVED_OR_PENDING_SYMBOL",
                    "LOCALITY_INTEGRITY_NO_DANGLING_REFS",
                ],
            }
        )

    flags = {
        "q2267_residual_set_present": len(residual) > 0,
        "residual_singleton": len(residual) == 1,
        "residual_matches_expected_core_symbol": residual == [expected],
        "single_obligation_packet_emitted": len(obligations) == 1,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2267_residual_set_present"]
        and flags["residual_singleton"]
        and flags["single_obligation_packet_emitted"]
    )

    verdict = (
        "RG_RESIDUAL_CORE_BLOCKER_DISCHARGE_SPEC_GATE_PASS_SINGLE_OBLIGATION_PACKET_READY"
        if core_ok
        else "RG_RESIDUAL_CORE_BLOCKER_DISCHARGE_SPEC_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "proof_object_qw2267_rg_effective_active_blocker_set_v2.json",
        "residual_core_blockers": residual,
        "obligations": obligations,
        "scope_boundary": {
            "non_axiomatic_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    spec_path = ROOT / "spec_qw2269_rg_residual_core_blocker_discharge_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "p2267": "proof_object_qw2267_rg_effective_active_blocker_set_v2.json",
            "spec": spec_path.name,
        },
        "residual_core_blockers": residual,
        "n_residual_core_blockers": len(residual),
        "n_obligations": len(obligations),
    }

    proof_path = ROOT / "proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "residual_core_blockers": residual,
        "n_obligations": len(obligations),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "EXECUTE_RG_RESIDUAL_O1_NON_AXIOMATIC_PROVIDER_DISCHARGE",
    }

    out_json = ROOT / "report_qw2269_rg_residual_core_blocker_discharge_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2269_RG_RESIDUAL_CORE_BLOCKER_DISCHARGE_SPEC_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2269: RG RESIDUAL CORE BLOCKER DISCHARGE SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- residual_core_blockers: `{residual}`",
                f"- n_obligations: `{len(obligations)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "n_residual_core_blockers": len(residual),
                "n_obligations": len(obligations),
            }
        )
    )


if __name__ == "__main__":
    main()
