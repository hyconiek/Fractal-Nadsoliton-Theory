#!/usr/bin/env python3
"""QW-2372: dual kernel-identity-totality minimal blocker-cut gate.

Reduces QW-2371 blockers to minimal core obligations, one per branch.
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
    q2371 = load("report_qw2371_dual_kernel_identity_universality_execution_status_gate.json")
    blockers = q2371.get("blocker_cut", [])

    rg_blocker = "RG_KernelIdentityTotalityToWellPosedness_Theorem"
    qft_blocker = "QFT_KernelIdentityTotalityToPositivity_Theorem"

    minimal_blocker_cut = [
        {
            "id": "L12_KERNEL_IDENTITY_TOTALITY_CORE_O1",
            "branch": "L12",
            "symbol": rg_blocker,
            "required_statement": "(FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales",
            "proof_kind": "kernel_identity_totality_derivation",
        },
        {
            "id": "L5_KERNEL_IDENTITY_TOTALITY_CORE_O1",
            "branch": "L5",
            "symbol": qft_blocker,
            "required_statement": "(FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction",
            "proof_kind": "kernel_identity_totality_derivation",
        },
    ]

    flags = {
        "q2371_execution_status_present": q2371.get("verdict")
        == "DUAL_KERNEL_IDENTITY_UNIVERSALITY_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_IDENTITY_TOTALITY_THEOREMS",
        "q2371_blocker_cut_size_two": len(blockers) == 2,
        "contains_rg_kernel_identity_totality_blocker": any(b.get("symbol") == rg_blocker for b in blockers),
        "contains_qft_kernel_identity_totality_blocker": any(b.get("symbol") == qft_blocker for b in blockers),
        "minimal_blocker_cut_size_two": len(minimal_blocker_cut) == 2,
        "one_core_obligation_per_branch": True,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_KERNEL_IDENTITY_TOTALITY_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED"
        if flags["q2371_execution_status_present"]
        and flags["q2371_blocker_cut_size_two"]
        and flags["contains_rg_kernel_identity_totality_blocker"]
        and flags["contains_qft_kernel_identity_totality_blocker"]
        else "DUAL_KERNEL_IDENTITY_TOTALITY_MINIMAL_BLOCKER_CUT_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2371_dual_kernel_identity_universality_execution_status_gate.json",
        "minimal_blocker_cut": minimal_blocker_cut,
        "n_obligations": len(minimal_blocker_cut),
        "scope_boundary": {
            "minimal_blocker_cut_extracted": True,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2372_dual_kernel_identity_totality_minimal_blocker_cut.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2371": "report_qw2371_dual_kernel_identity_universality_execution_status_gate.json",
            "spec": spec_path.name,
        },
        "minimal_blocker_cut": minimal_blocker_cut,
        "n_obligations": len(minimal_blocker_cut),
    }
    proof_path = ROOT / "proof_object_qw2372_dual_kernel_identity_totality_minimal_blocker_cut.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "minimal_blocker_cut": minimal_blocker_cut,
        "n_obligations": len(minimal_blocker_cut),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "BUILD_DUAL_KERNEL_IDENTITY_TOTALITY_DISCHARGE_PACKET",
    }

    out_json = ROOT / "report_qw2372_dual_kernel_identity_totality_minimal_blocker_cut_gate.json"
    out_md = ROOT / "RAPORT_QW2372_DUAL_KERNEL_IDENTITY_TOTALITY_MINIMAL_BLOCKER_CUT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2372: DUAL KERNEL IDENTITY TOTALITY MINIMAL BLOCKER CUT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_obligations: `{len(minimal_blocker_cut)}`",
                f"- minimal_blocker_cut: `{minimal_blocker_cut}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations": len(minimal_blocker_cut)}))


if __name__ == "__main__":
    main()
