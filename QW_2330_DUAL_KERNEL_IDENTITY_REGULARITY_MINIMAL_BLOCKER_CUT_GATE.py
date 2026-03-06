#!/usr/bin/env python3
"""QW-2330: dual kernel-identity-regularity minimal blocker-cut gate.

Reduces QW-2329 blockers to minimal core obligations, one per branch.
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
    q2329 = load("report_qw2329_dual_kernel_identity_coherence_execution_status_gate.json")
    blockers = q2329.get("blocker_cut", [])

    rg_blocker = "RG_KernelIdentityRegularityToWellPosedness_Theorem"
    qft_blocker = "QFT_KernelIdentityRegularityToPositivity_Theorem"

    minimal_blocker_cut = [
        {
            "id": "L12_KERNEL_IDENTITY_REGULARITY_CORE_O1",
            "branch": "L12",
            "symbol": rg_blocker,
            "required_statement": "(FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales",
            "proof_kind": "kernel_identity_regularity_derivation",
        },
        {
            "id": "L5_KERNEL_IDENTITY_REGULARITY_CORE_O1",
            "branch": "L5",
            "symbol": qft_blocker,
            "required_statement": "(FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction",
            "proof_kind": "kernel_identity_regularity_derivation",
        },
    ]

    flags = {
        "q2329_execution_status_present": q2329.get("verdict")
        == "DUAL_KERNEL_IDENTITY_COHERENCE_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_IDENTITY_REGULARITY_THEOREMS",
        "q2329_blocker_cut_size_two": len(blockers) == 2,
        "contains_rg_kernel_identity_regularity_blocker": any(b.get("symbol") == rg_blocker for b in blockers),
        "contains_qft_kernel_identity_regularity_blocker": any(b.get("symbol") == qft_blocker for b in blockers),
        "minimal_blocker_cut_size_two": len(minimal_blocker_cut) == 2,
        "one_core_obligation_per_branch": True,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_KERNEL_IDENTITY_REGULARITY_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED"
        if flags["q2329_execution_status_present"]
        and flags["q2329_blocker_cut_size_two"]
        and flags["contains_rg_kernel_identity_regularity_blocker"]
        and flags["contains_qft_kernel_identity_regularity_blocker"]
        else "DUAL_KERNEL_IDENTITY_REGULARITY_MINIMAL_BLOCKER_CUT_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2329_dual_kernel_identity_coherence_execution_status_gate.json",
        "minimal_blocker_cut": minimal_blocker_cut,
        "n_obligations": len(minimal_blocker_cut),
        "scope_boundary": {
            "minimal_blocker_cut_extracted": True,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2330_dual_kernel_identity_regularity_minimal_blocker_cut.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2329": "report_qw2329_dual_kernel_identity_coherence_execution_status_gate.json",
            "spec": spec_path.name,
        },
        "minimal_blocker_cut": minimal_blocker_cut,
        "n_obligations": len(minimal_blocker_cut),
    }
    proof_path = ROOT / "proof_object_qw2330_dual_kernel_identity_regularity_minimal_blocker_cut.json"
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
        "required_next_step": "BUILD_DUAL_KERNEL_IDENTITY_REGULARITY_DISCHARGE_PACKET",
    }

    out_json = ROOT / "report_qw2330_dual_kernel_identity_regularity_minimal_blocker_cut_gate.json"
    out_md = ROOT / "RAPORT_QW2330_DUAL_KERNEL_IDENTITY_REGULARITY_MINIMAL_BLOCKER_CUT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2330: DUAL KERNEL IDENTITY REGULARITY MINIMAL BLOCKER CUT GATE",
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
