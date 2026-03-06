#!/usr/bin/env python3
"""QW-2385: dual noncircular anchor obligation packet gate.

Translates QW-2384 cycle diagnosis into minimal anchor obligations required to
avoid false progress and break theorem-cycle recurrence.
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
    q2384 = load("report_qw2384_dual_kernel_identity_cycle_structure_gate.json")

    obligations = [
        {
            "id": "L12_NONCIRCULAR_ANCHOR_O1",
            "branch": "L12",
            "anchor_symbol": "RG_KernelIdentityLocalityToWellPosedness_Theorem",
            "acceptance_rules": [
                "proof_term must not be single-symbol exact to any theorem inside blocker SCC",
                "proof must be machine-checkable without `axiom` and without `_DerivedOrPending` tokens",
                "dependency graph delta must reduce SCC cyclic closure rank on L12 branch",
            ],
        },
        {
            "id": "L5_NONCIRCULAR_ANCHOR_O1",
            "branch": "L5",
            "anchor_symbol": "QFT_KernelIdentityLocalityToPositivity_Theorem",
            "acceptance_rules": [
                "proof_term must not be single-symbol exact to any theorem inside blocker SCC",
                "proof must be machine-checkable without `axiom` and without `_DerivedOrPending` tokens",
                "dependency graph delta must reduce SCC cyclic closure rank on L5 branch",
            ],
        },
    ]

    flags = {
        "q2384_cycle_structure_confirmed": q2384.get("verdict")
        == "DUAL_KERNEL_IDENTITY_CYCLE_STRUCTURE_GATE_PASS_STRUCTURAL_CYCLE_CONFIRMED",
        "n_obligations_eq_two": len(obligations) == 2,
        "one_obligation_per_branch": {o["branch"] for o in obligations} == {"L12", "L5"},
        "all_rules_require_noncircularity": all(
            any("not be single-symbol exact" in rule for rule in o["acceptance_rules"]) for o in obligations
        ),
        "all_rules_require_axiom_free_machine_check": all(
            any("without `axiom`" in rule for rule in o["acceptance_rules"]) for o in obligations
        ),
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_KERNEL_IDENTITY_ANCHOR_OBLIGATION_PACKET_GATE_PASS_PACKET_READY"
        if flags["q2384_cycle_structure_confirmed"]
        and flags["n_obligations_eq_two"]
        and flags["one_obligation_per_branch"]
        and flags["all_rules_require_noncircularity"]
        and flags["all_rules_require_axiom_free_machine_check"]
        else "DUAL_KERNEL_IDENTITY_ANCHOR_OBLIGATION_PACKET_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2384_dual_kernel_identity_cycle_structure_gate.json",
        "obligations": obligations,
        "scope_boundary": {
            "packet_ready": True,
            "anchor_evidence_present": False,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    spec_path = ROOT / "spec_qw2385_dual_kernel_identity_anchor_obligation_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2384": "report_qw2384_dual_kernel_identity_cycle_structure_gate.json",
            "spec": spec_path.name,
        },
        "obligations": obligations,
        "n_obligations": len(obligations),
    }

    proof_path = ROOT / "proof_object_qw2385_dual_kernel_identity_anchor_obligation_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "obligations": obligations,
        "n_obligations": len(obligations),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "RUN_DUAL_KERNEL_IDENTITY_ANCHOR_EVIDENCE_ADMISSION_GATE",
    }

    out_json = ROOT / "report_qw2385_dual_kernel_identity_anchor_obligation_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2385_DUAL_KERNEL_IDENTITY_ANCHOR_OBLIGATION_PACKET_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2385: DUAL KERNEL IDENTITY ANCHOR OBLIGATION PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_obligations: `{len(obligations)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations": len(obligations)}))


if __name__ == "__main__":
    main()
