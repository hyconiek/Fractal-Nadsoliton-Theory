#!/usr/bin/env python3
"""QW-2393: dual foundational noncyclic anchor-obligation packet gate.

Defines minimal noncyclic obligations for foundational layer after QW-2392.
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
    q2392 = load("report_qw2392_dual_foundational_chain_reuse_admission_gate.json")

    obligations = [
        {
            "id": "L12_FOUNDATIONAL_NONCYCLIC_ANCHOR_O1",
            "branch": "L12",
            "anchor_symbol": "RG_FundamentalActionToWellPosedness_Derivation",
            "acceptance_rules": [
                "proof_term must not be single-symbol exact to historical foundational chain step",
                "proof must be machine-checkable without `axiom` and without `_DerivedOrPending`",
                "dependency graph delta must not re-enter confirmed identity-cycle frontier",
            ],
        },
        {
            "id": "L5_FOUNDATIONAL_NONCYCLIC_ANCHOR_O1",
            "branch": "L5",
            "anchor_symbol": "QFT_FundamentalActionToPositivity_Derivation",
            "acceptance_rules": [
                "proof_term must not be single-symbol exact to historical foundational chain step",
                "proof must be machine-checkable without `axiom` and without `_DerivedOrPending`",
                "dependency graph delta must not re-enter confirmed identity-cycle frontier",
            ],
        },
    ]

    q2392_verdict = str(q2392.get("verdict", ""))
    flags = {
        "q2392_reuse_guard_active": bool(q2392.get("flags", {}).get("reuse_without_new_evidence_forbidden", False)),
        "q2392_admission_state_valid": q2392_verdict
        in {
            "DUAL_FOUNDATIONAL_CHAIN_REUSE_ADMISSION_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_FOUNDATIONAL_ANCHOR_CANDIDATES",
            "DUAL_FOUNDATIONAL_CHAIN_REUSE_ADMISSION_GATE_PASS_ADMITTED",
        },
        "n_obligations_eq_two": len(obligations) == 2,
        "one_obligation_per_branch": {o["branch"] for o in obligations} == {"L12", "L5"},
        "all_rules_require_noncircularity": all(
            any("must not be single-symbol exact" in r for r in o["acceptance_rules"]) for o in obligations
        ),
        "all_rules_require_axiom_free_machine_check": all(
            any("without `axiom`" in r for r in o["acceptance_rules"]) for o in obligations
        ),
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_FOUNDATIONAL_NONCYCLIC_ANCHOR_OBLIGATION_PACKET_GATE_PASS_PACKET_READY"
        if flags["q2392_reuse_guard_active"]
        and flags["q2392_admission_state_valid"]
        and flags["n_obligations_eq_two"]
        and flags["one_obligation_per_branch"]
        and flags["all_rules_require_noncircularity"]
        and flags["all_rules_require_axiom_free_machine_check"]
        else "DUAL_FOUNDATIONAL_NONCYCLIC_ANCHOR_OBLIGATION_PACKET_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2392_dual_foundational_chain_reuse_admission_gate.json",
        "obligations": obligations,
        "scope_boundary": {
            "packet_ready": True,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    spec_path = ROOT / "spec_qw2393_dual_foundational_noncyclic_anchor_obligation_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2392": "report_qw2392_dual_foundational_chain_reuse_admission_gate.json",
            "spec": spec_path.name,
        },
        "obligations": obligations,
        "n_obligations": len(obligations),
    }

    proof_path = ROOT / "proof_object_qw2393_dual_foundational_noncyclic_anchor_obligation_packet.json"
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
        "required_next_step": "AUTHOR_DUAL_FOUNDATIONAL_NONCYCLIC_ANCHOR_CANDIDATE_FILES",
    }

    out_json = ROOT / "report_qw2393_dual_foundational_noncyclic_anchor_obligation_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2393_DUAL_FOUNDATIONAL_NONCYCLIC_ANCHOR_OBLIGATION_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2393: DUAL FOUNDATIONAL NONCYCLIC ANCHOR OBLIGATION PACKET GATE",
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
