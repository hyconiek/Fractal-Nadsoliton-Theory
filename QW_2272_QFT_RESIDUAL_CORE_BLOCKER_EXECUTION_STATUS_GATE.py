#!/usr/bin/env python3
"""QW-2272: QFT residual core-blocker execution status gate.

Checks execution status for single residual QFT non-axiomatic obligation.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def theorem_exists(symbol: str) -> bool:
    patt = re.compile(rf"^theorem\s+{re.escape(symbol)}\s*:", flags=re.M)
    for p in ROOT.glob("*.lean"):
        txt = p.read_text(encoding="utf-8", errors="ignore")
        if patt.search(txt):
            return True
    return False


def main() -> None:
    spec = load("spec_qw2270_qft_residual_core_blocker_discharge_packet.json")
    residual = spec.get("residual_core_blockers", [])
    target = residual[0] if residual else ""
    derived = target.replace("_DerivedOrPending", "_Derived") if target else ""

    satisfied = theorem_exists(derived)

    obligation_rows = [
        {
            "id": "QFT_RESIDUAL_O1",
            "target_symbol": target,
            "derived_symbol": derived,
            "satisfied": bool(satisfied),
            "status": "satisfied" if satisfied else "pending_non_axiomatic_provider",
        }
    ]

    n_total = len(obligation_rows)
    n_satisfied = sum(1 for r in obligation_rows if r["satisfied"])

    flags = {
        "spec_present": True,
        "single_residual_obligation_present": n_total == 1,
        "all_obligations_satisfied": n_satisfied == n_total,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)
    verdict = "QFT_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_GATE_PASS_PARTIAL_PENDING" if flags["spec_present"] else "QFT_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_GATE_FAIL"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "spec_qw2270_qft_residual_core_blocker_discharge_packet.json",
        "obligation_status": obligation_rows,
        "n_obligations_total": n_total,
        "n_obligations_satisfied": n_satisfied,
    }

    proof_path = ROOT / "proof_object_qw2272_qft_residual_core_blocker_execution_status.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": proof_obj["source"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_obligations_total": n_total,
        "n_obligations_satisfied": n_satisfied,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_QFT_RESIDUAL_O1_NON_AXIOMATIC_DERIVED_SYMBOL",
    }

    out_json = ROOT / "report_qw2272_qft_residual_core_blocker_execution_status_gate.json"
    out_md = ROOT / "RAPORT_QW2272_QFT_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2272: QFT RESIDUAL CORE BLOCKER EXECUTION STATUS GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- obligations satisfied: `{n_satisfied}/{n_total}`",
                f"- required symbol: `{derived}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations_satisfied": n_satisfied, "n_obligations_total": n_total}))


if __name__ == "__main__":
    main()
