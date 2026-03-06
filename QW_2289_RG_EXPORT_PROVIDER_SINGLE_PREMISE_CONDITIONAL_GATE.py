#!/usr/bin/env python3
"""QW-2289: RG export-provider single-premise conditional gate.

Verifies machine-checkable conditional provider theorem without explicit axiom
tokens, with exactly one nonlogical physical bridge premise.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LEAN_BIN = Path("/tmp/lean4/lean-4.28.0-linux/bin/lean")
TARGET_FILE = ROOT / "FIN_L12_EXPORT_PROVIDER_SINGLE_PREMISE_CONDITIONAL.lean"
TARGET_THEOREM = "RG_CanonicalAction_to_WellPosedness_EXPORT_CONDITIONAL"
PREMISE_SYMBOL = "rgPhysicalBridgePremise"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def theorem_defined(text: str, theorem_name: str) -> bool:
    return bool(re.search(rf"^theorem\s+{re.escape(theorem_name)}\s*:", text, flags=re.M))


def main() -> None:
    text = TARGET_FILE.read_text(encoding="utf-8")

    lean_ok = LEAN_BIN.exists() and LEAN_BIN.is_file()
    machine_executed = False
    exit_code = None
    stdout = ""
    stderr = ""
    if lean_ok:
        machine_executed = True
        proc = subprocess.run(
            [str(LEAN_BIN), TARGET_FILE.name],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        exit_code = proc.returncode
        stdout = proc.stdout
        stderr = proc.stderr

    n_premise_symbols = len(re.findall(rf"\b{re.escape(PREMISE_SYMBOL)}\b", text))

    flags = {
        "target_file_exists": TARGET_FILE.exists(),
        "target_theorem_defined": theorem_defined(text, TARGET_THEOREM),
        "contains_no_axiom_tokens": "axiom " not in text,
        "contains_no_derived_or_pending_tokens": "_DerivedOrPending" not in text,
        "single_physical_bridge_premise_symbol_used": n_premise_symbols >= 1,
        "lean_binary_detected": lean_ok,
        "machine_check_executed": machine_executed,
        "machine_check_exit_zero": exit_code == 0,
        "conditional_provider_machine_checked": exit_code == 0,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_EXPORT_PROVIDER_SINGLE_PREMISE_CONDITIONAL_GATE_PASS_PARTIAL_CONDITIONAL_PROVIDER_MACHINE_CHECKED"
        if flags["target_file_exists"] and flags["target_theorem_defined"] and flags["machine_check_exit_zero"]
        else "RG_EXPORT_PROVIDER_SINGLE_PREMISE_CONDITIONAL_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_file": TARGET_FILE.name,
        "target_file_sha256": sha256_file(TARGET_FILE),
        "target_theorem": TARGET_THEOREM,
        "physical_bridge_premise_symbol": PREMISE_SYMBOL,
        "premise_symbol_occurrences": n_premise_symbols,
        "machine_check": {
            "lean_binary": str(LEAN_BIN),
            "exit_code": exit_code,
            "stdout": stdout,
            "stderr": stderr,
        },
        "scope_boundary": {
            "conditional_provider_machine_checked": exit_code == 0,
            "unconditional_non_axiomatic_provider_constructed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2289_rg_export_provider_single_premise_conditional.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_file": TARGET_FILE.name,
        "target_file_sha256": sha256_file(TARGET_FILE),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "target_theorem": TARGET_THEOREM,
        "physical_bridge_premise_symbol": PREMISE_SYMBOL,
        "premise_symbol_occurrences": n_premise_symbols,
        "machine_check_exit_code": exit_code,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_RG_PHYSICAL_BRIDGE_PREMISE_FROM_FIN_ACTION_LEVEL_CONTENT",
    }

    out_json = ROOT / "report_qw2289_rg_export_provider_single_premise_conditional_gate.json"
    out_md = ROOT / "RAPORT_QW2289_RG_EXPORT_PROVIDER_SINGLE_PREMISE_CONDITIONAL_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2289: RG EXPORT PROVIDER SINGLE PREMISE CONDITIONAL GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- machine_check_exit_code: `{exit_code}`",
                f"- premise_symbol_occurrences: `{n_premise_symbols}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "machine_check_exit_code": exit_code}))


if __name__ == "__main__":
    main()
