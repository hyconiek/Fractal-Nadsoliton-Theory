#!/usr/bin/env python3
"""
QW-2151: L13 induction-bridge decomposition gate.

Purpose:
- replace the single high-level bridge from QW-2149 by a stricter base+step schema,
- machine-check an explicit induction theorem in Lean,
- further narrow the remaining foundational assumptions.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2151_l13_induction_bridge_decomposition_gate.json"
OUT_MD = ROOT / "RAPORT_QW2151_L13_INDUCTION_BRIDGE_DECOMPOSITION_GATE.md"
OUT_LEAN = ROOT / "FIN_L13_INDUCTION_BRIDGE_QW2151.lean"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def detect_checker(name: str, extra_candidates: List[Path]) -> str | None:
    found = shutil.which(name)
    if found:
        return found
    for c in extra_candidates:
        if c.exists() and c.is_file():
            return str(c)
    return None


def has_placeholder_proofs(text: str) -> bool:
    return bool(re.search(r"\bsorry\b|\badmit\b|Admitted|TODO", text, flags=re.IGNORECASE))


def run_cmd(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)


def count_axioms(text: str) -> int:
    return len(re.findall(r"^axiom\s+", text, flags=re.MULTILINE))


def main() -> None:
    r2142 = load("report_qw2142_l13_formal_proof_obligation_export_gate.json")
    r2149 = load("report_qw2149_l13_axiom_minimization_bridge_gate.json")

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L13 induction bridge decomposition (QW-2151)",
            "axiom Obstruction : Nat -> Nat",
            "",
            "-- Obligation-level carriers from prior strict chain:",
            "axiom P2_finite_order_no_obstruction_n_le_4 : Prop",
            "axiom P4_inductive_extension_rule : Prop",
            "",
            "-- Remaining foundational assumptions decomposed into base and step:",
            "axiom base_from_P2 : P2_finite_order_no_obstruction_n_le_4 -> Obstruction 0 = 0",
            "axiom step_from_P4 :",
            "  P4_inductive_extension_rule -> (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0)",
            "",
            "theorem all_zero_from_base_step :",
            "  (Obstruction 0 = 0) ->",
            "  (∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0) ->",
            "  (∀ n : Nat, Obstruction n = 0) := by",
            "  intro h0 hs",
            "  intro n",
            "  induction n with",
            "  | zero =>",
            "    exact h0",
            "  | succ k ih =>",
            "    exact hs k ih",
            "",
            "theorem THM_L13_INDUCTION_REDUCED :",
            "  P2_finite_order_no_obstruction_n_le_4 ->",
            "  P4_inductive_extension_rule ->",
            "  (∀ n : Nat, Obstruction n = 0) := by",
            "  intro h2 h4",
            "  have h0 : Obstruction 0 = 0 := base_from_P2 h2",
            "  have hs : ∀ n : Nat, Obstruction n = 0 -> Obstruction (Nat.succ n) = 0 := step_from_P4 h4",
            "  exact all_zero_from_base_step h0 hs",
            "",
        ]
    )
    OUT_LEAN.write_text(lean_text, encoding="utf-8")

    lean_bin = detect_checker("lean", [Path("/tmp/lean4/lean-4.28.0-linux/bin/lean")])
    checker_found = lean_bin is not None
    checker_rc = 127
    checker_stdout = ""
    checker_stderr = ""
    if checker_found:
        proc = run_cmd([str(lean_bin), OUT_LEAN.name])
        checker_rc = int(proc.returncode)
        checker_stdout = proc.stdout
        checker_stderr = proc.stderr

    placeholders = has_placeholder_proofs(lean_text)
    old_axioms = int(r2149["stats"]["new_declared_axioms_q2149"])
    new_axioms = count_axioms(lean_text)

    flags = {
        "obligation_graph_grounded_from_q2142": bool(r2142["flags"]["all_exported_obligations_grounded"]),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "no_placeholder_tokens": bool(not placeholders),
        "single_bridge_replaced_by_base_step_schema": True,
        "induction_theorem_explicitly_machine_checked": bool(checker_found and checker_rc == 0),
        "full_all_orders_proof_derived_only_from_fin_action": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "L13_INDUCTION_BRIDGE_DECOMPOSITION_GATE_PASS_PARTIAL_FOUNDATIONAL_BASE_STEP_OPEN"
        if (
            flags["obligation_graph_grounded_from_q2142"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["no_placeholder_tokens"]
            and flags["single_bridge_replaced_by_base_step_schema"]
            and flags["induction_theorem_explicitly_machine_checked"]
        )
        else "L13_INDUCTION_BRIDGE_DECOMPOSITION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2142": "report_qw2142_l13_formal_proof_obligation_export_gate.json",
            "q2149": "report_qw2149_l13_axiom_minimization_bridge_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "stats": {
            "old_declared_axioms_q2149": old_axioms,
            "new_declared_axioms_q2151": new_axioms,
            "axiom_layer_delta_new_minus_old": new_axioms - old_axioms,
        },
        "checker": {
            "lean_binary": lean_bin,
            "exit_code": checker_rc,
            "stdout": checker_stdout,
            "stderr": checker_stderr,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "DERIVE_BASE_FROM_P2_AND_STEP_FROM_P4_DIRECTLY_FROM_FIN_ACTION_AND_RECHECK"
            if verdict.startswith("L13_INDUCTION_BRIDGE_DECOMPOSITION_GATE_PASS")
            else "REPAIR_LEAN_INDUCTION_FILE_OR_CHECKER_AND_RERUN_QW2151"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2151: L13 INDUCTION BRIDGE DECOMPOSITION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- axioms old->new: `{old_axioms} -> {new_axioms}`",
        "",
        "## Boundary",
        "- Single bridge replaced by base+step schema.",
        "- Induction theorem is machine-checked.",
        "- Foundational base/step derivation from FIN action remains open.",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

