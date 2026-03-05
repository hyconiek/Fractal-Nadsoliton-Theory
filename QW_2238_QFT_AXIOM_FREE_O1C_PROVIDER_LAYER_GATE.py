#!/usr/bin/env python3
"""
QW-2238: QFT axiom-free O1c provider-layer gate.

Purpose:
- build explicit provider theorems QFT_C1_1_DERIVED / QFT_C1_2_DERIVED,
- verify machine-check for provider layer,
- keep axiomatic provenance boundary explicit.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
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


def detect_lean() -> str | None:
    direct = shutil.which("lean")
    if direct:
        return direct
    fallback = Path("/tmp/lean4/lean-4.28.0-linux/bin/lean")
    return str(fallback) if fallback.exists() else None


def run(cmd: list[str]) -> tuple[int, str, str]:
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    return int(proc.returncode), proc.stdout, proc.stderr


def main() -> None:
    q2234 = load("report_qw2234_qft_axiom_free_o1c_theorem_discharge_spec_gate.json")
    q2236 = load("report_qw2236_qft_axiom_free_o1c_theorem_discharge_execution_gate.json")

    provider = ROOT / "FIN_L5_O1C_PROVIDER_LAYER.lean"
    provider.write_text(
        "\n".join(
            [
                "-- FIN Release 5.1: QW-2238 QFT O1c provider layer",
                "-- Explicit provider theorems with axiomatic provenance boundary retained.",
                "",
                "axiom FINActionComplete : Prop",
                "axiom ConstructiveNonPerturbativeScheme : Prop",
                "axiom PositivityToReconstruction : Prop",
                "axiom UnitarySMatrixAndScatteringCompleteness : Prop",
                "",
                "axiom PositivityToReconstruction_DerivedOrPending :",
                "  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction",
                "axiom UnitarySMatrixAndScatteringCompleteness_DerivedOrPending :",
                "  PositivityToReconstruction -> UnitarySMatrixAndScatteringCompleteness",
                "",
                "theorem QFT_C1_1_DERIVED :",
                "  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by",
                "  intro h",
                "  exact PositivityToReconstruction_DerivedOrPending h",
                "",
                "theorem QFT_C1_2_DERIVED :",
                "  PositivityToReconstruction -> UnitarySMatrixAndScatteringCompleteness := by",
                "  intro h",
                "  exact UnitarySMatrixAndScatteringCompleteness_DerivedOrPending h",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    lean_bin = detect_lean()
    checker_attempted = False
    checker_ok = False
    checker_stdout = ""
    checker_stderr = ""
    checker_rc = None
    if lean_bin:
        checker_attempted = True
        rc, so, se = run([lean_bin, provider.name])
        checker_rc, checker_stdout, checker_stderr = rc, so, se
        checker_ok = rc == 0

    text = provider.read_text(encoding="utf-8")
    provider_theorems_present = (
        "theorem QFT_C1_1_DERIVED" in text and "theorem QFT_C1_2_DERIVED" in text
    )
    derived_or_pending_present = "DerivedOrPending" in text

    proof_object = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": {
            "q2234": "report_qw2234_qft_axiom_free_o1c_theorem_discharge_spec_gate.json",
            "q2236": "report_qw2236_qft_axiom_free_o1c_theorem_discharge_execution_gate.json",
        },
        "provider_file": {
            "name": provider.name,
            "sha256": sha256_file(provider),
        },
        "checker": {
            "attempted": checker_attempted,
            "exit_code": checker_rc,
            "stdout": checker_stdout,
            "stderr": checker_stderr,
        },
        "scope_boundary": {
            "provider_layer_axiomatic": True,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2238_qft_o1c_provider_layer.json"
    proof_path.write_text(json.dumps(proof_object, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2234_discharge_spec_pass_present": q2234.get("verdict")
        == "QFT_AXIOM_FREE_O1C_THEOREM_DISCHARGE_SPEC_GATE_PASS_PARTIAL_PROOFS_PENDING",
        "q2236_execution_blocker_present": q2236.get("verdict")
        == "QFT_AXIOM_FREE_O1C_THEOREM_DISCHARGE_EXECUTION_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_SOURCE_THEOREMS",
        "provider_file_written": provider.exists(),
        "provider_theorems_present": bool(provider_theorems_present),
        "derived_or_pending_axioms_present": bool(derived_or_pending_present),
        "lean_checker_detected": bool(lean_bin),
        "lean_checker_attempted": checker_attempted,
        "lean_checker_exit_zero": bool(checker_ok),
        "proof_object_generated": proof_path.exists(),
        "proof_object_hash_recorded": True,
        "provider_layer_axiomatic_boundary_explicit": True,
        "c1_theorem_discharge_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2234_discharge_spec_pass_present"]
        and flags["q2236_execution_blocker_present"]
        and flags["provider_file_written"]
        and flags["provider_theorems_present"]
        and flags["derived_or_pending_axioms_present"]
        and flags["lean_checker_detected"]
        and flags["lean_checker_attempted"]
        and flags["lean_checker_exit_zero"]
        and flags["proof_object_generated"]
        and flags["proof_object_hash_recorded"]
        and flags["provider_layer_axiomatic_boundary_explicit"]
    )

    verdict = (
        "QFT_AXIOM_FREE_O1C_PROVIDER_LAYER_GATE_PASS_PARTIAL_AXIOMATIC_PROVIDER_OPEN"
        if core_ok
        else "QFT_AXIOM_FREE_O1C_PROVIDER_LAYER_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2234": "report_qw2234_qft_axiom_free_o1c_theorem_discharge_spec_gate.json",
            "q2236": "report_qw2236_qft_axiom_free_o1c_theorem_discharge_execution_gate.json",
        },
        "provider_file": provider.name,
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "RERUN_THEOREM_DISCHARGE_EXECUTION_WITH_PROVIDER_LAYER_AND_CLASSIFY_AXIOMATIC_PROVENANCE_BOUNDARY",
    }

    out_json = ROOT / "report_qw2238_qft_axiom_free_o1c_provider_layer_gate.json"
    out_md = ROOT / "RAPORT_QW2238_QFT_AXIOM_FREE_O1C_PROVIDER_LAYER_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2238: QFT AXIOM-FREE O1C PROVIDER LAYER GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Dodano jawne provider theorems `QFT_C1_1_DERIVED`, `QFT_C1_2_DERIVED`.",
                "- Boundary pozostaje jawna: provider layer nadal oparty o `DerivedOrPending` axioms.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2238_qft_axiom_free_o1c_provider_layer_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
