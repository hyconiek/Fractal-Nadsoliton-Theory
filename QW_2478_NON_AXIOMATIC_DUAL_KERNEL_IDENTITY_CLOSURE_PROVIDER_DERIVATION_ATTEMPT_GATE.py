#!/usr/bin/env python3
"""QW-2478: non-axiomatic dual kernel-identity-closure provider derivation attempt gate."""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

L12_ATTEMPT = ROOT / "FIN_L12_STRICT_NON_AXIOMATIC_KERNEL_IDENTITY_CLOSURE_PROVIDER_ATTEMPT.lean"
L5_ATTEMPT = ROOT / "FIN_L5_STRICT_NON_AXIOMATIC_KERNEL_IDENTITY_CLOSURE_PROVIDER_ATTEMPT.lean"

L12_EXPECTED_BLOCKER = "RG_KernelIdentityLocalityToWellPosedness_Theorem"
L5_EXPECTED_BLOCKER = "QFT_KernelIdentityLocalityToPositivity_Theorem"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def parse_unknowns(text: str) -> list[str]:
    pats = [
        re.compile(r"unknown identifier '([^']+)'"),
        re.compile(r"unknown constant '([^']+)'"),
        re.compile(r"Unknown identifier `([^`]+)`"),
    ]
    out: set[str] = set()
    for pat in pats:
        out.update(pat.findall(text))
    return sorted(out)


def file_is_axiom_token_free(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    has_axiom_token = bool(re.search(r"\baxiom\b", text))
    has_derived_pending = "_DerivedOrPending" in text
    return (not has_axiom_token) and (not has_derived_pending)


def run_lean(lean_bin: str, path: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["ELAN_HOME"] = str(ROOT / ".elan")
    env["HOME"] = str(ROOT / ".home_lean")
    env["PATH"] = f"{ROOT / '.elan/bin'}:{env.get('PATH', '')}"
    proc = subprocess.run([lean_bin, str(path.name)], cwd=ROOT, capture_output=True, text=True, check=False, env=env)
    merged = proc.stdout + "\n" + proc.stderr
    return {
        "file": path.name,
        "lean_bin": lean_bin,
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "unknown_identifiers": parse_unknowns(merged),
    }


def main() -> None:
    q2476 = load("report_qw2476_dual_kernel_identity_closure_provider_theorem_spec_gate.json")
    q2477 = load("report_qw2477_dual_kernel_identity_closure_provider_counterexample_search_gate.json")
    q2444 = load("report_qw2444_lean_runtime_discovery_gate.json")
    lean_bin = q2444.get("selected_runtime")

    l12_code = "\n".join(
        [
            "-- QW-2478 strict non-axiomatic kernel-identity-closure provider attempt (L12/RG)",
            "set_option autoImplicit false",
            "variable (FINActionComplete RGConstructiveMap RGGlobalWellPosednessAllScales : Prop)",
            "",
            "theorem RG_KernelIdentityClosureToWellPosedness_Theorem_NON_AXIOMATIC :",
            "  (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales := by",
            "  exact RG_KernelIdentityLocalityToWellPosedness_Theorem",
            "",
        ]
    )
    L12_ATTEMPT.write_text(l12_code, encoding="utf-8")

    l5_code = "\n".join(
        [
            "-- QW-2478 strict non-axiomatic kernel-identity-closure provider attempt (L5/QFT)",
            "set_option autoImplicit false",
            "variable (FINActionComplete ConstructiveNonPerturbativeScheme PositivityToReconstruction : Prop)",
            "",
            "theorem QFT_KernelIdentityClosureToPositivity_Theorem_NON_AXIOMATIC :",
            "  (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction := by",
            "  exact QFT_KernelIdentityLocalityToPositivity_Theorem",
            "",
        ]
    )
    L5_ATTEMPT.write_text(l5_code, encoding="utf-8")

    l12_axiom_free = file_is_axiom_token_free(L12_ATTEMPT)
    l5_axiom_free = file_is_axiom_token_free(L5_ATTEMPT)

    l12_run: dict[str, Any] | None = None
    l5_run: dict[str, Any] | None = None
    if lean_bin:
        l12_run = run_lean(str(lean_bin), L12_ATTEMPT)
        l5_run = run_lean(str(lean_bin), L5_ATTEMPT)

    flags = {
        "q2476_theorem_spec_ready": q2476.get("verdict")
        == "DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY",
        "q2477_counterexample_search_ready": q2477.get("verdict")
        == "DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_COUNTEREXAMPLE_SEARCH_GATE_PASS_PARTIAL_NO_STRICT_COUNTEREXAMPLE_IN_BOUNDED_DOMAIN",
        "q2444_runtime_available": q2444.get("verdict") == "LEAN_RUNTIME_DISCOVERY_GATE_PASS_RUNTIME_AVAILABLE",
        "l12_attempt_axiom_token_free": l12_axiom_free,
        "l5_attempt_axiom_token_free": l5_axiom_free,
        "execution_attempt_performed": l12_run is not None and l5_run is not None,
        "l12_machine_check_exit_nonzero": bool(l12_run and l12_run["exit_code"] != 0),
        "l5_machine_check_exit_nonzero": bool(l5_run and l5_run["exit_code"] != 0),
        "l12_expected_deeper_blocker_isolated": bool(
            l12_run and L12_EXPECTED_BLOCKER in set(l12_run.get("unknown_identifiers", []))
        ),
        "l5_expected_deeper_blocker_isolated": bool(
            l5_run and L5_EXPECTED_BLOCKER in set(l5_run.get("unknown_identifiers", []))
        ),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if not flags["q2444_runtime_available"]:
        verdict = "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"
        required_next_step = "ATTACH_RUNTIME_AND_RERUN_QW2478"
    elif (
        flags["q2476_theorem_spec_ready"]
        and flags["q2477_counterexample_search_ready"]
        and flags["execution_attempt_performed"]
        and flags["l12_attempt_axiom_token_free"]
        and flags["l5_attempt_axiom_token_free"]
        and flags["l12_machine_check_exit_nonzero"]
        and flags["l5_machine_check_exit_nonzero"]
        and flags["l12_expected_deeper_blocker_isolated"]
        and flags["l5_expected_deeper_blocker_isolated"]
    ):
        verdict = "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_IDENTITY_LOCALITY_PROVIDER_THEOREMS"
        required_next_step = "EXTRACT_DUAL_KERNEL_IDENTITY_LOCALITY_PROVIDER_MINIMAL_BLOCKER_CUT"
    elif (
        flags["q2476_theorem_spec_ready"]
        and flags["q2477_counterexample_search_ready"]
        and flags["execution_attempt_performed"]
        and flags["l12_attempt_axiom_token_free"]
        and flags["l5_attempt_axiom_token_free"]
        and (l12_run and l12_run["exit_code"] == 0)
        and (l5_run and l5_run["exit_code"] == 0)
    ):
        verdict = "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_LOCAL_DISCHARGE_READY_FOR_INTEGRATION"
        required_next_step = "INTEGRATE_LOCAL_NON_AXIOMATIC_KERNEL_IDENTITY_CLOSURE_PROVIDERS_IN_GLOBAL_CHAIN"
    else:
        verdict = "NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_FAIL"
        required_next_step = "REPAIR_QW2478_PIPELINE_AND_RERUN"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2476": "report_qw2476_dual_kernel_identity_closure_provider_theorem_spec_gate.json",
            "q2477": "report_qw2477_dual_kernel_identity_closure_provider_counterexample_search_gate.json",
            "q2444": "report_qw2444_lean_runtime_discovery_gate.json",
            "l12_attempt": L12_ATTEMPT.name,
            "l5_attempt": L5_ATTEMPT.name,
        },
        "runs": {"l12": l12_run, "l5": l5_run},
        "scope_boundary": {
            "local_non_axiomatic_attempt_only": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2478_non_axiomatic_dual_kernel_identity_closure_provider_derivation_attempt.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "unknown_identifiers": {
            "l12": (l12_run or {}).get("unknown_identifiers", []),
            "l5": (l5_run or {}).get("unknown_identifiers", []),
        },
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2478_non_axiomatic_dual_kernel_identity_closure_provider_derivation_attempt_gate.json"
    out_md = ROOT / "RAPORT_QW2478_NON_AXIOMATIC_DUAL_KERNEL_IDENTITY_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2478: NON AXIOMATIC DUAL KERNEL IDENTITY CLOSURE PROVIDER DERIVATION ATTEMPT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unknown_l12: `{out['unknown_identifiers']['l12']}`",
                f"- unknown_l5: `{out['unknown_identifiers']['l5']}`",
                "",
                "## Wniosek rygorystyczny",
                "- Attempt wykonany dopiero po theorem-spec i counterexample-search.",
                "- Attempt jest axiom-token-free i bez `_DerivedOrPending`.",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "unknown_identifiers": out["unknown_identifiers"]}))


if __name__ == "__main__":
    main()
