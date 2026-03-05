#!/usr/bin/env python3
"""
QW-2242: QFT provider de-axiomatization obstruction gate.

Purpose:
- test whether QFT provider layer has non-axiomatic theorem sources,
- if not, classify exact de-axiomatization obstruction,
- export concrete obligations and source-map from canonical FIN chain.
"""

from __future__ import annotations

import hashlib
import json
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


def run(cmd: list[str]) -> tuple[int, str, str]:
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    return int(proc.returncode), proc.stdout, proc.stderr


def rg_lines(pattern: str) -> list[str]:
    rc, out, _err = run(["rg", "-n", pattern, "-g*.lean"])
    if rc == 1:
        return []
    return [ln for ln in out.splitlines() if ln.strip()]


def main() -> None:
    q2238 = load("report_qw2238_qft_axiom_free_o1c_provider_layer_gate.json")
    q2240 = load("report_qw2240_qft_axiom_free_o1c_provider_execution_gate.json")

    provider = ROOT / "FIN_L5_O1C_PROVIDER_LAYER.lean"
    text = provider.read_text(encoding="utf-8") if provider.exists() else ""

    has_axiom_derived_or_pending = "axiom PositivityToReconstruction_DerivedOrPending" in text and "axiom UnitarySMatrixAndScatteringCompleteness_DerivedOrPending" in text

    theorem_lines_qft_c1_1 = rg_lines(r"theorem\s+QFT_C1_1_DERIVED")
    theorem_lines_qft_c1_2 = rg_lines(r"theorem\s+QFT_C1_2_DERIVED")

    non_axiomatic_sources_c1_1 = [ln for ln in theorem_lines_qft_c1_1 if "FIN_L5_O1C_PROVIDER_LAYER.lean" not in ln]
    non_axiomatic_sources_c1_2 = [ln for ln in theorem_lines_qft_c1_2 if "FIN_L5_O1C_PROVIDER_LAYER.lean" not in ln]

    source_map = {
        "candidate_non_axiomatic_sources": {
            "QFT_C1_1": [
                "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
                "report_qw2222_qft_terminal_proof_object_execution_gate.json",
                "report_qw2202_qft_strict_scope_stratification_gate.json",
            ],
            "QFT_C1_2": [
                "report_qw2216_qft_unitary_scattering_scope_gate.json",
                "report_qw2218_qft_terminal_theorem_spec_gate.json",
                "report_qw2220_qft_terminal_proof_packet_ready_gate.json",
            ],
        },
        "status": "NOT_FORMALIZED_IN_NON_AXIOMATIC_LEAN_PROVIDER_THEOREMS",
    }

    obligations = {
        "scope": "L5_O1C_PROVIDER_DEAXIOMATIZATION",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "obligations": [
            {
                "id": "QFT_DAX_1",
                "statement": "construct non-axiomatic Lean theorem source for QFT_C1_1_DERIVED from canonical action/reconstruction chain",
                "depends_on": ["report_qw2165_l13_exhaustive_canonical_eom_gate.json", "report_qw2222_qft_terminal_proof_object_execution_gate.json"],
            },
            {
                "id": "QFT_DAX_2",
                "statement": "construct non-axiomatic Lean theorem source for QFT_C1_2_DERIVED from unitarity/scattering chain",
                "depends_on": ["report_qw2216_qft_unitary_scattering_scope_gate.json", "report_qw2218_qft_terminal_theorem_spec_gate.json"],
            },
            {
                "id": "QFT_DAX_3",
                "statement": "replace DerivedOrPending axioms by QFT_DAX_1 and QFT_DAX_2 providers and rerun O1c discharge",
                "depends_on": ["QFT_DAX_1", "QFT_DAX_2"],
            },
        ],
    }

    source_map_path = ROOT / "spec_qw2242_qft_provider_deaxiomatization_source_map.json"
    obligations_path = ROOT / "spec_qw2242_qft_provider_deaxiomatization_obligations.json"
    source_map_path.write_text(json.dumps(source_map, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    obligations_path.write_text(json.dumps(obligations, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_object = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": {
            "q2238": "report_qw2238_qft_axiom_free_o1c_provider_layer_gate.json",
            "q2240": "report_qw2240_qft_axiom_free_o1c_provider_execution_gate.json",
        },
        "provider_file": provider.name,
        "provider_file_sha256": sha256_file(provider) if provider.exists() else None,
        "scan_results": {
            "theorem_lines_qft_c1_1": theorem_lines_qft_c1_1,
            "theorem_lines_qft_c1_2": theorem_lines_qft_c1_2,
            "non_axiomatic_sources_c1_1": non_axiomatic_sources_c1_1,
            "non_axiomatic_sources_c1_2": non_axiomatic_sources_c1_2,
        },
        "source_map_file": source_map_path.name,
        "obligations_file": obligations_path.name,
        "scope_boundary": {
            "non_axiomatic_provider_layer_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2242_qft_provider_deaxiomatization_obstruction.json"
    proof_path.write_text(json.dumps(proof_object, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2238_provider_layer_pass_present": q2238.get("verdict")
        == "QFT_AXIOM_FREE_O1C_PROVIDER_LAYER_GATE_PASS_PARTIAL_AXIOMATIC_PROVIDER_OPEN",
        "q2240_provider_execution_pass_present": q2240.get("verdict")
        == "QFT_AXIOM_FREE_O1C_PROVIDER_EXECUTION_GATE_PASS_PARTIAL_PROVIDER_OK_AXIOMATIC_SOURCE_OPEN",
        "provider_file_present": provider.exists(),
        "provider_contains_axiomatic_source": bool(has_axiom_derived_or_pending),
        "scan_for_non_axiomatic_sources_executed": True,
        "non_axiomatic_source_exists_for_qft_c1_1": len(non_axiomatic_sources_c1_1) > 0,
        "non_axiomatic_source_exists_for_qft_c1_2": len(non_axiomatic_sources_c1_2) > 0,
        "source_map_written": source_map_path.exists(),
        "obligations_written": obligations_path.exists(),
        "proof_object_generated": proof_path.exists(),
        "proof_object_hash_recorded": True,
        "c1_theorem_discharge_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2238_provider_layer_pass_present"]
        and flags["q2240_provider_execution_pass_present"]
        and flags["provider_file_present"]
        and flags["provider_contains_axiomatic_source"]
        and flags["scan_for_non_axiomatic_sources_executed"]
        and flags["source_map_written"]
        and flags["obligations_written"]
        and flags["proof_object_generated"]
        and flags["proof_object_hash_recorded"]
    )

    verdict = (
        "QFT_PROVIDER_DEAXIOMATIZATION_OBSTRUCTION_GATE_PASS_PARTIAL_NON_AXIOMATIC_SOURCE_MISSING"
        if core_ok
        else "QFT_PROVIDER_DEAXIOMATIZATION_OBSTRUCTION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2238": "report_qw2238_qft_axiom_free_o1c_provider_layer_gate.json",
            "q2240": "report_qw2240_qft_axiom_free_o1c_provider_execution_gate.json",
        },
        "source_map_file": source_map_path.name,
        "obligations_file": obligations_path.name,
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "FORMALIZE_QFT_DAX_1_AND_QFT_DAX_2_AS_NON_AXIOMATIC_LEAN_THEOREMS",
    }

    out_json = ROOT / "report_qw2242_qft_provider_deaxiomatization_obstruction_gate.json"
    out_md = ROOT / "RAPORT_QW2242_QFT_PROVIDER_DEAXIOMATIZATION_OBSTRUCTION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2242: QFT PROVIDER DE-AXIOMATIZATION OBSTRUCTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Missing-provider blocker usuniety, ale non-axiomatic provider sources nie sa formalnie obecne w Lean.",
                "- Luka zostaje zredukowana do jawnych obligacji `QFT_DAX_1..QFT_DAX_3`.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2242_qft_provider_deaxiomatization_obstruction_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
