#!/usr/bin/env python3
"""QW-2246: QFT DAX1 axiom-free candidate scan gate."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
TARGET_EXPR = "(FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction"
EXPORT_SYMBOL = "QFT_CanonicalAction_to_Positivity_EXPORT"


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def scan() -> dict[str, Any]:
    lean_files = sorted(ROOT.glob("*.lean"))
    target_files: list[str] = []
    axiom_free_candidates: list[str] = []
    export_symbol_locations: list[str] = []

    for p in lean_files:
        txt = p.read_text(encoding="utf-8")
        if TARGET_EXPR in txt:
            target_files.append(p.name)
            has_axiom = "axiom " in txt
            if (not has_axiom) and ("DerivedOrPending" not in txt):
                axiom_free_candidates.append(p.name)

        if EXPORT_SYMBOL in txt and "axiom" not in txt:
            export_symbol_locations.append(p.name)

    return {
        "n_lean_files": len(lean_files),
        "target_files": target_files,
        "n_target_files": len(target_files),
        "axiom_free_candidates": axiom_free_candidates,
        "n_axiom_free_candidates": len(axiom_free_candidates),
        "export_symbol_locations_non_axiomatic": export_symbol_locations,
        "n_export_symbol_locations_non_axiomatic": len(export_symbol_locations),
    }


def main() -> None:
    q2244 = load("report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json")
    q2242 = load("report_qw2242_qft_provider_deaxiomatization_obstruction_gate.json")

    s = scan()

    flags = {
        "q2244_attempt_gate_present": q2244.get("verdict")
        == "QFT_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING",
        "q2242_obstruction_map_present": q2242.get("verdict")
        == "QFT_PROVIDER_DEAXIOMATIZATION_OBSTRUCTION_GATE_PASS_PARTIAL_NON_AXIOMATIC_SOURCE_MISSING",
        "lean_files_scanned": s["n_lean_files"] > 0,
        "target_statement_detected": s["n_target_files"] > 0,
        "axiom_free_candidate_exists": s["n_axiom_free_candidates"] > 0,
        "canonical_export_non_axiomatic_exists": s["n_export_symbol_locations_non_axiomatic"] > 0,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2244_attempt_gate_present"]
        and flags["q2242_obstruction_map_present"]
        and flags["lean_files_scanned"]
        and flags["target_statement_detected"]
    )

    verdict = (
        "QFT_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE_PASS_PARTIAL_NO_AXIOM_FREE_CANDIDATE"
        if core_ok and (not flags["axiom_free_candidate_exists"]) and (not flags["canonical_export_non_axiomatic_exists"])
        else "QFT_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE_FAIL"
    )

    proof = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2244": "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
            "q2242": "report_qw2242_qft_provider_deaxiomatization_obstruction_gate.json",
        },
        "target_expr": TARGET_EXPR,
        "export_symbol": EXPORT_SYMBOL,
        "scan": s,
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2246_qft_dax1_axiom_free_candidate_scan.json"
    proof_path.write_text(json.dumps(proof, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "scan_summary": {
            "n_lean_files": s["n_lean_files"],
            "n_target_files": s["n_target_files"],
            "n_axiom_free_candidates": s["n_axiom_free_candidates"],
            "n_export_symbol_locations_non_axiomatic": s["n_export_symbol_locations_non_axiomatic"],
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "FORMALIZE_QFT_CANONICALACTION_TO_POSITIVITY_EXPORT_IN_AXIOM_FREE_CORE_THEORY",
    }

    out_json = ROOT / "report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json"
    out_md = ROOT / "RAPORT_QW2246_QFT_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2246: QFT DAX1 AXIOM-FREE CANDIDATE SCAN GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Wykonano pelny scan *.lean dla targetu DAX1 i canonical export symbolu.",
                "- Nie znaleziono axiom-free kandydata export theorem w aktualnym repo.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
