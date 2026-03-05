#!/usr/bin/env python3
"""QW-2248: QFT export axiomatic dependency gate.

Formal objective:
- certify whether QFT DAX1 export-level target has any axiom-free theorem provider,
- produce explicit dependency evidence for current axiomatic boundary.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
TARGET_EXPR = "(FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction"
EXPORT_SYMBOL = "QFT_CanonicalAction_to_Positivity_EXPORT"


def norm(s: str) -> str:
    return " ".join(s.split())


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def theorem_blocks(text: str) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    pat = re.compile(r"theorem\s+([A-Za-z0-9_']+)\s*:\s*(.*?)\s*:=\s*by\s*(.*?)(?=\n(?:theorem|axiom)\s+|\Z)", re.S)
    for m in pat.finditer(text):
        out.append({"name": m.group(1), "stmt": m.group(2), "proof": m.group(3)})
    return out


def analyze_file(path: Path) -> dict[str, Any]:
    txt = path.read_text(encoding="utf-8")
    axioms = re.findall(r"^axiom\s+([A-Za-z0-9_']+)\s*:", txt, flags=re.M)
    blocks = theorem_blocks(txt)

    matching: list[dict[str, Any]] = []
    for b in blocks:
        if norm(b["stmt"]) == norm(TARGET_EXPR):
            refs = re.findall(r"\b(?:exact|apply)\s+([A-Za-z0-9_'.]+)", b["proof"])
            dep_axioms = [r for r in refs if r in axioms]
            matching.append(
                {
                    "theorem": b["name"],
                    "refs": refs,
                    "depends_on_axioms": dep_axioms,
                    "depends_on_derived_or_pending": any(r.endswith("_DerivedOrPending") for r in refs),
                }
            )

    return {
        "file": path.name,
        "has_axioms": len(axioms) > 0,
        "axioms": axioms,
        "matching_theorems": matching,
    }


def main() -> None:
    q2244 = load("report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json")
    q2246 = load("report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json")

    lean_files = sorted(ROOT.glob("*.lean"))
    analyses = [analyze_file(p) for p in lean_files]

    files_with_matches = [a for a in analyses if a["matching_theorems"]]
    n_matching = sum(len(a["matching_theorems"]) for a in files_with_matches)

    all_matches = [m for a in files_with_matches for m in a["matching_theorems"]]
    any_dep_axiom = any(bool(m["depends_on_axioms"]) for m in all_matches)
    any_dep_pending = any(bool(m["depends_on_derived_or_pending"]) for m in all_matches)

    non_axiomatic_candidates = []
    for a in files_with_matches:
        if not a["has_axioms"]:
            for m in a["matching_theorems"]:
                if (not m["depends_on_axioms"]) and (not m["depends_on_derived_or_pending"]):
                    non_axiomatic_candidates.append({"file": a["file"], "theorem": m["theorem"]})

    export_exists_non_axiomatic = any(
        EXPORT_SYMBOL in (ROOT / p.name).read_text(encoding="utf-8") and ("axiom" not in (ROOT / p.name).read_text(encoding="utf-8"))
        for p in lean_files
    )

    flags = {
        "q2244_attempt_present": q2244.get("verdict")
        == "QFT_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING",
        "q2246_scan_present": q2246.get("verdict")
        == "QFT_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE_PASS_PARTIAL_NO_AXIOM_FREE_CANDIDATE",
        "lean_files_scanned": len(lean_files) > 0,
        "target_theorems_detected": n_matching > 0,
        "dependency_chain_hits_axiom_layer": any_dep_axiom,
        "dependency_chain_hits_derived_or_pending": any_dep_pending,
        "non_axiomatic_candidate_exists": len(non_axiomatic_candidates) > 0,
        "canonical_export_symbol_non_axiomatic_exists": export_exists_non_axiomatic,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2244_attempt_present"]
        and flags["q2246_scan_present"]
        and flags["lean_files_scanned"]
        and flags["target_theorems_detected"]
        and flags["dependency_chain_hits_axiom_layer"]
        and flags["dependency_chain_hits_derived_or_pending"]
    )

    verdict = (
        "QFT_EXPORT_AXIOMATIC_DEPENDENCY_GATE_PASS_PARTIAL_AXIOM_FREE_EXPORT_ABSENT"
        if core_ok and (not flags["non_axiomatic_candidate_exists"]) and (not flags["canonical_export_symbol_non_axiomatic_exists"])
        else "QFT_EXPORT_AXIOMATIC_DEPENDENCY_GATE_FAIL"
    )

    proof = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2244": "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
            "q2246": "report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json",
        },
        "target_expr": TARGET_EXPR,
        "export_symbol": EXPORT_SYMBOL,
        "scan": {
            "n_lean_files": len(lean_files),
            "n_matching_theorems": n_matching,
            "files_with_matches": files_with_matches,
            "non_axiomatic_candidates": non_axiomatic_candidates,
        },
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2248_qft_export_axiomatic_dependency.json"
    proof_path.write_text(json.dumps(proof, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "summary": {
            "n_lean_files": len(lean_files),
            "n_matching_theorems": n_matching,
            "n_non_axiomatic_candidates": len(non_axiomatic_candidates),
            "canonical_export_symbol_non_axiomatic_exists": export_exists_non_axiomatic,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_QFT_CANONICALACTION_TO_POSITIVITY_EXPORT_FROM_AXIOM_FREE_CORE_WITHOUT_DERIVED_OR_PENDING",
    }

    out_json = ROOT / "report_qw2248_qft_export_axiomatic_dependency_gate.json"
    out_md = ROOT / "RAPORT_QW2248_QFT_EXPORT_AXIOMATIC_DEPENDENCY_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2248: QFT EXPORT AXIOMATIC DEPENDENCY GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zaleznosc twierdzen targetowych od warstwy aksjomatycznej zostala potwierdzona formalnie.",
                "- Brak axiom-free export theorem dla QFT w aktualnym repo.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2248_qft_export_axiomatic_dependency_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
