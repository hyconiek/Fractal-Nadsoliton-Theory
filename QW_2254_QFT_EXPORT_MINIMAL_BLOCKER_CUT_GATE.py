#!/usr/bin/env python3
"""QW-2254: QFT export minimal blocker cut gate."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
TARGET_EXPR = "(FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction"


def norm(s: str) -> str:
    return " ".join(s.split())


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def theorem_blocks(text: str) -> list[dict[str, str]]:
    pat = re.compile(r"theorem\s+([A-Za-z0-9_']+)\s*:\s*(.*?)\s*:=\s*by\s*(.*?)(?=\n(?:theorem|axiom)\s+|\Z)", re.S)
    out = []
    for m in pat.finditer(text):
        out.append({"name": m.group(1), "stmt": m.group(2), "proof": m.group(3)})
    return out


def main() -> None:
    q2248 = load("report_qw2248_qft_export_axiomatic_dependency_gate.json")
    q2252 = load("report_qw2252_qft_export_obligation_execution_status_gate.json")

    lean_files = sorted(ROOT.glob("*.lean"))
    theorem_instances = []

    for p in lean_files:
        txt = p.read_text(encoding="utf-8")
        axioms = set(re.findall(r"^axiom\s+([A-Za-z0-9_']+)\s*:", txt, flags=re.M))
        for b in theorem_blocks(txt):
            if norm(b["stmt"]) != norm(TARGET_EXPR):
                continue
            refs = re.findall(r"\b(?:exact|apply)\s+([A-Za-z0-9_'.]+)", b["proof"])
            blockers = [
                r
                for r in refs
                if (r in axioms) or r.endswith("_DerivedOrPending") or r.endswith("Witness")
            ]
            theorem_instances.append(
                {
                    "file": p.name,
                    "theorem": b["name"],
                    "refs": refs,
                    "blockers": blockers,
                }
            )

    all_blockers = []
    for t in theorem_instances:
        all_blockers.extend(t["blockers"])

    unique_blockers = sorted(set(all_blockers))
    blocker_sets = [set(t["blockers"]) for t in theorem_instances if t["blockers"]]
    common_blockers = sorted(set.intersection(*blocker_sets)) if blocker_sets else []
    minimal_cut = common_blockers if common_blockers else unique_blockers

    flags = {
        "q2248_dependency_certificate_present": q2248.get("verdict")
        == "QFT_EXPORT_AXIOMATIC_DEPENDENCY_GATE_PASS_PARTIAL_AXIOM_FREE_EXPORT_ABSENT",
        "q2252_execution_status_present": q2252.get("verdict")
        == "QFT_EXPORT_OBLIGATION_EXECUTION_STATUS_GATE_PASS_PARTIAL_EXPORT_PENDING",
        "lean_files_scanned": len(lean_files) > 0,
        "target_theorem_instances_detected": len(theorem_instances) > 0,
        "blockers_detected": len(unique_blockers) > 0,
        "minimal_blocker_cut_computed": len(minimal_cut) > 0,
        "single_core_blocker_identified": len(common_blockers) == 1,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2248_dependency_certificate_present"]
        and flags["q2252_execution_status_present"]
        and flags["lean_files_scanned"]
        and flags["target_theorem_instances_detected"]
        and flags["blockers_detected"]
        and flags["minimal_blocker_cut_computed"]
    )

    verdict = (
        "QFT_EXPORT_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_BLOCKER_CUT_COMPUTED"
        if core_ok
        else "QFT_EXPORT_MINIMAL_BLOCKER_CUT_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2248": "report_qw2248_qft_export_axiomatic_dependency_gate.json",
            "q2252": "report_qw2252_qft_export_obligation_execution_status_gate.json",
        },
        "target_expr": TARGET_EXPR,
        "n_theorem_instances": len(theorem_instances),
        "theorem_instances": theorem_instances,
        "unique_blockers": unique_blockers,
        "common_blockers": common_blockers,
        "minimal_blocker_cut": minimal_cut,
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2254_qft_export_minimal_blocker_cut.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_theorem_instances": len(theorem_instances),
        "n_unique_blockers": len(unique_blockers),
        "n_common_blockers": len(common_blockers),
        "minimal_blocker_cut": minimal_cut,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "ELIMINATE_MINIMAL_BLOCKER_CUT_BY_PROVING_NON_AXIOMATIC_QFT_EXPORT",
    }

    out_json = ROOT / "report_qw2254_qft_export_minimal_blocker_cut_gate.json"
    out_md = ROOT / "RAPORT_QW2254_QFT_EXPORT_MINIMAL_BLOCKER_CUT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2254: QFT EXPORT MINIMAL BLOCKER CUT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- minimal_blocker_cut: `{minimal_cut}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "pass_count": pass_count,
                "total_flags": total_flags,
                "n_unique_blockers": len(unique_blockers),
                "minimal_blocker_cut": minimal_cut,
            }
        )
    )


if __name__ == "__main__":
    main()
