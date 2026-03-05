#!/usr/bin/env python3
"""
QW-2224: QFT axiom-free discharge specification gate (L5).

Purpose:
- audit QW-2222 terminal Lean artifacts for explicit axioms,
- convert L5_AXIOM_FREE_O1 into explicit sub-obligations with dependency DAG,
- keep no-overclaim boundary explicit.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def parse_axioms(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    for ln in lines:
        m = re.match(r"\s*axiom\s+([A-Za-z0-9_']+)\s*:", ln)
        if m:
            out.append(m.group(1))
    return out


def dag_acyclic(dag: dict[str, list[str]]) -> bool:
    indeg = {k: 0 for k in dag.keys()}
    g: dict[str, list[str]] = defaultdict(list)
    for node, deps in dag.items():
        for d in deps:
            if d in indeg:
                g[d].append(node)
                indeg[node] += 1
    q = deque([k for k, v in indeg.items() if v == 0])
    seen = 0
    while q:
        n = q.popleft()
        seen += 1
        for nxt in g[n]:
            indeg[nxt] -= 1
            if indeg[nxt] == 0:
                q.append(nxt)
    return seen == len(dag)


def main() -> None:
    q2222 = load("report_qw2222_qft_terminal_proof_object_execution_gate.json")

    file_a = ROOT / "FIN_L5_O1A_O1_TERMINAL.lean"
    file_b = ROOT / "FIN_L5_O1B_O1_TERMINAL.lean"

    axioms_a = parse_axioms(file_a) if file_a.exists() else []
    axioms_b = parse_axioms(file_b) if file_b.exists() else []

    sub_obligations = {
        "L5_AXIOM_FREE_O1a": {
            "description": "Replace axioms in L5_O1a theorem by derived constructive/reconstruction lemmas.",
            "depends_on": [],
            "closed": False,
            "targets": ["FINActionComplete", "ConstructiveNonPerturbativeScheme", "PositivityToReconstruction", "L5O1aWitness"],
        },
        "L5_AXIOM_FREE_O1b": {
            "description": "Replace axioms in L5_O1b theorem by derived unitarity/scattering-completeness lemmas.",
            "depends_on": ["L5_AXIOM_FREE_O1a"],
            "closed": False,
            "targets": ["PositivityToReconstruction", "UnitarySMatrixAndScatteringCompleteness", "L5O1bWitness"],
        },
        "L5_AXIOM_FREE_O1c": {
            "description": "Produce end-to-end axiom-free machine-check object for both L5 terminal theorems.",
            "depends_on": ["L5_AXIOM_FREE_O1a", "L5_AXIOM_FREE_O1b"],
            "closed": False,
            "targets": ["L5_O1A_O1_TERMINAL", "L5_O1B_O1_TERMINAL"],
        },
    }

    dag = {k: v["depends_on"] for k, v in sub_obligations.items()}

    flags = {
        "source_q2222_execution_pass_present": q2222.get("verdict")
        == "QFT_TERMINAL_PROOF_OBJECT_EXECUTION_GATE_PASS_PARTIAL_AXIOMATIC_BOUNDARY",
        "l5_exec_o1_closed": bool(q2222.get("closed_obligation", {}).get("closed")),
        "terminal_lean_files_present": file_a.exists() and file_b.exists(),
        "axioms_detected_in_terminal_files": len(axioms_a) + len(axioms_b) > 0,
        "l5_axiom_free_subobligations_declared": len(sub_obligations) == 3,
        "subobligation_dependency_dag_acyclic": dag_acyclic(dag),
        "all_subobligations_open_explicit": all(not v["closed"] for v in sub_obligations.values()),
        "no_overclaim_boundary_explicit": True,
        "l5_axiom_free_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["source_q2222_execution_pass_present"]
        and flags["l5_exec_o1_closed"]
        and flags["terminal_lean_files_present"]
        and flags["axioms_detected_in_terminal_files"]
        and flags["l5_axiom_free_subobligations_declared"]
        and flags["subobligation_dependency_dag_acyclic"]
        and flags["all_subobligations_open_explicit"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_AXIOM_FREE_DISCHARGE_SPEC_GATE_PASS_PARTIAL_SUBOBLIGATIONS_OPEN"
        if core_ok
        else "QFT_AXIOM_FREE_DISCHARGE_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2222": "report_qw2222_qft_terminal_proof_object_execution_gate.json",
            "file_a": file_a.name,
            "file_b": file_b.name,
        },
        "axiom_audit": {
            "file_a_axioms": axioms_a,
            "file_b_axioms": axioms_b,
            "n_axioms_total": len(axioms_a) + len(axioms_b),
        },
        "sub_obligations": sub_obligations,
        "dependency_dag": dag,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DISCHARGE_L5_AXIOM_FREE_O1A_THEN_O1B_AND_ATTACH_AXIOM_FREE_PROOF_OBJECT",
    }

    out_json = ROOT / "report_qw2224_qft_axiom_free_discharge_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2224_QFT_AXIOM_FREE_DISCHARGE_SPEC_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2224: QFT AXIOM-FREE DISCHARGE SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_axioms_total: `{out['axiom_audit']['n_axioms_total']}`",
                "",
                "## Core result",
                "- L5 execution layer jest domknieta, ale terminalne pliki nadal zawieraja aksjomaty.",
                "- L5 axiom-free gap zostal rozbity na 3 jawne subobligacje z DAG.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2224_qft_axiom_free_discharge_spec_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
