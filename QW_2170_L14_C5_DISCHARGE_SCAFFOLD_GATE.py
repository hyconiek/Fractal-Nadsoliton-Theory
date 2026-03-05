#!/usr/bin/env python3
"""
QW-2170: L14 C5 discharge scaffold gate.

Purpose:
- decompose final open C5 step into explicit sub-obligations,
- show finite-grid/continuum-support part is closed in strict chain,
- keep explicit terminal open boundary on exact continuum theorem from complete
  FIN action.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2170_l14_c5_discharge_scaffold_gate.json"
OUT_MD = ROOT / "RAPORT_QW2170_L14_C5_DISCHARGE_SCAFFOLD_GATE.md"
OUT_LEAN = ROOT / "FIN_L14_C5_DISCHARGE_SCAFFOLD_QW2170.lean"
OUT_PACKET = ROOT / "proof_packet_qw2170_l14_c5_discharge_scaffold.json"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def detect_checker(name: str, extra_candidates: List[Path]) -> str | None:
    found = shutil.which(name)
    if found:
        return found
    for c in extra_candidates:
        if c.exists() and c.is_file():
            return str(c)
    return None


def run_cmd(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)


def has_placeholder(text: str) -> bool:
    return bool(re.search(r"\bsorry\b|\badmit\b|Admitted|TODO", text, flags=re.IGNORECASE))


def topo_sort(nodes: List[str], edges: List[List[str]]) -> List[str]:
    indeg = {n: 0 for n in nodes}
    adj = {n: [] for n in nodes}
    for a, b in edges:
        adj[a].append(b)
        indeg[b] += 1
    q = [n for n in nodes if indeg[n] == 0]
    order = []
    while q:
        x = q.pop(0)
        order.append(x)
        for y in adj[x]:
            indeg[y] -= 1
            if indeg[y] == 0:
                q.append(y)
    return order


def main() -> None:
    r2168 = load("report_qw2168_l14_final_continuum_theorem_packet_gate.json")
    r2166 = load("report_qw2166_l14_exhaustive_canonical_hessian_gate.json")
    r2148 = load("report_qw2148_continuum_dg_delta_extrapolation_gate.json")
    r2141 = load("report_qw2141_continuum_weak_distribution_proxy_gate.json")

    c5a_closed = bool(
        r2166["flags"]["exhaustive_continuum_bundle_closed_on_canonical_hessian_level"]
        and r2148["flags"]["weak_distribution_proxy_regularized_error_monotone_with_volume"]
        and r2148["flags"]["periodic_proxy_continuum_support_established"]
        and r2141["flags"]["exact_pairing_identity_all_cases"]
    )

    c5b_closed = False
    c5_closed = c5a_closed and c5b_closed

    obligations = [
        {
            "id": "C5a",
            "name": "finite_grid_to_continuum_support_from_strict_chain",
            "status": "satisfied" if c5a_closed else "open",
            "source": "QW-2166/QW-2148/QW-2141",
        },
        {
            "id": "C5b",
            "name": "exact_distributional_limit_directly_from_complete_fin_action",
            "status": "open",
            "source": "terminal theorem",
        },
        {
            "id": "C5",
            "name": "final_continuum_theorem_from_complete_fin_action",
            "status": "open",
            "source": "C5a + C5b",
        },
    ]

    nodes = [o["id"] for o in obligations]
    edges = [["C5a", "C5"], ["C5b", "C5"]]
    topo = topo_sort(nodes, edges)
    acyclic = len(topo) == len(nodes)

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L14 C5 discharge scaffold (QW-2170)",
            "",
            "theorem l14_c5_discharge_scaffold",
            "  {C5a C5b : Prop}",
            "  (ha : C5a) (hb : C5b) :",
            "  C5a ∧ C5b := by",
            "  exact And.intro ha hb",
            "",
            "theorem l14_c5_composition",
            "  {C5a C5b C5 : Prop}",
            "  (hcomp : C5a -> C5b -> C5)",
            "  (ha : C5a) (hb : C5b) : C5 := by",
            "  exact hcomp ha hb",
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

    placeholders = has_placeholder(lean_text)

    packet = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "packet_name": "QW2170_L14_C5_DISCHARGE_SCAFFOLD",
        "inputs": {
            "q2168": "report_qw2168_l14_final_continuum_theorem_packet_gate.json",
            "q2166": "report_qw2166_l14_exhaustive_canonical_hessian_gate.json",
            "q2148": "report_qw2148_continuum_dg_delta_extrapolation_gate.json",
            "q2141": "report_qw2141_continuum_weak_distribution_proxy_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "obligations": obligations,
        "dependency_graph": {"nodes": nodes, "edges": edges, "topological_order": topo},
        "decomposition": {
            "c5a_finite_continuum_support_closed": c5a_closed,
            "c5b_exact_distribution_limit_closed": c5b_closed,
            "c5_closed": c5_closed,
        },
        "hashes": {"lean_sha256": sha256_file(OUT_LEAN)},
    }
    OUT_PACKET.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "q2168_final_packet_present": bool(r2168["flags"]["final_obligation_vector_declared"]),
        "q2166_exhaustive_hessian_layer_present": bool(
            r2166["flags"]["exhaustive_continuum_bundle_closed_on_canonical_hessian_level"]
        ),
        "c5a_finite_continuum_support_closed": bool(c5a_closed),
        "c5_decomposition_declared": True,
        "c5_dependency_graph_acyclic": bool(acyclic),
        "c5_composition_schema_declared": True,
        "no_placeholder_tokens_in_scaffold_theorem": bool(not placeholders),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "proof_packet_manifest_written": OUT_PACKET.exists(),
        "terminal_c5b_exact_distribution_limit_closed": False,
        "full_continuum_theorem_from_complete_fin_action_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L14_C5_DISCHARGE_SCAFFOLD_GATE_PASS_PARTIAL_TERMINAL_BOUND_OPEN"
        if (
            flags["q2168_final_packet_present"]
            and flags["q2166_exhaustive_hessian_layer_present"]
            and flags["c5a_finite_continuum_support_closed"]
            and flags["c5_decomposition_declared"]
            and flags["c5_dependency_graph_acyclic"]
            and flags["c5_composition_schema_declared"]
            and flags["no_placeholder_tokens_in_scaffold_theorem"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["proof_packet_manifest_written"]
        )
        else "L14_C5_DISCHARGE_SCAFFOLD_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2168": "report_qw2168_l14_final_continuum_theorem_packet_gate.json",
            "q2166": "report_qw2166_l14_exhaustive_canonical_hessian_gate.json",
            "q2148": "report_qw2148_continuum_dg_delta_extrapolation_gate.json",
            "q2141": "report_qw2141_continuum_weak_distribution_proxy_gate.json",
            "proof_packet": OUT_PACKET.name,
            "lean_file": OUT_LEAN.name,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "checker": {
            "lean_binary": lean_bin,
            "exit_code": checker_rc,
            "stdout": checker_stdout,
            "stderr": checker_stderr,
        },
        "verdict": verdict,
        "required_next_step": (
            "PROVE_C5B_EXACT_DISTRIBUTIONAL_LIMIT_FROM_COMPLETE_FIN_ACTION"
            if verdict.startswith("L14_C5_DISCHARGE_SCAFFOLD_GATE_PASS")
            else "REPAIR_QW2170_AND_RERUN"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2170: L14 C5 DISCHARGE SCAFFOLD GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Boundary",
        "- C5 is decomposed into C5a (closed support) + C5b (terminal open exact limit).",
        "- No hidden closure: final continuum theorem remains open until C5b is proven.",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
