#!/usr/bin/env python3
"""
QW-2169: L13 F5 discharge scaffold gate.

Purpose:
- decompose final open F5 step into explicit sub-obligations,
- show finite/order-supported part is closed in strict chain,
- keep explicit terminal open boundary on uniform all-orders tail bound from
  complete FIN action.
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
OUT_JSON = ROOT / "report_qw2169_l13_f5_discharge_scaffold_gate.json"
OUT_MD = ROOT / "RAPORT_QW2169_L13_F5_DISCHARGE_SCAFFOLD_GATE.md"
OUT_LEAN = ROOT / "FIN_L13_F5_DISCHARGE_SCAFFOLD_QW2169.lean"
OUT_PACKET = ROOT / "proof_packet_qw2169_l13_f5_discharge_scaffold.json"


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
    r2167 = load("report_qw2167_l13_final_all_orders_theorem_packet_gate.json")
    r2165 = load("report_qw2165_l13_exhaustive_canonical_eom_gate.json")
    r2138 = load("report_qw2138_interacting_microcausality_proof_completion_gate.json")
    r2136 = load("report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json")

    f5a_closed = bool(
        r2165["flags"]["exhaustive_bundle_closed_on_canonical_action_template"]
        and r2138["flags"]["all_8_obligations_satisfied"]
        and r2138["flags"]["high_order_remainder_control_n80"]
        and r2136["flags"]["induction_schema_declared_explicitly"]
        and r2136["flags"]["finite_counterterm_basis_condition_holds"]
    )

    # Terminal theoretical boundary (honest-open).
    f5b_closed = False
    f5_closed = f5a_closed and f5b_closed

    obligations = [
        {
            "id": "F5a",
            "name": "finite_order_plus_induction_scaffold_support_from_strict_chain",
            "status": "satisfied" if f5a_closed else "open",
            "source": "QW-2165/QW-2138/QW-2136",
        },
        {
            "id": "F5b",
            "name": "uniform_all_orders_tail_bound_directly_from_complete_fin_action",
            "status": "open",
            "source": "terminal theorem",
        },
        {
            "id": "F5",
            "name": "all_orders_variational_completion_from_complete_fin_action",
            "status": "open",
            "source": "F5a + F5b",
        },
    ]

    nodes = [o["id"] for o in obligations]
    edges = [["F5a", "F5"], ["F5b", "F5"]]
    topo = topo_sort(nodes, edges)
    acyclic = len(topo) == len(nodes)

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L13 F5 discharge scaffold (QW-2169)",
            "",
            "theorem l13_f5_discharge_scaffold",
            "  {F5a F5b : Prop}",
            "  (ha : F5a) (hb : F5b) :",
            "  F5a ∧ F5b := by",
            "  exact And.intro ha hb",
            "",
            "theorem l13_f5_composition",
            "  {F5a F5b F5 : Prop}",
            "  (hcomp : F5a -> F5b -> F5)",
            "  (ha : F5a) (hb : F5b) : F5 := by",
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
        "packet_name": "QW2169_L13_F5_DISCHARGE_SCAFFOLD",
        "inputs": {
            "q2167": "report_qw2167_l13_final_all_orders_theorem_packet_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "q2138": "report_qw2138_interacting_microcausality_proof_completion_gate.json",
            "q2136": "report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "obligations": obligations,
        "dependency_graph": {"nodes": nodes, "edges": edges, "topological_order": topo},
        "decomposition": {
            "f5a_finite_induction_support_closed": f5a_closed,
            "f5b_uniform_tail_bound_closed": f5b_closed,
            "f5_closed": f5_closed,
        },
        "hashes": {"lean_sha256": sha256_file(OUT_LEAN)},
    }
    OUT_PACKET.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "q2167_final_packet_present": bool(r2167["flags"]["final_obligation_vector_declared"]),
        "q2165_exhaustive_canonical_eom_present": bool(r2165["flags"]["exhaustive_bundle_closed_on_canonical_action_template"]),
        "f5a_finite_induction_support_closed": bool(f5a_closed),
        "f5_decomposition_declared": True,
        "f5_dependency_graph_acyclic": bool(acyclic),
        "f5_composition_schema_declared": True,
        "no_placeholder_tokens_in_scaffold_theorem": bool(not placeholders),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "proof_packet_manifest_written": OUT_PACKET.exists(),
        "terminal_f5b_uniform_tail_bound_closed": False,
        "full_all_orders_theorem_from_complete_fin_action_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L13_F5_DISCHARGE_SCAFFOLD_GATE_PASS_PARTIAL_TERMINAL_BOUND_OPEN"
        if (
            flags["q2167_final_packet_present"]
            and flags["q2165_exhaustive_canonical_eom_present"]
            and flags["f5a_finite_induction_support_closed"]
            and flags["f5_decomposition_declared"]
            and flags["f5_dependency_graph_acyclic"]
            and flags["f5_composition_schema_declared"]
            and flags["no_placeholder_tokens_in_scaffold_theorem"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["proof_packet_manifest_written"]
        )
        else "L13_F5_DISCHARGE_SCAFFOLD_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2167": "report_qw2167_l13_final_all_orders_theorem_packet_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "q2138": "report_qw2138_interacting_microcausality_proof_completion_gate.json",
            "q2136": "report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json",
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
            "PROVE_F5B_UNIFORM_TAIL_BOUND_FROM_COMPLETE_FIN_ACTION"
            if verdict.startswith("L13_F5_DISCHARGE_SCAFFOLD_GATE_PASS")
            else "REPAIR_QW2169_AND_RERUN"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2169: L13 F5 DISCHARGE SCAFFOLD GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Boundary",
        "- F5 is decomposed into F5a (closed support) + F5b (terminal open uniform bound).",
        "- No hidden closure: final all-orders theorem remains open until F5b is proven.",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
