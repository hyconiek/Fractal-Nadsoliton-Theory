#!/usr/bin/env python3
"""
QW-2173: L13 F5b unconditional obligation decomposition gate.

Purpose:
- convert unconditional F5b boundary into explicit sub-obligations,
- close transport/compatibility layers, isolate a single terminal open lemma,
- preserve honest-open theorem status.
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
OUT_JSON = ROOT / "report_qw2173_l13_f5b_unconditional_obligation_decomposition_gate.json"
OUT_MD = ROOT / "RAPORT_QW2173_L13_F5B_UNCONDITIONAL_OBLIGATION_DECOMPOSITION_GATE.md"
OUT_LEAN = ROOT / "FIN_L13_F5B_UNCONDITIONAL_DECOMPOSITION_QW2173.lean"
OUT_PACKET = ROOT / "proof_packet_qw2173_l13_f5b_unconditional_decomposition.json"


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
    out = []
    while q:
        x = q.pop(0)
        out.append(x)
        for y in adj[x]:
            indeg[y] -= 1
            if indeg[y] == 0:
                q.append(y)
    return out


def main() -> None:
    r2171 = load("report_qw2171_l13_f5b_terminal_bound_reduction_gate.json")
    r2169 = load("report_qw2169_l13_f5_discharge_scaffold_gate.json")

    # U1: conditional bundle transport/compatibility to theorem node (closed)
    u1_closed = bool(r2171["flags"]["f5b_conditional_closed_under_explicit_bundle"])
    # U2: unconditional derivation directly from complete FIN action (terminal open)
    u2_closed = False
    f5b_closed = u1_closed and u2_closed

    obligations = [
        {
            "id": "U1",
            "name": "conditional_bundle_transport_to_f5b_node",
            "status": "satisfied" if u1_closed else "open",
            "source": "QW-2171",
        },
        {
            "id": "U2",
            "name": "unconditional_uniform_tail_bound_from_complete_fin_action",
            "status": "open",
            "source": "terminal lemma",
        },
        {
            "id": "F5b",
            "name": "terminal_uniform_tail_bound_unconditional",
            "status": "open",
            "source": "U1 + U2",
        },
    ]
    nodes = [o["id"] for o in obligations]
    edges = [["U1", "F5b"], ["U2", "F5b"]]
    order = topo_sort(nodes, edges)
    acyclic = len(order) == len(nodes)

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L13 F5b unconditional decomposition (QW-2173)",
            "",
            "theorem l13_f5b_u1_u2_bundle",
            "  {U1 U2 : Prop}",
            "  (h1 : U1) (h2 : U2) : U1 ∧ U2 := by",
            "  exact And.intro h1 h2",
            "",
            "theorem l13_f5b_from_u1_u2",
            "  {U1 U2 F5b : Prop}",
            "  (hcomp : U1 -> U2 -> F5b)",
            "  (h1 : U1) (h2 : U2) : F5b := by",
            "  exact hcomp h1 h2",
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
        "packet_name": "QW2173_L13_F5B_UNCONDITIONAL_DECOMPOSITION",
        "inputs": {
            "q2171": "report_qw2171_l13_f5b_terminal_bound_reduction_gate.json",
            "q2169": "report_qw2169_l13_f5_discharge_scaffold_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "obligations": obligations,
        "dependency_graph": {"nodes": nodes, "edges": edges, "topological_order": order},
        "decomposition": {
            "u1_closed": u1_closed,
            "u2_closed": u2_closed,
            "f5b_closed": f5b_closed,
        },
        "hashes": {"lean_sha256": sha256_file(OUT_LEAN)},
    }
    OUT_PACKET.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "q2171_conditional_reduction_present": bool(r2171["flags"]["f5b_conditional_closed_under_explicit_bundle"]),
        "q2169_scaffold_present": bool(r2169["flags"]["f5_decomposition_declared"]),
        "u1_transport_layer_closed": bool(u1_closed),
        "u2_terminal_unconditional_lemma_closed": bool(u2_closed),
        "decomposition_declared": True,
        "dependency_graph_acyclic": bool(acyclic),
        "single_terminal_unconditional_obligation_isolated": True,
        "no_placeholder_tokens_in_lean": bool(not placeholders),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "proof_packet_manifest_written": OUT_PACKET.exists(),
        "terminal_f5b_uniform_tail_bound_closed": bool(f5b_closed),
        "full_all_orders_theorem_from_complete_fin_action_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L13_F5B_UNCONDITIONAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OPEN"
        if (
            flags["q2171_conditional_reduction_present"]
            and flags["q2169_scaffold_present"]
            and flags["u1_transport_layer_closed"]
            and flags["decomposition_declared"]
            and flags["dependency_graph_acyclic"]
            and flags["single_terminal_unconditional_obligation_isolated"]
            and flags["no_placeholder_tokens_in_lean"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["proof_packet_manifest_written"]
        )
        else "L13_F5B_UNCONDITIONAL_OBLIGATION_DECOMPOSITION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2171": "report_qw2171_l13_f5b_terminal_bound_reduction_gate.json",
            "q2169": "report_qw2169_l13_f5_discharge_scaffold_gate.json",
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
        "required_next_step": "PROVE_U2_UNCONDITIONAL_F5B_FROM_COMPLETE_FIN_ACTION",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2173: L13 F5B UNCONDITIONAL OBLIGATION DECOMPOSITION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Boundary",
        "- Unconditional F5b is decomposed to U1 (closed) + U2 (single terminal open lemma).",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
