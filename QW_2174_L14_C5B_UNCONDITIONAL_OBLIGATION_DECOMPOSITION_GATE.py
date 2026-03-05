#!/usr/bin/env python3
"""
QW-2174: L14 C5b unconditional obligation decomposition gate.

Purpose:
- convert unconditional C5b boundary into explicit sub-obligations,
- close proxy-transfer layer, isolate single terminal exact-limit lemma,
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
OUT_JSON = ROOT / "report_qw2174_l14_c5b_unconditional_obligation_decomposition_gate.json"
OUT_MD = ROOT / "RAPORT_QW2174_L14_C5B_UNCONDITIONAL_OBLIGATION_DECOMPOSITION_GATE.md"
OUT_LEAN = ROOT / "FIN_L14_C5B_UNCONDITIONAL_DECOMPOSITION_QW2174.lean"
OUT_PACKET = ROOT / "proof_packet_qw2174_l14_c5b_unconditional_decomposition.json"


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
    r2172 = load("report_qw2172_l14_c5b_terminal_limit_reduction_gate.json")
    r2170 = load("report_qw2170_l14_c5_discharge_scaffold_gate.json")

    # V1: conditional proxy-transfer layer to C5b (closed)
    v1_closed = bool(r2172["flags"]["c5b_conditional_closed_under_explicit_bundle"])
    # V2: exact unconditional limit from complete action (terminal open)
    v2_closed = False
    c5b_closed = v1_closed and v2_closed

    obligations = [
        {
            "id": "V1",
            "name": "conditional_proxy_transfer_to_c5b_node",
            "status": "satisfied" if v1_closed else "open",
            "source": "QW-2172",
        },
        {
            "id": "V2",
            "name": "unconditional_exact_distribution_limit_from_complete_fin_action",
            "status": "open",
            "source": "terminal lemma",
        },
        {
            "id": "C5b",
            "name": "terminal_exact_distribution_limit_unconditional",
            "status": "open",
            "source": "V1 + V2",
        },
    ]
    nodes = [o["id"] for o in obligations]
    edges = [["V1", "C5b"], ["V2", "C5b"]]
    order = topo_sort(nodes, edges)
    acyclic = len(order) == len(nodes)

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L14 C5b unconditional decomposition (QW-2174)",
            "",
            "theorem l14_c5b_v1_v2_bundle",
            "  {V1 V2 : Prop}",
            "  (h1 : V1) (h2 : V2) : V1 ∧ V2 := by",
            "  exact And.intro h1 h2",
            "",
            "theorem l14_c5b_from_v1_v2",
            "  {V1 V2 C5b : Prop}",
            "  (hcomp : V1 -> V2 -> C5b)",
            "  (h1 : V1) (h2 : V2) : C5b := by",
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
        "packet_name": "QW2174_L14_C5B_UNCONDITIONAL_DECOMPOSITION",
        "inputs": {
            "q2172": "report_qw2172_l14_c5b_terminal_limit_reduction_gate.json",
            "q2170": "report_qw2170_l14_c5_discharge_scaffold_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "obligations": obligations,
        "dependency_graph": {"nodes": nodes, "edges": edges, "topological_order": order},
        "decomposition": {
            "v1_closed": v1_closed,
            "v2_closed": v2_closed,
            "c5b_closed": c5b_closed,
        },
        "hashes": {"lean_sha256": sha256_file(OUT_LEAN)},
    }
    OUT_PACKET.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "q2172_conditional_reduction_present": bool(r2172["flags"]["c5b_conditional_closed_under_explicit_bundle"]),
        "q2170_scaffold_present": bool(r2170["flags"]["c5_decomposition_declared"]),
        "v1_proxy_transfer_layer_closed": bool(v1_closed),
        "v2_terminal_unconditional_lemma_closed": bool(v2_closed),
        "decomposition_declared": True,
        "dependency_graph_acyclic": bool(acyclic),
        "single_terminal_unconditional_obligation_isolated": True,
        "no_placeholder_tokens_in_lean": bool(not placeholders),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "proof_packet_manifest_written": OUT_PACKET.exists(),
        "terminal_c5b_exact_distribution_limit_closed": bool(c5b_closed),
        "full_continuum_theorem_from_complete_fin_action_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L14_C5B_UNCONDITIONAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OPEN"
        if (
            flags["q2172_conditional_reduction_present"]
            and flags["q2170_scaffold_present"]
            and flags["v1_proxy_transfer_layer_closed"]
            and flags["decomposition_declared"]
            and flags["dependency_graph_acyclic"]
            and flags["single_terminal_unconditional_obligation_isolated"]
            and flags["no_placeholder_tokens_in_lean"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["proof_packet_manifest_written"]
        )
        else "L14_C5B_UNCONDITIONAL_OBLIGATION_DECOMPOSITION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2172": "report_qw2172_l14_c5b_terminal_limit_reduction_gate.json",
            "q2170": "report_qw2170_l14_c5_discharge_scaffold_gate.json",
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
        "required_next_step": "PROVE_V2_UNCONDITIONAL_C5B_FROM_COMPLETE_FIN_ACTION",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2174: L14 C5B UNCONDITIONAL OBLIGATION DECOMPOSITION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Boundary",
        "- Unconditional C5b is decomposed to V1 (closed) + V2 (single terminal open lemma).",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
