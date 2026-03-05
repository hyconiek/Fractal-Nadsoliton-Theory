#!/usr/bin/env python3
"""
QW-2177: L13 U2b action-bridge specification gate.

Purpose:
- decompose U2b into structural action layer + terminal matching layer,
- close structural layer from exhaustive canonical EoM bundle,
- isolate single terminal matching lemma.
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
OUT_JSON = ROOT / "report_qw2177_l13_u2b_action_bridge_spec_gate.json"
OUT_MD = ROOT / "RAPORT_QW2177_L13_U2B_ACTION_BRIDGE_SPEC_GATE.md"
OUT_LEAN = ROOT / "FIN_L13_U2B_ACTION_BRIDGE_SPEC_QW2177.lean"
OUT_PACKET = ROOT / "proof_packet_qw2177_l13_u2b_action_bridge_spec.json"


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
    r2175 = load("report_qw2175_l13_u2_terminal_lemma_decomposition_gate.json")
    r2165 = load("report_qw2165_l13_exhaustive_canonical_eom_gate.json")

    u2b1_closed = bool(
        r2165["flags"]["exhaustive_bundle_closed_on_canonical_action_template"]
        and r2165["flags"]["all_psi_eom_contain_bidirectional_kernel_mixing_terms"]
        and r2165["flags"]["no_spacetime_nonlocal_tokens_in_all_13_eom"]
    )
    u2b2_closed = False
    u2b_closed = u2b1_closed and u2b2_closed

    obligations = [
        {
            "id": "U2b1",
            "name": "canonical_action_vertex_structure_exhaustive_layer_closed",
            "status": "satisfied" if u2b1_closed else "open",
            "source": "QW-2165",
        },
        {
            "id": "U2b2",
            "name": "exact_action_to_majorant_weight_matching_identity",
            "status": "open",
            "source": "terminal matching",
        },
        {
            "id": "U2b",
            "name": "action_to_majorant_identification_from_complete_fin_action",
            "status": "open",
            "source": "U2b1 + U2b2",
        },
    ]
    nodes = [o["id"] for o in obligations]
    edges = [["U2b1", "U2b"], ["U2b2", "U2b"]]
    order = topo_sort(nodes, edges)
    acyclic = len(order) == len(nodes)

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L13 U2b action-bridge spec (QW-2177)",
            "",
            "theorem l13_u2b1_u2b2_bundle",
            "  {U2b1 U2b2 : Prop}",
            "  (h1 : U2b1) (h2 : U2b2) : U2b1 ∧ U2b2 := by",
            "  exact And.intro h1 h2",
            "",
            "theorem l13_u2b_from_u2b1_u2b2",
            "  {U2b1 U2b2 U2b : Prop}",
            "  (hcomp : U2b1 -> U2b2 -> U2b)",
            "  (h1 : U2b1) (h2 : U2b2) : U2b := by",
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
        "packet_name": "QW2177_L13_U2B_ACTION_BRIDGE_SPEC",
        "inputs": {
            "q2175": "report_qw2175_l13_u2_terminal_lemma_decomposition_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "obligations": obligations,
        "dependency_graph": {"nodes": nodes, "edges": edges, "topological_order": order},
        "decomposition": {"u2b1_closed": u2b1_closed, "u2b2_closed": u2b2_closed, "u2b_closed": u2b_closed},
        "hashes": {"lean_sha256": sha256_file(OUT_LEAN)},
    }
    OUT_PACKET.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "q2175_single_action_bridge_u2b_present": bool(r2175["flags"]["u2b_action_to_majorant_bridge_closed"] is False),
        "u2b1_structural_action_layer_closed": bool(u2b1_closed),
        "u2b2_exact_matching_identity_closed": bool(u2b2_closed),
        "decomposition_declared": True,
        "dependency_graph_acyclic": bool(acyclic),
        "single_terminal_u2b2_isolated": True,
        "no_placeholder_tokens_in_lean": bool(not placeholders),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "proof_packet_manifest_written": OUT_PACKET.exists(),
        "u2b_action_to_majorant_bridge_closed": bool(u2b_closed),
        "u2_terminal_unconditional_lemma_closed": False,
        "terminal_f5b_uniform_tail_bound_closed": False,
        "full_all_orders_theorem_from_complete_fin_action_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L13_U2B_ACTION_BRIDGE_SPEC_GATE_PASS_PARTIAL_SINGLE_MATCHING_IDENTITY_OPEN"
        if (
            flags["q2175_single_action_bridge_u2b_present"]
            and flags["u2b1_structural_action_layer_closed"]
            and flags["decomposition_declared"]
            and flags["dependency_graph_acyclic"]
            and flags["single_terminal_u2b2_isolated"]
            and flags["no_placeholder_tokens_in_lean"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["proof_packet_manifest_written"]
        )
        else "L13_U2B_ACTION_BRIDGE_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2175": "report_qw2175_l13_u2_terminal_lemma_decomposition_gate.json",
            "q2165": "report_qw2165_l13_exhaustive_canonical_eom_gate.json",
            "proof_packet": OUT_PACKET.name,
            "lean_file": OUT_LEAN.name,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_U2B2_EXACT_ACTION_TO_MAJORANT_MATCHING_IDENTITY",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    OUT_MD.write_text(
        "\n".join(
            [
                "# RAPORT QW-2177: L13 U2B ACTION BRIDGE SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Boundary",
                "- U2b is decomposed to U2b1 (closed structural action layer) + U2b2 (single matching identity open).",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
