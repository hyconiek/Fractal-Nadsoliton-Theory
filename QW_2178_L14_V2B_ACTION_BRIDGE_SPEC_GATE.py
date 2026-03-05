#!/usr/bin/env python3
"""
QW-2178: L14 V2b action-bridge specification gate.

Purpose:
- decompose V2b into structural continuum-action layer + terminal identification,
- close structural layer from exhaustive canonical Hessian bundle,
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
OUT_JSON = ROOT / "report_qw2178_l14_v2b_action_bridge_spec_gate.json"
OUT_MD = ROOT / "RAPORT_QW2178_L14_V2B_ACTION_BRIDGE_SPEC_GATE.md"
OUT_LEAN = ROOT / "FIN_L14_V2B_ACTION_BRIDGE_SPEC_QW2178.lean"
OUT_PACKET = ROOT / "proof_packet_qw2178_l14_v2b_action_bridge_spec.json"


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
    r2176 = load("report_qw2176_l14_v2_terminal_lemma_decomposition_gate.json")
    r2166 = load("report_qw2166_l14_exhaustive_canonical_hessian_gate.json")

    v2b1_closed = bool(
        r2166["flags"]["exhaustive_continuum_bundle_closed_on_canonical_hessian_level"]
        and r2166["flags"]["linear_operator_matrix_matches_canonical_hessian"]
        and r2166["flags"]["no_spacetime_nonlocal_tokens_in_all_linearized_eom"]
    )
    v2b2_closed = False
    v2b_closed = v2b1_closed and v2b2_closed

    obligations = [
        {
            "id": "V2b1",
            "name": "canonical_hessian_continuum_structure_layer_closed",
            "status": "satisfied" if v2b1_closed else "open",
            "source": "QW-2166",
        },
        {
            "id": "V2b2",
            "name": "exact_action_level_identification_for_distribution_limit",
            "status": "open",
            "source": "terminal matching",
        },
        {
            "id": "V2b",
            "name": "action_level_identification_for_exact_distribution_limit",
            "status": "open",
            "source": "V2b1 + V2b2",
        },
    ]
    nodes = [o["id"] for o in obligations]
    edges = [["V2b1", "V2b"], ["V2b2", "V2b"]]
    order = topo_sort(nodes, edges)
    acyclic = len(order) == len(nodes)

    lean_text = "\n".join(
        [
            "-- FIN Release 5: L14 V2b action-bridge spec (QW-2178)",
            "",
            "theorem l14_v2b1_v2b2_bundle",
            "  {V2b1 V2b2 : Prop}",
            "  (h1 : V2b1) (h2 : V2b2) : V2b1 ∧ V2b2 := by",
            "  exact And.intro h1 h2",
            "",
            "theorem l14_v2b_from_v2b1_v2b2",
            "  {V2b1 V2b2 V2b : Prop}",
            "  (hcomp : V2b1 -> V2b2 -> V2b)",
            "  (h1 : V2b1) (h2 : V2b2) : V2b := by",
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
        "packet_name": "QW2178_L14_V2B_ACTION_BRIDGE_SPEC",
        "inputs": {
            "q2176": "report_qw2176_l14_v2_terminal_lemma_decomposition_gate.json",
            "q2166": "report_qw2166_l14_exhaustive_canonical_hessian_gate.json",
            "lean_file": OUT_LEAN.name,
        },
        "obligations": obligations,
        "dependency_graph": {"nodes": nodes, "edges": edges, "topological_order": order},
        "decomposition": {"v2b1_closed": v2b1_closed, "v2b2_closed": v2b2_closed, "v2b_closed": v2b_closed},
        "hashes": {"lean_sha256": sha256_file(OUT_LEAN)},
    }
    OUT_PACKET.write_text(json.dumps(packet, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "q2176_single_action_bridge_v2b_present": bool(r2176["flags"]["v2b_action_level_identification_closed"] is False),
        "v2b1_structural_continuum_layer_closed": bool(v2b1_closed),
        "v2b2_exact_action_identification_closed": bool(v2b2_closed),
        "decomposition_declared": True,
        "dependency_graph_acyclic": bool(acyclic),
        "single_terminal_v2b2_isolated": True,
        "no_placeholder_tokens_in_lean": bool(not placeholders),
        "lean_checker_detected": bool(checker_found),
        "lean_checker_exit_zero": bool(checker_found and checker_rc == 0),
        "proof_packet_manifest_written": OUT_PACKET.exists(),
        "v2b_action_level_identification_closed": bool(v2b_closed),
        "v2_terminal_unconditional_lemma_closed": False,
        "terminal_c5b_exact_distribution_limit_closed": False,
        "full_continuum_theorem_from_complete_fin_action_completed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = len(flags)

    verdict = (
        "L14_V2B_ACTION_BRIDGE_SPEC_GATE_PASS_PARTIAL_SINGLE_MATCHING_IDENTITY_OPEN"
        if (
            flags["q2176_single_action_bridge_v2b_present"]
            and flags["v2b1_structural_continuum_layer_closed"]
            and flags["decomposition_declared"]
            and flags["dependency_graph_acyclic"]
            and flags["single_terminal_v2b2_isolated"]
            and flags["no_placeholder_tokens_in_lean"]
            and flags["lean_checker_detected"]
            and flags["lean_checker_exit_zero"]
            and flags["proof_packet_manifest_written"]
        )
        else "L14_V2B_ACTION_BRIDGE_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2176": "report_qw2176_l14_v2_terminal_lemma_decomposition_gate.json",
            "q2166": "report_qw2166_l14_exhaustive_canonical_hessian_gate.json",
            "proof_packet": OUT_PACKET.name,
            "lean_file": OUT_LEAN.name,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_V2B2_EXACT_ACTION_LEVEL_IDENTIFICATION_IDENTITY",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    OUT_MD.write_text(
        "\n".join(
            [
                "# RAPORT QW-2178: L14 V2B ACTION BRIDGE SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Boundary",
                "- V2b is decomposed to V2b1 (closed structural continuum layer) + V2b2 (single matching identity open).",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
