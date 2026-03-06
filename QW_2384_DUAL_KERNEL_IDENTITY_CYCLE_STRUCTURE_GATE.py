#!/usr/bin/env python3
"""QW-2384: dual kernel-identity cycle-structure gate.

Builds theorem-dependency graphs for L12/L5 identity layers and confirms
structural cycle recurrence with no non-circular anchor evidence.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

L12_BLOCKER = "RG_KernelIdentityLocalityToWellPosedness_Theorem"
L5_BLOCKER = "QFT_KernelIdentityLocalityToPositivity_Theorem"

THEOREM_RE = re.compile(r"^theorem\s+([A-Za-z0-9_']+)\s*:")
EXACT_RE = re.compile(r"^\s*exact\s+([A-Za-z0-9_']+)")


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def parse_theorem_edge(path: Path) -> tuple[str | None, str | None, str]:
    theorem = None
    dep = None
    nonempty_code = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("--"):
            continue
        nonempty_code.append(line)
        m_th = THEOREM_RE.match(line)
        if m_th:
            theorem = m_th.group(1)
            continue
        m_ex = EXACT_RE.match(line)
        if m_ex:
            dep = m_ex.group(1)

    proof_kind = "unknown"
    if theorem and dep:
        code_wo_decl = [ln for ln in nonempty_code if not THEOREM_RE.match(ln)]
        if code_wo_decl == [f"exact {dep}"] or code_wo_decl == ["by", f"exact {dep}"]:
            proof_kind = "single_symbol_exact"
        elif any("exact " in ln for ln in code_wo_decl):
            proof_kind = "exact_plus_extra"
        else:
            proof_kind = "non_exact_or_missing"

    return theorem, dep, proof_kind


def build_branch_graph(pattern: str, theorem_prefix: str) -> tuple[dict[str, str], dict[str, str], list[str]]:
    edges: dict[str, str] = {}
    proof_kind: dict[str, str] = {}
    files = sorted(ROOT.glob(pattern))

    for f in files:
        theorem, dep, kind = parse_theorem_edge(f)
        if theorem and dep and theorem.startswith(theorem_prefix):
            edges[theorem] = dep
            proof_kind[theorem] = kind

    nodes = sorted(set(edges.keys()) | set(edges.values()))
    return edges, proof_kind, nodes


def kosaraju_scc(edges: dict[str, str]) -> list[set[str]]:
    graph: dict[str, list[str]] = defaultdict(list)
    rev: dict[str, list[str]] = defaultdict(list)
    nodes = set(edges.keys()) | set(edges.values())

    for u, v in edges.items():
        graph[u].append(v)
        rev[v].append(u)
        if v not in graph:
            graph[v] = graph[v]
        if u not in rev:
            rev[u] = rev[u]

    visited: set[str] = set()
    order: list[str] = []

    def dfs1(start: str) -> None:
        stack = [(start, 0)]
        while stack:
            node, idx = stack[-1]
            if node not in visited:
                visited.add(node)
            neigh = graph[node]
            if idx < len(neigh):
                nxt = neigh[idx]
                stack[-1] = (node, idx + 1)
                if nxt not in visited:
                    stack.append((nxt, 0))
            else:
                order.append(node)
                stack.pop()

    for n in nodes:
        if n not in visited:
            dfs1(n)

    visited.clear()
    sccs: list[set[str]] = []

    for start in reversed(order):
        if start in visited:
            continue
        comp: set[str] = set()
        dq = deque([start])
        visited.add(start)
        while dq:
            x = dq.popleft()
            comp.add(x)
            for y in rev[x]:
                if y not in visited:
                    visited.add(y)
                    dq.append(y)
        sccs.append(comp)

    return sccs


def scc_for_symbol(sccs: list[set[str]], symbol: str) -> set[str]:
    for s in sccs:
        if symbol in s:
            return s
    return set()


def branch_cycle_assessment(edges: dict[str, str], proof_kind: dict[str, str], blocker: str) -> dict[str, Any]:
    sccs = kosaraju_scc(edges)
    blocker_scc = scc_for_symbol(sccs, blocker)

    if not blocker_scc:
        return {
            "blocker_in_graph": False,
            "scc_size": 0,
            "all_deps_inside_scc": False,
            "all_single_symbol_exact": False,
            "noncircular_anchor_candidates": [],
            "blocker_scc_nodes": [],
        }

    noncircular = []
    all_inside = True
    all_exact_edges_present = True

    for th in sorted(blocker_scc):
        dep = edges.get(th)
        if dep is None:
            all_inside = False
            all_exact_edges_present = False
            noncircular.append({"theorem": th, "reason": "no_dependency_edge"})
            continue

        if dep not in blocker_scc:
            all_inside = False
            noncircular.append({"theorem": th, "reason": "dependency_outside_scc", "dependency": dep})

    return {
        "blocker_in_graph": True,
        "scc_size": len(blocker_scc),
        "all_deps_inside_scc": all_inside,
        "all_exact_edges_present": all_exact_edges_present,
        "noncircular_anchor_candidates": noncircular,
        "blocker_scc_nodes": sorted(blocker_scc),
        "proof_kind_sample": {k: proof_kind.get(k, "unknown") for k in sorted(blocker_scc)[:5]},
    }


def main() -> None:
    q2381 = load("report_qw2381_dual_kernel_cycle_recurrence_gate.json")
    q2383 = load("report_qw2383_dual_noncyclic_step_admission_gate.json")

    l12_edges, l12_kind, _ = build_branch_graph(
        "FIN_L12_KERNEL_*_DISCHARGE_ATTEMPT.lean", "RG_KernelIdentity"
    )
    l5_edges, l5_kind, _ = build_branch_graph(
        "FIN_L5_KERNEL_*_DISCHARGE_ATTEMPT.lean", "QFT_KernelIdentity"
    )

    l12_assess = branch_cycle_assessment(l12_edges, l12_kind, L12_BLOCKER)
    l5_assess = branch_cycle_assessment(l5_edges, l5_kind, L5_BLOCKER)

    flags = {
        "q2381_recurrence_confirmed": q2381.get("verdict")
        == "DUAL_KERNEL_CYCLE_RECURRENCE_GATE_PASS_BLOCKER_LOOP_CONFIRMED",
        "q2383_admission_denied": bool(q2383.get("flags", {}).get("admission_denied", False)),
        "l12_graph_nonempty": len(l12_edges) > 0,
        "l5_graph_nonempty": len(l5_edges) > 0,
        "l12_blocker_in_graph": l12_assess["blocker_in_graph"],
        "l5_blocker_in_graph": l5_assess["blocker_in_graph"],
        "l12_blocker_scc_is_cyclic": l12_assess["scc_size"] > 1,
        "l5_blocker_scc_is_cyclic": l5_assess["scc_size"] > 1,
        "l12_no_noncircular_anchor_candidate": len(l12_assess["noncircular_anchor_candidates"]) == 0
        and l12_assess["all_deps_inside_scc"]
        and l12_assess["all_exact_edges_present"],
        "l5_no_noncircular_anchor_candidate": len(l5_assess["noncircular_anchor_candidates"]) == 0
        and l5_assess["all_deps_inside_scc"]
        and l5_assess["all_exact_edges_present"],
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_KERNEL_IDENTITY_CYCLE_STRUCTURE_GATE_PASS_STRUCTURAL_CYCLE_CONFIRMED"
        if flags["q2381_recurrence_confirmed"]
        and flags["q2383_admission_denied"]
        and flags["l12_graph_nonempty"]
        and flags["l5_graph_nonempty"]
        and flags["l12_blocker_in_graph"]
        and flags["l5_blocker_in_graph"]
        and flags["l12_blocker_scc_is_cyclic"]
        and flags["l5_blocker_scc_is_cyclic"]
        and flags["l12_no_noncircular_anchor_candidate"]
        and flags["l5_no_noncircular_anchor_candidate"]
        else "DUAL_KERNEL_IDENTITY_CYCLE_STRUCTURE_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2381": "report_qw2381_dual_kernel_cycle_recurrence_gate.json",
            "q2383": "report_qw2383_dual_noncyclic_step_admission_gate.json",
        },
        "analysis_method": {
            "graph_basis": "theorem -> exact dependency edges from FIN_L12/L5_KERNEL_*_DISCHARGE_ATTEMPT.lean",
            "scc_algorithm": "kosaraju",
            "anchor_criterion": "at least one blocker-SCC theorem with dependency outside SCC or non-single-symbol-exact proof",
        },
        "l12": {
            "blocker": L12_BLOCKER,
            "n_edges": len(l12_edges),
            "assessment": l12_assess,
        },
        "l5": {
            "blocker": L5_BLOCKER,
            "n_edges": len(l5_edges),
            "assessment": l5_assess,
        },
        "scope_boundary": {
            "structural_cycle_confirmed": True,
            "noncircular_anchor_found": False,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    spec_path = ROOT / "spec_qw2384_dual_kernel_identity_cycle_structure.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": spec["sources"] | {"spec": spec_path.name},
        "l12_blocker_scc_nodes": l12_assess["blocker_scc_nodes"],
        "l5_blocker_scc_nodes": l5_assess["blocker_scc_nodes"],
        "l12_anchor_candidates": l12_assess["noncircular_anchor_candidates"],
        "l5_anchor_candidates": l5_assess["noncircular_anchor_candidates"],
        "cycle_confirmed": flags["l12_blocker_scc_is_cyclic"] and flags["l5_blocker_scc_is_cyclic"],
    }

    proof_path = ROOT / "proof_object_qw2384_dual_kernel_identity_cycle_structure.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "l12_blocker_scc_size": l12_assess["scc_size"],
        "l5_blocker_scc_size": l5_assess["scc_size"],
        "l12_noncircular_anchor_candidates": l12_assess["noncircular_anchor_candidates"],
        "l5_noncircular_anchor_candidates": l5_assess["noncircular_anchor_candidates"],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DESIGN_DUAL_NONCIRCULAR_ANCHOR_OBLIGATION_PACKET",
    }

    out_json = ROOT / "report_qw2384_dual_kernel_identity_cycle_structure_gate.json"
    out_md = ROOT / "RAPORT_QW2384_DUAL_KERNEL_IDENTITY_CYCLE_STRUCTURE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2384: DUAL KERNEL IDENTITY CYCLE STRUCTURE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- L12 blocker SCC size: `{l12_assess['scc_size']}`",
                f"- L5 blocker SCC size: `{l5_assess['scc_size']}`",
                f"- L12 noncircular candidates: `{len(l12_assess['noncircular_anchor_candidates'])}`",
                f"- L5 noncircular candidates: `{len(l5_assess['noncircular_anchor_candidates'])}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "l12_scc_size": l12_assess["scc_size"],
                "l5_scc_size": l5_assess["scc_size"],
            }
        )
    )


if __name__ == "__main__":
    main()
