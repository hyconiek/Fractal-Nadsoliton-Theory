#!/usr/bin/env python3
from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P686 = (
    GENERATED
    / "p686_current_strict_t173_w_break_rooted_directed_state_full_transition_edge_compatibility_audit_probe.json"
)

OUT = (
    GENERATED
    / "p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P686.exists():
        artifact = {
            "stage": "P687",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P686.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p686 = load_json(IN_P686)
    edges_obj = (((p686.get("edge_audit") or {}).get("edges")) or {})
    if not isinstance(edges_obj, dict) or not edges_obj:
        artifact = {
            "stage": "P687",
            "status": "NOT_COMPUTABLE_INVALID_P686_EDGE_AUDIT_SHAPE",
            "as_of": AS_OF,
            "error": "P686 must contain edge_audit.edges dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    # Build an undirected sign graph from P686 edges:
    # s_ij = best_match_sign in {+1,-1} such that O_ij u_i ≈ s_ij u_j.
    nodes: set[str] = set()
    undirected_sign: dict[tuple[str, str], int] = {}
    edge_errors: list[dict[str, Any]] = []

    for edge_id, edge in edges_obj.items():
        if not isinstance(edge, dict):
            edge_errors.append({"edge": edge_id, "error": "edge entry not a dict"})
            continue
        src = edge.get("src")
        dst = edge.get("dst")
        s = edge.get("best_match_sign")
        ok = bool(edge.get("edge_compatible_up_to_sign") is True)
        if not ok:
            edge_errors.append({"edge": edge_id, "error": "edge not compatible up to sign"})
            continue
        if not isinstance(src, str) or not isinstance(dst, str):
            edge_errors.append({"edge": edge_id, "error": "missing src/dst"})
            continue
        if s not in (1, -1):
            edge_errors.append({"edge": edge_id, "error": f"best_match_sign not in {{+1,-1}} (got {s!r})"})
            continue

        a, b = (src, dst) if src <= dst else (dst, src)
        nodes.add(a)
        nodes.add(b)
        key = (a, b)
        if key in undirected_sign and undirected_sign[key] != int(s):
            edge_errors.append(
                {
                    "edge": edge_id,
                    "error": f"inconsistent duplicate sign for undirected edge {key}: existing={undirected_sign[key]}, new={int(s)}",
                }
            )
            continue
        undirected_sign[key] = int(s)

    if edge_errors:
        artifact = {
            "stage": "P687",
            "status": "NOT_COMPUTABLE_P686_EDGE_ERRORS_PRESENT",
            "as_of": AS_OF,
            "edge_errors": edge_errors,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    # Solve for per-chart signs t_i in {±1} such that for every overlap edge (i,j):
    # s_ij = t_i * t_j  (equivalently t_j = s_ij * t_i).
    assigned: dict[str, int] = {}
    conflicts: list[dict[str, Any]] = []

    # adjacency for BFS: node -> list of (neighbor, s_ij)
    adj: dict[str, list[tuple[str, int]]] = {n: [] for n in sorted(nodes)}
    for (a, b), s in undirected_sign.items():
        adj[a].append((b, s))
        adj[b].append((a, s))

    for root in sorted(nodes):
        if root in assigned:
            continue
        assigned[root] = 1
        queue = [root]
        while queue:
            u = queue.pop(0)
            t_u = assigned[u]
            for v, s_uv in adj.get(u, []):
                expected_t_v = int(s_uv) * int(t_u)
                if v not in assigned:
                    assigned[v] = expected_t_v
                    queue.append(v)
                    continue
                if int(assigned[v]) != expected_t_v:
                    # Found inconsistency: s_uv != t_u * t_v.
                    conflicts.append(
                        {
                            "u": u,
                            "v": v,
                            "s_uv": int(s_uv),
                            "t_u": int(t_u),
                            "t_v": int(assigned[v]),
                            "expected_t_v": expected_t_v,
                            "witness_equation": "s_uv != t_u * t_v",
                        }
                    )

    sign_system_solvable = len(conflicts) == 0

    # Provide one explicit triangle witness with negative sign product if obstructed.
    triangle_witness = None
    if not sign_system_solvable:
        for a, b, c in itertools.combinations(sorted(nodes), 3):
            s_ab = undirected_sign.get((a, b) if a <= b else (b, a))
            s_bc = undirected_sign.get((b, c) if b <= c else (c, b))
            s_ac = undirected_sign.get((a, c) if a <= c else (c, a))
            if s_ab is None or s_bc is None or s_ac is None:
                continue
            prod = int(s_ab) * int(s_bc) * int(s_ac)
            if prod == -1:
                triangle_witness = {
                    "nodes": [a, b, c],
                    "edge_signs": {
                        f"{a}-{b}": int(s_ab),
                        f"{b}-{c}": int(s_bc),
                        f"{a}-{c}": int(s_ac),
                    },
                    "product": -1,
                    "meaning": "Negative Z2 holonomy witness on a 3-cycle: no global per-chart sign relift can make all overlap edges sign-consistent under fixed axis-only transition representatives.",
                }
                break

    status = "PASS_SIGN_COHERENCE_SOLVABLE" if sign_system_solvable else "PASS_SIGN_COHERENCE_OBSTRUCTED"

    artifact = {
        "stage": "P687",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "decide_whether_edgewise_sign_flips_from_P686_can_be_removed_by_a_per_chart_Z2_sign_relift_keeping_exported_transition_operators_fixed__no_false_pass",
        "inputs": {
            "p686_full_edge_audit_ref": str(IN_P686.relative_to(REPO)),
        },
        "status": status,
        "sign_graph": {
            "nodes": sorted(nodes),
            "edge_count": len(undirected_sign),
            "edge_signs_undirected": {f"{a}-{b}": s for (a, b), s in sorted(undirected_sign.items())},
        },
        "result": {
            "sign_system_solvable": sign_system_solvable,
            "chart_sign_solution_if_solvable": assigned if sign_system_solvable else None,
            "conflicts": conflicts,
            "triangle_witness_if_obstructed": triangle_witness,
        },
        "counts_as_strict_physical_orientation_datum": False,
        "hard_limits": [
            "Does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P687",
        "status": status,
        "sign_system_solvable": sign_system_solvable,
        "edge_count": len(undirected_sign),
        "triangle_witness_present": bool(triangle_witness is not None),
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

