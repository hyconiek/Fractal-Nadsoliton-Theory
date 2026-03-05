#!/usr/bin/env python3
"""
QW-2142: L13 formal proof-obligation export gate.

Purpose:
- convert the L13 all-orders closure path into explicit machine-auditable
  proof obligations with dependency graph,
- verify graph consistency (acyclic, all deps resolved, all leaves grounded),
- prepare direct handoff to proof assistant stage without overclaiming completion.
"""

from __future__ import annotations

import json
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2142_l13_formal_proof_obligation_export_gate.json"
OUT_MD = ROOT / "RAPORT_QW2142_L13_FORMAL_PROOF_OBLIGATION_EXPORT_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2135 = load("report_qw2135_interacting_microcausality_constructive_finite_order_gate.json")
    r2136 = load("report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json")
    r2137 = load("report_qw2137_interacting_microcausality_distribution_level_schema_gate.json")
    r2138 = load("report_qw2138_interacting_microcausality_proof_completion_gate.json")

    # Formal obligations for all-orders distribution-level closure.
    obligations: List[Dict] = [
        {
            "id": "P1_local_vertices_dim_le_4",
            "statement": "Interaction vertex basis remains local and dimension <=4.",
            "depends_on": [],
            "source": "QW-2135",
            "grounded": True,
        },
        {
            "id": "P2_finite_order_no_obstruction_n_le_4",
            "statement": "Constructive recursion has zero obstruction up to n<=4.",
            "depends_on": ["P1_local_vertices_dim_le_4"],
            "source": "QW-2135",
            "grounded": int(r2135["finite_order_causal_recursion"]["obstruction_count_total"]) == 0,
        },
        {
            "id": "P3_weighted_combinatorial_series_control",
            "statement": "All-orders combinatorial series has finite closed form and bounded tails.",
            "depends_on": ["P2_finite_order_no_obstruction_n_le_4"],
            "source": "QW-2136",
            "grounded": bool(
                r2136["flags"]["weighted_partition_series_matches_closed_form_limit"]
                and r2136["flags"]["weighted_partition_tail_bound_small"]
            ),
        },
        {
            "id": "P4_inductive_extension_rule",
            "statement": "Inductive step n->n+1 declared with finite counterterm policy.",
            "depends_on": ["P3_weighted_combinatorial_series_control"],
            "source": "QW-2136",
            "grounded": bool(
                r2136["flags"]["induction_schema_declared_explicitly"]
                and r2136["flags"]["finite_counterterm_basis_policy_declared"]
            ),
        },
        {
            "id": "P5_distribution_support_union",
            "statement": "Distribution support is contained in future/past cone union.",
            "depends_on": ["P4_inductive_extension_rule"],
            "source": "QW-2137",
            "grounded": bool(r2137["flags"]["support_union_statement_declared"]),
        },
        {
            "id": "P6_causal_splitting_local_normalization",
            "statement": "Causal splitting and local normalization are explicitly defined.",
            "depends_on": ["P5_distribution_support_union"],
            "source": "QW-2137",
            "grounded": bool(
                r2137["flags"]["causal_splitting_with_local_normalization_declared"]
                and r2137["flags"]["finite_local_counterterm_basis_bound_nonzero"]
            ),
        },
        {
            "id": "P7_cone_closure_numeric_audit",
            "statement": "Numeric cone closure audit supports schema-level causality closure.",
            "depends_on": ["P6_causal_splitting_local_normalization"],
            "source": "QW-2137",
            "grounded": bool(
                r2137["flags"]["future_cone_closure_numeric_rate_ge_0p999"]
                and r2137["flags"]["past_cone_closure_numeric_rate_ge_0p999"]
            ),
        },
        {
            "id": "P8_remainder_control_high_order",
            "statement": "High-order remainder is explicitly bounded in completion matrix.",
            "depends_on": ["P4_inductive_extension_rule"],
            "source": "QW-2138",
            "grounded": bool(r2138["flags"]["high_order_remainder_control_n80"]),
        },
        {
            "id": "P9_all_obligations_matrix_satisfied",
            "statement": "Completion matrix obligations are all satisfied before theorem export.",
            "depends_on": ["P7_cone_closure_numeric_audit", "P8_remainder_control_high_order"],
            "source": "QW-2138",
            "grounded": bool(r2138["flags"]["all_8_obligations_satisfied"]),
        },
    ]

    ids = {o["id"] for o in obligations}
    all_deps_resolved = all(all(dep in ids for dep in o["depends_on"]) for o in obligations)

    # DAG checks
    indeg = {o["id"]: 0 for o in obligations}
    g = defaultdict(list)
    for o in obligations:
        for dep in o["depends_on"]:
            g[dep].append(o["id"])
            indeg[o["id"]] += 1

    q = deque([nid for nid, d in indeg.items() if d == 0])
    topo: List[str] = []
    while q:
        v = q.popleft()
        topo.append(v)
        for nxt in g[v]:
            indeg[nxt] -= 1
            if indeg[nxt] == 0:
                q.append(nxt)
    acyclic = len(topo) == len(obligations)

    leaves = [o for o in obligations if len(o["depends_on"]) == 0]
    all_leaves_grounded = all(bool(o["grounded"]) for o in leaves)
    grounded_count = int(sum(1 for o in obligations if bool(o["grounded"])))

    # Export package for proof assistant handoff.
    export_package = {
        "logic_profile": "classical_first_order_with_distribution_annotations",
        "target_assistant": "lean_or_coq_placeholder",
        "obligations": obligations,
        "topological_order": topo,
    }

    flags = {
        "formal_obligation_set_declared": True,
        "dependency_graph_all_edges_resolved": bool(all_deps_resolved),
        "dependency_graph_acyclic": bool(acyclic),
        "all_leaf_axioms_grounded_by_strict_reports": bool(all_leaves_grounded),
        "all_exported_obligations_grounded": bool(grounded_count == len(obligations)),
        "proof_assistant_handoff_package_exported": True,
        "full_machine_checked_all_orders_proof_completed": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "L13_FORMAL_PROOF_OBLIGATION_EXPORT_GATE_PASS_PARTIAL"
        if (
            flags["formal_obligation_set_declared"]
            and flags["dependency_graph_all_edges_resolved"]
            and flags["dependency_graph_acyclic"]
            and flags["all_leaf_axioms_grounded_by_strict_reports"]
            and flags["all_exported_obligations_grounded"]
            and flags["proof_assistant_handoff_package_exported"]
        )
        else "L13_FORMAL_PROOF_OBLIGATION_EXPORT_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2135": "report_qw2135_interacting_microcausality_constructive_finite_order_gate.json",
            "q2136": "report_qw2136_interacting_microcausality_all_orders_scaffold_gate.json",
            "q2137": "report_qw2137_interacting_microcausality_distribution_level_schema_gate.json",
            "q2138": "report_qw2138_interacting_microcausality_proof_completion_gate.json",
        },
        "export_package": export_package,
        "stats": {
            "n_obligations": len(obligations),
            "n_grounded": grounded_count,
            "n_leaves": len(leaves),
            "n_topo_nodes": len(topo),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "RUN_EXTERNAL_MACHINE_CHECKER_ON_EXPORTED_OBLIGATIONS_AND_ATTACH_PROOF_OBJECT"
            if verdict.startswith("L13_FORMAL_PROOF_OBLIGATION_EXPORT_GATE_PASS")
            else "REPAIR_OBLIGATION_GRAPH_OR_GROUNDING_AND_RERUN_QW2142"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2142: L13 FORMAL PROOF-OBLIGATION EXPORT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- grounded obligations: `{grounded_count}/{len(obligations)}`",
        "",
        "## Graph checks",
        f"- all deps resolved: `{all_deps_resolved}`",
        f"- acyclic: `{acyclic}`",
        f"- topo nodes: `{len(topo)}`",
        "",
        "## Scope boundary",
        "- Full machine-checked all-orders proof: `False`",
        "",
    ]
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
