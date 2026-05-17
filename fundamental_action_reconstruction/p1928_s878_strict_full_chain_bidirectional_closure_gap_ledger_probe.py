#!/usr/bin/env python3
"""P1928 S878 strict full-chain bidirectional closure gap ledger probe."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def chain_row(step: str, exported_by: str, state: str, strict_basis: str) -> dict:
    return {
        "step": step,
        "exported_by": exported_by,
        "state": state,
        "strict_basis": strict_basis,
    }


def main() -> None:
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1927 = load("p1927_s877_strict_c1_gr_counterexample_certification_probe.json")

    renorm_open_items = [
        "Explicit regulator-and-scheme complete one-loop/two-loop export for all active strict sectors",
        "Counterterm universality theorem over the shared strict operator basis (not only selected channels)",
    ]
    unitarity_open_items = [
        "Cutkosky completeness witness for all declared channels with explicit discontinuity equalities",
        "Optical-theorem closure theorem on the same common basis used in renormalization exports",
    ]
    background_open_items = [
        "Repaired background-independence witness that survives CE_curv_pair_v1 admissible branch",
        "Transport-invariance theorem FRW <-> Bianchi-I under strict premises with no proxy shortcuts",
    ]
    selector_open_items = [
        "QW-2191 discharge via explicit strict symmetry-breaking/selector premise or internal selector source theorem",
    ]

    out = {
        "packet_id": "P1928",
        "stage_id": "S878",
        "status": "OPEN_STRICT_CORE_CLOSURE_GAPS_WITH_FULL_CHAIN_LEDGER",
        "route": "strict_only",
        "strict_chain_forward_kernel_to_eom": [
            chain_row(
                "K_strict_gate(d) definition",
                "QW-2049 lineage",
                "EXPORTED_WORKING_KERNEL",
                "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            ),
            chain_row(
                "effective coefficients from strict kernel",
                "P1874/P1907 family",
                "PARTIAL_EXPORTED",
                "strict-effective coefficient dictionary active, not yet globally witness-closed",
            ),
            chain_row(
                "full non-skeleton L_total = L_SM + L_GR + interactions",
                "P1907",
                "EXPORTED_NON_SKELETON_REGISTRY",
                "term registry includes sector decomposition and coupling placements",
            ),
            chain_row(
                "covariant EOM and Einstein residual transport",
                "P1874/P1907/P1927",
                "EXPORTED_WITH_ACTIVE_BLOCKERS",
                "EOM map exported, but strict closure still blocked by unresolved witness classes",
            ),
        ],
        "strict_chain_reverse_eom_to_kernel": [
            chain_row(
                "EOM residual classes -> required coefficient equalities",
                "P1907/P1927",
                "PARTIAL_RECONSTRUCTION",
                "some channel-level equalities exported",
            ),
            chain_row(
                "coefficient equalities -> admissible strict parameter bands",
                "P1910-P1917 family",
                "PARTIAL_RECONSTRUCTION",
                "selected channels solved; global basis-wide theorem missing",
            ),
            chain_row(
                "parameter bands -> K_strict_gate admissibility",
                "QW-2049 + follow-up probes",
                "WORKING_NOT_CLOSED",
                "operational admissibility retained; strict-core uniqueness/selector still open",
            ),
        ],
        "full_lagrangian_non_skeleton_anchor": {
            "present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "source": "P1907",
            "note": "Current chain uses full Lagrangian route, not skeleton-only route.",
        },
        "counterexample_pressure_from_p1927": {
            "background_independence_branch": p1927.get("strict_core_closure_statusvector", {}).get("background_independence", "OPEN"),
            "impact": "Any final strict-core PASS is forbidden until CE_curv_pair_v1 branch is either excluded by strict premise or neutralized by a stronger witness theorem.",
        },
        "strict_core_missing_exports_witnesses_theorems": {
            "renormalization": renorm_open_items,
            "unitarity": unitarity_open_items,
            "background_independence": background_open_items,
            "selector_qw2191": selector_open_items,
        },
        "toe_gap_count": {
            "renormalization_open": len(renorm_open_items),
            "unitarity_open": len(unitarity_open_items),
            "background_independence_open": len(background_open_items),
            "selector_open": len(selector_open_items),
            "total_minimum_open_obligations": len(renorm_open_items)
            + len(unitarity_open_items)
            + len(background_open_items)
            + len(selector_open_items),
        },
        "strict_false_pass_guard": "No theorem-level strict-core closure claim is allowed while any item in toe_gap_count remains open.",
        "next_honest_step": "Export P1929: explicit background-independence repair candidate (or admissibility-exclusion theorem) tied to CE_curv_pair_v1, then re-evaluate the full strict statusvector.",
        "lay_explanation": "Mamy pełny tor od kernela strict do równań ruchu i z powrotem, ale nadal brakuje 7 twardych dowodów/eksportów. To nie porażka — to precyzyjna mapa braków, która mówi dokładnie co domknąć, żeby nie udawać sukcesu ToE.",
    }

    out_path = GEN / "p1928_s878_strict_full_chain_bidirectional_closure_gap_ledger_probe.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
