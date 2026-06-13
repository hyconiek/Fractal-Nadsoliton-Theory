#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.json"
MD = GEN / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.md"
P2663 = GEN / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NODES = ["nadsoliton", "light", "matter", "observer", "compression"]
EDGES = [
    ("nadsoliton", "light"),
    ("light", "observer"),
    ("observer", "matter"),
    ("matter", "nadsoliton"),
    ("light", "compression"),
    ("compression", "observer"),
]
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
CYCLE_SQUARE = [("nadsoliton", "light"), ("light", "observer"), ("observer", "matter"), ("matter", "nadsoliton")]
CYCLE_TRIANGLE = [("light", "observer"), ("light", "compression"), ("compression", "observer")]
LOCAL_COEFFICIENTS = [0.25, 1.0, 3.0]
THETA_COEFFICIENTS = [-2.0, -0.5, 0.5, 2.0]
NEGATIVE_EXPORT_FLAGS = [
    "nonexact_sector_selector_exported_unconditionally",
    "preferred_cycle_functional_exported_as_dynamics",
    "boundary_phase_bit_target_exported_unconditionally",
    "intrinsic_entropy_level_exported",
    "bit_to_action_map_sourced_unconditionally",
    "bit_to_length_map_sourced_unconditionally",
    "unique_beta_representative_selected_unconditionally",
    "entropy_arrow_discharges_qw_2191",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
    "bridge_completion_exported",
    "role_transfer_revalidated",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, ".",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
        "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
    ], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "sector_selector_content": "sector selector|non-exact sector|holonomy sector|preferred cycle",
        "variational_phase_content": "variational|action functional|theta|phase source|boundary-phase",
        "entropy_bit_target_content": "N_bits|bit target|log\\(2\\)|entropy target",
        "nonclosure_guard_content": "QW-2191|role-bearing L_total|bridge completion|role transfer|ToE closure|beta source",
    }
    return {"tool": "rg", "mode": "content-first boundary-phase sector-selector no-go audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def edge_value(cochain: tuple[int, ...], edge: tuple[str, str]) -> int:
    if edge in EDGE_INDEX:
        return cochain[EDGE_INDEX[edge]]
    rev = (edge[1], edge[0])
    return cochain[EDGE_INDEX[rev]]


def holonomy(cochain: tuple[int, ...], cycle: list[tuple[str, str]]) -> int:
    return sum(edge_value(cochain, edge) for edge in cycle) % 2


def triangle_flat(cochain: tuple[int, ...]) -> bool:
    return holonomy(cochain, CYCLE_TRIANGLE) == 0


def coboundary(vertex_bits: dict[str, int]) -> tuple[int, ...]:
    return tuple((vertex_bits[a] + vertex_bits[b]) % 2 for a, b in EDGES)


def flat_cochains() -> list[tuple[int, ...]]:
    return [c for c in itertools.product([0, 1], repeat=len(EDGES)) if triangle_flat(c)]


def gauge_orbit(cochain: tuple[int, ...]) -> set[tuple[int, ...]]:
    orbit = set()
    for bits in itertools.product([0, 1], repeat=len(NODES)):
        delta = coboundary(dict(zip(NODES, bits, strict=True)))
        orbit.add(tuple((x + y) % 2 for x, y in zip(cochain, delta, strict=True)))
    return orbit


def local_even_variational_energy(cochain: tuple[int, ...], flatness_weight: float, edge_weight: float) -> float:
    # The audited no-go class includes only positive local penalties visible before a sector source is added.
    # It can enforce filled-triangle flatness and penalize edge excitations, but it contains no linear term in the nontrivial square holonomy.
    triangle_penalty = flatness_weight * holonomy(cochain, CYCLE_TRIANGLE)
    edge_penalty = edge_weight * sum(cochain) / len(EDGES)
    return triangle_penalty + edge_penalty


def theta_energy(cochain: tuple[int, ...], theta: float) -> float:
    # This term is deliberately classified as declared: it chooses the non-exact square sector by putting in a theta/source coupling.
    return -theta * holonomy(cochain, CYCLE_SQUARE)


def sector_selector_witness() -> dict[str, Any]:
    flats = flat_cochains()
    exact = {coboundary(dict(zip(NODES, bits, strict=True))) for bits in itertools.product([0, 1], repeat=len(NODES))}
    sector_counts = {0: 0, 1: 0}
    exact_sector_counts = {0: 0, 1: 0}
    for cochain in flats:
        sector_counts[holonomy(cochain, CYCLE_SQUARE)] += 1
        if cochain in exact:
            exact_sector_counts[holonomy(cochain, CYCLE_SQUARE)] += 1

    gauge_representatives = []
    unseen = set(flats)
    while unseen:
        seed = min(unseen)
        orbit = gauge_orbit(seed) & set(flats)
        gauge_representatives.append({
            "representative": list(seed),
            "orbit_size_inside_flat_slice": len(orbit),
            "square_holonomy_sector": holonomy(seed, CYCLE_SQUARE),
            "contains_exact_coboundary": any(c in exact for c in orbit),
        })
        unseen -= orbit

    local_rows = []
    for flatness_weight in LOCAL_COEFFICIENTS:
        for edge_weight in LOCAL_COEFFICIENTS:
            energies = [(local_even_variational_energy(c, flatness_weight, edge_weight), c, holonomy(c, CYCLE_SQUARE)) for c in flats]
            minimum = min(e for e, _, _ in energies)
            minimizers = [(c, sector) for e, c, sector in energies if abs(e - minimum) < 1e-12]
            sectors = sorted({sector for _, sector in minimizers})
            local_rows.append({
                "flatness_weight": flatness_weight,
                "edge_weight": edge_weight,
                "minimum_energy": minimum,
                "minimizer_count": len(minimizers),
                "minimizer_square_holonomy_sectors": sectors,
                "selects_nonexact_one_bit_sector": sectors == [1],
            })

    theta_rows = []
    for theta in THETA_COEFFICIENTS:
        energies = [(theta_energy(c, theta), c, holonomy(c, CYCLE_SQUARE)) for c in flats]
        minimum = min(e for e, _, _ in energies)
        sectors = sorted({sector for e, _, sector in energies if abs(e - minimum) < 1e-12})
        theta_rows.append({
            "declared_theta": theta,
            "minimizer_square_holonomy_sectors": sectors,
            "selects_nonexact_one_bit_sector": sectors == [1],
            "is_internal_source": False,
        })

    return {
        "statement": "P2664 audits the next P2663 gap: can an internal variational boundary-phase rule select the non-exact one-bit square sector?  The finite witness separates three facts: the flat complex has two holonomy sectors, exact/gauge coboundaries remain in the zero sector, and positive local flatness/edge penalties do not select the non-exact one-bit sector.  A theta-like holonomy source can select it, but only by declaring the missing sector source.",
        "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
        "no_sub_nadsoliton_information_layer": True,
        "triangle_flat_cochain_count": len(flats),
        "square_holonomy_sector_counts_after_flatness": sector_counts,
        "exact_coboundary_sector_counts": exact_sector_counts,
        "gauge_orbit_representatives": gauge_representatives,
        "gauge_orbit_count_inside_flat_slice": len(gauge_representatives),
        "local_even_variational_rows": local_rows,
        "theta_source_rows": theta_rows,
        "local_even_class_selects_nonexact_sector": any(row["selects_nonexact_one_bit_sector"] for row in local_rows),
        "declared_theta_can_select_nonexact_sector": any(row["selects_nonexact_one_bit_sector"] for row in theta_rows),
        "theta_selection_is_internal_source": False,
    }


def upstream_consistency() -> dict[str, Any]:
    p2663 = load_json(P2663)
    return {
        "p2663_nonexact_bit_carrier_found": p2663.get("closure_decision", {}).get("nonexact_bit_carrier_found") is True,
        "p2663_unique_target_not_derived": p2663.get("closure_decision", {}).get("unique_n_bits_derived_without_sector_choice") is False,
        "p2663_boundary_phase_bit_target_not_exported": p2663.get("closure_decision", {}).get("boundary_phase_bit_target_exported_now") is False,
        "p2663_qw2191_not_discharged": p2663.get("closure_decision", {}).get("qw2191_discharged_now") is False,
    }


def source_candidate_matrix(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"candidate": "positive_local_flatness_edge_action", "computational_result": "does not minimize uniquely in non-exact sector", "passes_as_sector_selector_now": False, "verdict": "blocked: local positive penalties prefer consistency/small support, not the one-bit holonomy sector"},
        {"candidate": "gauge_exact_boundary_phase_dynamics", "computational_result": f"exact sector counts {witness['exact_coboundary_sector_counts']}", "passes_as_sector_selector_now": False, "verdict": "blocked: exact coboundary dynamics stays in zero square holonomy"},
        {"candidate": "declared_theta_holonomy_source", "computational_result": "can select sector for chosen sign", "passes_as_sector_selector_now": False, "verdict": "false pass: theta/holonomy coupling is precisely the missing source unless derived from nadsoliton dynamics"},
        {"candidate": "bridge_completed_sector_selector_theorem_target", "computational_result": "not present in current audit inputs", "passes_as_sector_selector_now": False, "verdict": "open theorem target: derive the theta sign or non-exact sector from bridge-completed nadsoliton dynamics"},
    ]


def closure_decision(witness: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_sector_selector_now"]]
    return {
        "decision": "P2664_BOUNDARY_PHASE_SECTOR_SELECTOR_VARIATIONAL_NO_GO_AUDIT__CARRIER_REMAINS_UNSOURCED",
        "professorial_verdict": "P2664 is the requested more proof-grade/computational continuation of P2663.  It shows that the non-exact one-bit holonomy carrier is not enough: the audited internal positive local boundary-phase action class does not select the non-exact sector, and gauge-exact dynamics remains zero-holonomy.  A theta-like holonomy source would select a sector for a chosen sign, but that sign/source is exactly the missing premise.  Therefore the current repo still has no boundary-phase bit-target source, no intrinsic UV unit, no beta source, no QW-2191 discharge, no L_total reopening, and no ToE closure.",
        "next_honest_step": "Search specifically for a bridge-completed nadsoliton term that derives a theta/holonomy-source sign or an equivalent non-exact sector selector.  The acceptance test is now sharp: the term must be gauge-invariant, not a declared entropy level, and must uniquely prefer square holonomy one over zero before P2662 is rerun.  If no such term exists, promote P2664 to a no-go for local positive boundary-phase variational selectors and keep the entropy/cocycle UV anchor conditional.",
        "passing_sector_selector_candidates": passing,
        "nonexact_bit_carrier_still_available": witness["square_holonomy_sector_counts_after_flatness"].get(1, 0) > 0,
        "local_even_class_selects_nonexact_sector": witness["local_even_class_selects_nonexact_sector"],
        "declared_theta_can_select_nonexact_sector": witness["declared_theta_can_select_nonexact_sector"],
        "theta_selection_is_internal_source": witness["theta_selection_is_internal_source"],
        "boundary_phase_sector_selector_exported_now": False,
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["sector_selector_witness"]
    d = payload["closure_decision"]
    lines = [
        "# P2664/S1614 boundary-phase sector-selector variational no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first audit",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += [
        "",
        "## Computational witness",
        w["statement"],
        f"Triangle-flat cochains: `{w['triangle_flat_cochain_count']}`.",
        f"Square holonomy sector counts after flatness: `{w['square_holonomy_sector_counts_after_flatness']}`.",
        f"Exact coboundary sector counts: `{w['exact_coboundary_sector_counts']}`.",
        f"Gauge orbit count inside flat slice: `{w['gauge_orbit_count_inside_flat_slice']}`.",
        f"Local even variational class selects non-exact sector? `{w['local_even_class_selects_nonexact_sector']}`.",
        f"Declared theta source can select non-exact sector? `{w['declared_theta_can_select_nonexact_sector']}`.",
        f"Theta selection is internal source? `{w['theta_selection_is_internal_source']}`.",
        "",
        "## Source candidate matrix",
        "| candidate | computational result | selector now? | verdict |",
        "| --- | --- | ---: | --- |",
    ]
    for row in payload["source_candidate_matrix"]:
        lines.append(f"| `{row['candidate']}` | {row['computational_result']} | `{row['passes_as_sector_selector_now']}` | {row['verdict']} |")
    lines += [
        "",
        "## Verdict",
        d["professorial_verdict"],
        f"Decision: `{d['decision']}`.",
        f"Passing sector-selector candidates: `{d['passing_sector_selector_candidates']}`.",
        f"Local even class selects non-exact sector? `{d['local_even_class_selects_nonexact_sector']}`.",
        f"Declared theta can select non-exact sector? `{d['declared_theta_can_select_nonexact_sector']}`.",
        f"Theta selection is internal source? `{d['theta_selection_is_internal_source']}`.",
        f"Boundary-phase bit target exported now? `{d['boundary_phase_bit_target_exported_now']}`.",
        f"Beta source exported now? `{d['beta_source_exported_now']}`.",
        f"QW-2191 discharged now? `{d['qw2191_discharged_now']}`.",
        f"Role-bearing L_total now? `{d['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{d['toe_closure_now']}`.",
        "",
        "## Next honest step",
        d["next_honest_step"],
        "",
        "## Negative exports",
    ]
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    witness = sector_selector_witness()
    matrix = source_candidate_matrix(witness)
    decision = closure_decision(witness, matrix)
    payload: dict[str, Any] = {
        "status": "P2664_BOUNDARY_PHASE_SECTOR_SELECTOR_VARIATIONAL_NO_GO_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {"P2663": sha256_file(P2663), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream_consistency(),
        "typed_chain_complex": {"nodes": NODES, "edges": EDGES, "square_cycle": CYCLE_SQUARE, "filled_triangle_cycle": CYCLE_TRIANGLE, "ontology_order": "nadsoliton -> light -> matter -> emergent observer", "no_sub_nadsoliton_information_layer": True},
        "sector_selector_witness": witness,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2664/S1614 boundary-phase sector-selector variational no-go guard", "## P2664/S1614 boundary-phase sector-selector variational no-go guard\n\n`P2664/S1614` audits whether the P2663 non-exact one-bit holonomy carrier is selected by an internal boundary-phase variational rule.  The audited positive local flatness/edge action class does not select the non-exact square sector, and gauge-exact dynamics remains zero-holonomy; a theta-like holonomy source can select a sector only by declaring the missing source/sign.  Therefore this is a sharper no-go for local positive boundary-phase selectors, not an intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2664/S1614 boundary-phase sector-selector Ltotal guard", "## P2664/S1614 boundary-phase sector-selector Ltotal guard\n\n`P2664/S1614` keeps `L_total` closed to a boundary-phase entropy source term.  The only audited way to force non-exact one-bit holonomy is a theta-like source/sign, and that coefficient is not internally derived.  A future Lagrangian term must first pass the P2664 acceptance test: gauge-invariant, non-declared, and uniquely selecting square holonomy one over zero from bridge-completed nadsoliton dynamics.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
