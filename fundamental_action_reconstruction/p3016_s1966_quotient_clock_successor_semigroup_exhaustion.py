#!/usr/bin/env python3
"""P3016/S1966: quotient clock-successor semigroup exhaustion.

This follows P3015 by attacking exactly one missing successor/evolution atom:
can the U(12)-orbit quotient used by the unit-compatible time observable carry a
strict directed successor/semigroup that represents the label successor d -> d+1?

The answer is a finite exhaustive no-go.  A quotient successor S would have to
assign one target orbit to each source orbit while satisfying

    S([d]) = [d+1]

for all twelve labels.  Four source orbits impose conflicting target orbits, so
there is no total quotient map satisfying all constraints.  Exhausting all
6^6 maps confirms that the maximum satisfiable label constraints is 6/12, not
12/12.  Thus the missing successor/evolution atom remains unexported.
"""
from __future__ import annotations

import hashlib, itertools, json, math
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3015_s1965_unit_action_compatible_time_observable_functor_obstruction import OUT as P3015

OUT = GEN / "p3016_s1966_quotient_clock_successor_semigroup_exhaustion.json"
MD = GEN / "p3016_s1966_quotient_clock_successor_semigroup_exhaustion.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = [1, 5, 7, 11]


def normalize_label(x: int) -> int:
    r = x % N
    return N if r == 0 else r


def act(unit: int, d: int) -> int:
    return normalize_label(unit * d)


def orbit_of(d: int) -> tuple[int, ...]:
    return tuple(sorted({act(u, d) for u in UNITS}))


def orbit_name(orbit: tuple[int, ...]) -> str:
    return "{" + ",".join(str(x) for x in orbit) + "}"


def build_successor_exhaustion() -> dict[str, Any]:
    orbits = sorted({orbit_of(d) for d in range(1, N + 1)}, key=lambda o: (len(o), o))
    orbit_index = {orbit: i for i, orbit in enumerate(orbits)}
    label_constraints = []
    constraints_by_source: dict[int, list[int]] = defaultdict(list)
    for d in range(1, N + 1):
        src = orbit_index[orbit_of(d)]
        tgt = orbit_index[orbit_of(normalize_label(d + 1))]
        label_constraints.append({
            "d": d,
            "source_orbit": orbit_name(orbits[src]),
            "successor_label": normalize_label(d + 1),
            "target_orbit": orbit_name(orbits[tgt]),
            "source_index": src,
            "target_index": tgt,
        })
        constraints_by_source[src].append(tgt)

    source_rows = []
    forced_total_map_possible = True
    local_best_total = 0
    for src in range(len(orbits)):
        targets = constraints_by_source[src]
        counts = Counter(targets)
        conflict = len(counts) > 1
        forced_total_map_possible = forced_total_map_possible and not conflict
        local_best_total += max(counts.values())
        source_rows.append({
            "source_orbit": orbit_name(orbits[src]),
            "target_multiset": {orbit_name(orbits[k]): v for k, v in sorted(counts.items())},
            "unique_target_required": not conflict,
            "best_local_satisfied_constraints": max(counts.values()),
            "total_local_constraints": len(targets),
        })

    best_maps = []
    max_satisfied = -1
    total_maps = len(orbits) ** len(orbits)
    for image_tuple in itertools.product(range(len(orbits)), repeat=len(orbits)):
        satisfied = sum(1 for row in label_constraints if image_tuple[row["source_index"]] == row["target_index"])
        if satisfied > max_satisfied:
            max_satisfied = satisfied
            best_maps = [image_tuple]
        elif satisfied == max_satisfied and len(best_maps) < 12:
            best_maps.append(image_tuple)

    best_map_rows = [
        {orbit_name(orbits[i]): orbit_name(orbits[tgt]) for i, tgt in enumerate(image_tuple)}
        for image_tuple in best_maps
    ]

    obligations = [
        {"obligation": "unit_quotient_domain", "satisfied": True, "detail": "domain is the six U(12)-orbits from P3015"},
        {"obligation": "total_successor_map_on_quotient", "satisfied": True, "detail": "candidate space exhausts all total functions from six orbits to six orbits"},
        {"obligation": "represents_all_label_successors", "satisfied": max_satisfied == N, "detail": f"best total quotient maps satisfy {max_satisfied}/{N} label successor constraints"},
        {"obligation": "nonconflicting_source_orbit_targets", "satisfied": forced_total_map_possible, "detail": "each source orbit must impose one target orbit"},
        {"obligation": "directed_semigroup_time_evolution_export", "satisfied": False, "detail": "no strict composition/evolution law is exported after the successor constraints fail"},
        {"obligation": "eom_hamiltonian_installation", "satisfied": False, "detail": "successor exhaustion supplies no unit-bearing action/EOM/Hamiltonian source"},
    ]
    return {
        "object": "QuotientClockSuccessorSemigroup_ExhaustionNoGoMatrix",
        "orbits": [orbit_name(o) for o in orbits],
        "label_constraints": label_constraints,
        "source_conflict_rows": source_rows,
        "conflicting_source_orbit_count": sum(1 for row in source_rows if not row["unique_target_required"]),
        "total_candidate_maps_exhausted": total_maps,
        "max_satisfied_label_constraints": max_satisfied,
        "required_label_constraints": N,
        "best_maps_sample": best_map_rows,
        "local_best_upper_bound": local_best_total,
        "proof_obligations": obligations,
        "accepted_as_directed_successor_semigroup": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3015_path: Any) -> dict[str, Any]:
    matrix = build_successor_exhaustion()
    return {
        "status": "P3016_QUOTIENT_CLOCK_SUCCESSOR_SEMIGROUP_EXHAUSTION_NO_GO",
        "input_hashes": {"P3015": hashlib.sha256(p3015_path.read_bytes()).hexdigest() if p3015_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "orbit_count": len(matrix["orbits"]),
            "label_constraint_count": len(matrix["label_constraints"]),
            "conflicting_source_orbit_count": matrix["conflicting_source_orbit_count"],
            "total_candidate_maps_exhausted": matrix["total_candidate_maps_exhausted"],
            "max_satisfied_label_constraints": matrix["max_satisfied_label_constraints"],
            "required_label_constraints": matrix["required_label_constraints"],
            "accepted_as_directed_successor_semigroup": matrix["accepted_as_directed_successor_semigroup"],
        },
        "decision": {
            "breakthrough": "The missing successor/evolution atom from P3015 was attacked directly by exhausting all total maps on the six-orbit U(12) quotient.  No quotient successor can represent d -> d+1 for all labels; the best maps satisfy only 6/12 label constraints, so no strict directed time semigroup is exported.",
            "negative_export_flags": {k: False for k in ["directed_successor_semigroup_exported", "time_observable_generator_exported", "time_arrow_exported", "unit_bearing_eom_hamiltonian_exported", "ltotal_exported", "observed_physics_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay quotient successor maps or orbit-average time readouts.  The next proof-grade move should pivot to exactly one unit-bearing action/EOM source for one already typed observable, or introduce a genuinely new strict time-order object outside the U(12)-orbit quotient; keep selector, bridge, role-transfer, L_total, observed-physics, and ToE closure blocked until that object passes finite acceptance.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3016/S1966 quotient clock-successor semigroup exhaustion", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- orbit count: `{c['orbit_count']}`",
        f"- label constraints: `{c['label_constraint_count']}`",
        f"- conflicting source orbits: `{c['conflicting_source_orbit_count']}`",
        f"- total candidate maps exhausted: `{c['total_candidate_maps_exhausted']}`",
        f"- max/required satisfied constraints: `{c['max_satisfied_label_constraints']}/{c['required_label_constraints']}`",
        f"- accepted as directed successor semigroup: `{c['accepted_as_directed_successor_semigroup']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3015)
    payload = build_payload(P3015)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3016/S1966 quotient clock-successor semigroup exhaustion", "## P3016/S1966 quotient clock-successor semigroup exhaustion\n\n`P3016/S1966` attacks the P3015 missing successor/evolution atom directly.  On the six U(12)-orbits it exhausts all `6^6 = 46656` total quotient maps and asks whether any map can represent the label successor `d -> d+1` for all twelve labels.  Four source orbits impose conflicting target orbits, and the best total maps satisfy only `6/12` label-successor constraints.  Therefore no strict unit-compatible directed successor/semigroup is exported from this quotient, and no time observable generator, time arrow, EOM/Hamiltonian, `L_total`, observed-physics export, bridge/role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3016/S1966 quotient successor exhaustion `L_total` guard", "## P3016/S1966 quotient successor exhaustion `L_total` guard\n\n`P3016/S1966` adds no `L_total` term.  Exhausting all total maps on the U(12)-orbit quotient proves that the label successor `d -> d+1` cannot descend to a strict directed successor/semigroup on that quotient; the best maps satisfy only `6/12` constraints and provide no unit-bearing action, variational density, EOM, Hamiltonian, selector source, bridge closure, or role-transfer theorem.\n")
    append_once(AGENTS, "Current quotient clock-successor semigroup guardrail (P3016/S1966, 2026-06-22)", "## Current quotient clock-successor semigroup guardrail (P3016/S1966, 2026-06-22)\n\n- P3016 exhausts all `6^6` total maps on the P3015 U(12)-orbit quotient to test whether a strict unit-compatible directed successor/semigroup can represent `d -> d+1`.\n- The finite result is bounded no-go: four source orbits impose conflicting target orbits and the best total maps satisfy only `6/12` label-successor constraints.\n- Do not replay quotient successor maps or orbit-average time readouts as time-arrow, EOM/Hamiltonian, `L_total`, observed-physics, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should pivot to exactly one unit-bearing action/EOM source for one already typed observable, or introduce a genuinely new strict time-order object outside the U(12)-orbit quotient.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
