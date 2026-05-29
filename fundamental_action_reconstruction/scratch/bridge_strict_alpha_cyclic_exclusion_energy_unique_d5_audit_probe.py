#!/usr/bin/env python3
"""Scratch probe: cyclic exclusion energy uniquely selects A5/d5, conditionally.

Previous audits found:

1. cyclic nearest-neighbour adjacency can conditionally name the {1,11} kernel
   of chi_11;
2. nearest-neighbour minimization d1 alone does not globally select A5/d5;
3. a lexicographic min-d1/max-d5 tie-break does select A5/d5, but imports a
   fifth-shell preference.

This probe asks for a cleaner finite variational certificate.  It tests the
orientation-free cyclic exclusion energy

    E_excl(S) = d1(S) + d6(S),

where d1 counts nearest-neighbour contacts and d6 counts antipodal contacts in
the folded Z_12 distance histogram of a 5-support S.

Exact result: among all C(12,5)=792 supports, E_excl has minimum 0, with 12
support winners forming exactly one dihedral orbit.  That orbit has histogram
[0,3,2,1,4,0], i.e. the A5/d5 orbit.  This is a stronger finite selector
certificate than d1 alone.

No false pass: the probe does not derive the exclusion energy from strict
nadsoliton geometry.  It only shows that if a strict source exports the
orientation-free rule "avoid adjacent and antipodal contacts", then the A5/d5
support orbit is selected without separately importing a directed unit-axis bit.
QW-2191 remains open until that source premise is proven or explicitly added.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_cyclic_exclusion_energy_unique_d5_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_cyclic_exclusion_energy_unique_d5_audit_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1_SUPPORT = tuple(range(5))
A5_SUPPORT = (0, 5, 10, 3, 8)
A1_NAME = "A1_k1_contiguous"
A5_NAME = "A5_k5_d5"
SHELL_NAMES = ["d1_nearest", "d2", "d3", "d4", "d5", "d6_antipodal"]


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def canonical_support(support: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(value % N for value in support))


def dihedral_orbit(support: tuple[int, ...]) -> set[tuple[int, ...]]:
    support_set = set(support)
    orbit = set()
    for shift in range(N):
        orbit.add(tuple(sorted((value + shift) % N for value in support_set)))
        orbit.add(tuple(sorted((-value + shift) % N for value in support_set)))
    return orbit


def orbit_representatives(supports: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    remaining = set(supports)
    rows = []
    while remaining:
        seed = min(remaining)
        orbit = dihedral_orbit(seed)
        members = sorted(remaining & orbit)
        remaining -= orbit
        representative = members[0]
        hist = distance_histogram(representative)
        rows.append(
            {
                "representative": list(representative),
                "orbit_size": len(members),
                "distance_histogram_d1_to_d6": list(hist),
                "exclusion_energy_d1_plus_d6": hist[0] + hist[5],
                "is_A5_d5_orbit": canonical_support(A5_SUPPORT) in orbit,
                "is_A1_contiguous_orbit": canonical_support(A1_SUPPORT) in orbit,
            }
        )
    return sorted(rows, key=lambda row: (row["exclusion_energy_d1_plus_d6"], row["representative"]))


def exclusion_energy(support: tuple[int, ...]) -> int:
    hist = distance_histogram(support)
    return hist[0] + hist[5]


def energy_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    energy_counter = Counter(exclusion_energy(support) for support in supports)
    minimum = min(energy_counter)
    winners = [support for support in supports if exclusion_energy(support) == minimum]
    winner_histograms = Counter(distance_histogram(support) for support in winners)
    winner_orbits = orbit_representatives(winners)
    return {
        "energy_definition": "E_excl = d1_nearest + d6_antipodal",
        "energy_distribution": [
            {"energy": energy, "support_count": count}
            for energy, count in sorted(energy_counter.items())
        ],
        "minimum_energy": minimum,
        "winner_support_count": len(winners),
        "winner_histogram_class_count": len(winner_histograms),
        "winner_dihedral_orbit_count": len(winner_orbits),
        "winner_histograms": [
            {"distance_histogram_d1_to_d6": list(hist), "support_count": count}
            for hist, count in sorted(winner_histograms.items())
        ],
        "winner_orbits": winner_orbits,
        "selects_A5_d5_orbit": len(winner_orbits) == 1 and winner_orbits[0]["is_A5_d5_orbit"],
        "selects_A1_contiguous_orbit": any(row["is_A1_contiguous_orbit"] for row in winner_orbits),
    }


def d1_zero_refinement_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    d1_zero = [support for support in supports if distance_histogram(support)[0] == 0]
    d1_zero_hist_counter = Counter(distance_histogram(support) for support in d1_zero)
    return {
        "d1_zero_support_count": len(d1_zero),
        "d1_zero_histogram_class_count": len(d1_zero_hist_counter),
        "d1_zero_dihedral_orbit_count": len(orbit_representatives(d1_zero)),
        "d1_zero_histogram_rows": [
            {
                "distance_histogram_d1_to_d6": list(hist),
                "d6_antipodal_count": hist[5],
                "support_count": count,
                "survives_d1_plus_d6_zero": hist[5] == 0,
                "is_A5_d5_histogram": list(hist) == list(distance_histogram(A5_SUPPORT)),
            }
            for hist, count in sorted(d1_zero_hist_counter.items())
        ],
    }


def binary_shell_mask_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    selecting_masks = []
    for mask in range(1, 1 << len(SHELL_NAMES)):
        shell_indices = [index for index in range(len(SHELL_NAMES)) if mask & (1 << index)]
        scored = [
            (sum(distance_histogram(support)[index] for index in shell_indices), support)
            for support in supports
        ]
        minimum = min(score for score, _support in scored)
        winners = [support for score, support in scored if score == minimum]
        winner_orbits = orbit_representatives(winners)
        selects_a5 = len(winner_orbits) == 1 and winner_orbits[0]["is_A5_d5_orbit"]
        if selects_a5:
            selecting_masks.append(
                {
                    "penalized_shell_indices_1_based": [index + 1 for index in shell_indices],
                    "penalized_shell_names": [SHELL_NAMES[index] for index in shell_indices],
                    "minimum_energy": minimum,
                    "winner_support_count": len(winners),
                    "winner_dihedral_orbit_count": len(winner_orbits),
                }
            )

    selecting_sets = [set(row["penalized_shell_indices_1_based"]) for row in selecting_masks]
    minimal_rows = []
    for row, shell_set in zip(selecting_masks, selecting_sets):
        if not any(other < shell_set for other in selecting_sets):
            minimal_rows.append(row)

    return {
        "binary_masks_tested": (1 << len(SHELL_NAMES)) - 1,
        "A5_selecting_binary_mask_count": len(selecting_masks),
        "A5_selecting_binary_masks": selecting_masks,
        "inclusion_minimal_A5_selecting_masks": minimal_rows,
        "contains_d1_plus_d6_mask": any(
            row["penalized_shell_indices_1_based"] == [1, 6]
            for row in selecting_masks
        ),
        "honest_status": "classification of simple binary shell penalties only; it is not a derivation of the penalty source",
    }


def pairwise_A1_A5_audit() -> dict[str, Any]:
    a1_hist = distance_histogram(A1_SUPPORT)
    a5_hist = distance_histogram(A5_SUPPORT)
    return {
        "A1_name": A1_NAME,
        "A1_support": list(A1_SUPPORT),
        "A1_histogram_d1_to_d6": list(a1_hist),
        "A1_exclusion_energy_d1_plus_d6": a1_hist[0] + a1_hist[5],
        "A5_name": A5_NAME,
        "A5_support": list(canonical_support(A5_SUPPORT)),
        "A5_histogram_d1_to_d6": list(a5_hist),
        "A5_exclusion_energy_d1_plus_d6": a5_hist[0] + a5_hist[5],
        "exclusion_energy_prefers_A5_over_A1_pairwise": (a5_hist[0] + a5_hist[5]) < (a1_hist[0] + a1_hist[5]),
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "enumeration_domain": "All 5-subsets of Z_12 are enumerated exactly: C(12,5)=792.",
        "energy_definition": "E_excl(S)=d1(S)+d6(S), an orientation-free penalty for nearest-neighbour and antipodal contacts.",
        "unique_orbit_fact": "The exact minimum is E_excl=0; the 12 winning supports form one dihedral orbit with histogram [0,3,2,1,4,0].",
        "d1_refinement": "The earlier d1=0 condition leaves 36 supports in 3 dihedral orbits; adding d6=0 removes the two non-A5 d1=0 histogram classes.",
        "binary_mask_classification": "Among 63 nonempty binary shell penalties, 7 select A5/d5; the inclusion-minimal selecting masks are {d1,d4} and {d1,d6}.",
        "missing_source": "The exclusion energy is an exact conditional selector witness, not a derivation of the energy from strict nadsoliton geometry.",
    }


def main() -> None:
    supports = all_supports()
    all_histograms = {distance_histogram(support) for support in supports}
    pairwise = pairwise_A1_A5_audit()
    energy = energy_audit(supports)
    d1_refinement = d1_zero_refinement_audit(supports)
    mask_audit = binary_shell_mask_audit(supports)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_CYCLIC_EXCLUSION_ENERGY_UNIQUE_D5_AUDIT_PROBE__CONDITIONAL_NOT_A_THEOREM",
        "status": "d1_plus_d6_exclusion_energy-uniquely-selects-A5-d5-orbit-conditionally-source-not-derived",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(all_histograms),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "pairwise_A1_A5_audit": pairwise,
        "cyclic_exclusion_energy_audit": energy,
        "d1_zero_refinement_audit": d1_refinement,
        "binary_shell_mask_audit": mask_audit,
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "The orientation-free exclusion energy d1+d6 uniquely selects the A5/d5 dihedral orbit among all 5-subsets of Z_12.",
            "improvement_over_d1_only": "Nearest-neighbour exclusion d1=0 leaves three orbits; adding antipodal exclusion d6=0 leaves exactly the A5/d5 orbit.",
            "relation_to_previous_tie_break": "This replaces the earlier min-d1/max-d5 tie-break with a cleaner shell-exclusion witness, but it still requires a source for why d1 and d6 are the penalized shells.",
            "honest_limit": "The probe does not prove that strict nadsoliton geometry exports E_excl=d1+d6; until then the selector is conditional/non-strict.",
        },
        "ontology_guardrail": {
            "allowed_reading": "Cyclic exclusion may be investigated as a candidate structure of the nadsoliton itself.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to provide the exclusion energy or unit-axis bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives cyclic exclusion energy, d1+d6 polarity, chi_11, or the unit-axis bit from strict nadsoliton geometry.",
            "The d1+d6 exclusion-energy selector is a conditional finite witness, not a strict-core theorem.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, cyclic-adjacency source theorem, variational tie-break source theorem, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Audit existing strict-side nad12/sigma artifacts for a real source of the orientation-free d1+d6 exclusion energy; if absent, record E_excl as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha cyclic-exclusion energy unique-d5 audit probe\n\n"
        "Status: d1+d6 exclusion energy uniquely selects A5/d5 conditionally; source premise not derived here.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(all_histograms)}`.\n"
        f"- A1 exclusion energy d1+d6: `{pairwise['A1_exclusion_energy_d1_plus_d6']}`.\n"
        f"- A5 exclusion energy d1+d6: `{pairwise['A5_exclusion_energy_d1_plus_d6']}`.\n"
        f"- Minimum exclusion energy: `{energy['minimum_energy']}`.\n"
        f"- Winner support count: `{energy['winner_support_count']}`; orbit count: `{energy['winner_dihedral_orbit_count']}`.\n"
        f"- Exclusion energy selects A5/d5: `{energy['selects_A5_d5_orbit']}`.\n"
        f"- d1=0 support count before d6 refinement: `{d1_refinement['d1_zero_support_count']}`.\n"
        f"- Binary shell masks tested: `{mask_audit['binary_masks_tested']}`; A5-selecting masks: `{mask_audit['A5_selecting_binary_mask_count']}`.\n"
        f"- Inclusion-minimal A5-selecting masks: `{mask_audit['inclusion_minimal_A5_selecting_masks']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no strict exclusion-energy source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
