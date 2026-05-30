#!/usr/bin/env python3
"""Scratch probe: full-Aut block amplitude versus chi_11 polarity obstruction.

The previous conditional selector used the signed D12 shell imbalance
I=h_d5-h_d1 and found that the A1/A11 <-> A5/A7 generator is the unique maximum
among the 13 D12 chi_11 two-cycles.  This probe separates two facts that must not
be conflated:

1. the unoriented amplitude |h_d5-h_d1| is stable on full-Aut blocks and can
   locate the branch full-Aut block uniquely;
2. the sign of h_d5-h_d1 is exactly the missing chi_11 polarity and is erased
   when the two D12 components are glued back into one full-Aut orbit.

Thus full-Aut support data can at most select the block that contains the branch
pair.  It still cannot choose A5 over A1, or export the chi_11 sign, without an
extra unit-axis / shell-orientation premise.

No false pass: this is a block-location certificate, not a strict-core selector
closure.  It does not derive shell labels, does not discharge QW-2191, and does
not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_full_aut_block_amplitude_chi11_polarity_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_full_aut_block_amplitude_chi11_polarity_obstruction_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
DIHEDRAL_UNITS = [1, 11]
BRANCH_MODES = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def affine_image(support: tuple[int, ...], shift: int, unit: int) -> tuple[int, ...]:
    return tuple(sorted((shift + unit * value) % N for value in support))


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def affine_orbit(support: tuple[int, ...], units: list[int]) -> frozenset[tuple[int, ...]]:
    return frozenset(affine_image(support, shift, unit) for shift in range(N) for unit in units)


def orbit_partition(supports: list[tuple[int, ...]], units: list[int]) -> list[frozenset[tuple[int, ...]]]:
    remaining = set(supports)
    orbits: list[frozenset[tuple[int, ...]]] = []
    while remaining:
        seed = min(remaining)
        orbit = affine_orbit(seed, units)
        orbits.append(orbit)
        remaining -= orbit
    return sorted(orbits, key=lambda orbit: min(orbit))


def index_by_support(orbits: list[frozenset[tuple[int, ...]]]) -> dict[tuple[int, ...], int]:
    return {support: index for index, orbit in enumerate(orbits) for support in orbit}


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def shell_imbalance(support: tuple[int, ...]) -> int:
    histogram = distance_histogram(support)
    return histogram[4] - histogram[0]


def cyclic_gaps(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def gap_necklace(support: tuple[int, ...]) -> tuple[int, ...]:
    gaps = cyclic_gaps(support)
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_pair(name: str) -> str:
    return "d5_pair_A5_A7" if name in {"A5_k5", "A7_k7"} else "contiguous_pair_A1_A11"


def full_orbit_rows(
    full_orbits: list[frozenset[tuple[int, ...]]],
    d_orbits: list[frozenset[tuple[int, ...]]],
    d_index: dict[tuple[int, ...], int],
) -> list[dict[str, Any]]:
    branch_support_set = {unit_support(mode) for mode in BRANCH_MODES}
    rows = []
    for full_index, full_orbit in enumerate(full_orbits):
        d_components = sorted({d_index[support] for support in full_orbit})
        component_rows = []
        for component in d_components:
            representative = min(d_orbits[component])
            imbalance = shell_imbalance(representative)
            component_rows.append(
                {
                    "dihedral_component_index": component,
                    "representative_support": list(representative),
                    "gap_necklace": list(gap_necklace(representative)),
                    "histogram_d1_to_d6": list(distance_histogram(representative)),
                    "d5_minus_d1": imbalance,
                }
            )
        signed_values = [row["d5_minus_d1"] for row in component_rows]
        rows.append(
            {
                "full_affine_orbit_index": full_index,
                "full_orbit_size": len(full_orbit),
                "dihedral_component_indices": d_components,
                "dihedral_component_count": len(d_components),
                "component_rows": component_rows,
                "signed_d5_minus_d1_values_by_component": signed_values,
                "unoriented_abs_d5_minus_d1_amplitude": max(abs(value) for value in signed_values),
                "has_polarity_pair": any(value > 0 for value in signed_values) and any(value < 0 for value in signed_values),
                "contains_branch_supports": bool(full_orbit & branch_support_set),
            }
        )
    return rows


def branch_rows(d_index: dict[tuple[int, ...], int], full_index: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        rows.append(
            {
                "name": name,
                "support": list(support),
                "branch_pair": branch_pair(name),
                "required_chi11_value": 1 if branch_pair(name) == "d5_pair_A5_A7" else -1,
                "dihedral_orbit_index": d_index[support],
                "full_affine_orbit_index": full_index[support],
                "d5_minus_d1": shell_imbalance(support),
            }
        )
    return rows


def exact_proof_certificate(max_amp: int, max_count: int, distribution: dict[int, int]) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are enumerated and partitioned into D12 and full affine Aut orbits.",
        "amplitude_invariance": "The quantity |h_d5-h_d1| is unchanged when unit 5 swaps d1 and d5, so it is a full-Aut block amplitude rather than a polarity.",
        "unique_block_location": f"Across the 25 full-Aut blocks, max amplitude is {max_amp}, attained by {max_count} block; amplitude distribution is {distribution}.",
        "polarity_obstruction": "Inside the unique branch block the two D12 components have signed values -4 and +4, so full-Aut gluing keeps only the block and erases the chi_11 sign.",
        "selector_boundary": "The block can be located without choosing polarity, but A5 over A1 still requires a unit-axis or shell-orientation premise.",
        "strict_limit": "This is not a strict-core chi_11 source theorem and does not discharge QW-2191.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    d_orbits = orbit_partition(supports, DIHEDRAL_UNITS)
    full_orbits = orbit_partition(supports, UNITS)
    d_index = index_by_support(d_orbits)
    full_index = index_by_support(full_orbits)
    rows = full_orbit_rows(full_orbits, d_orbits, d_index)
    distribution = dict(sorted(Counter(row["unoriented_abs_d5_minus_d1_amplitude"] for row in rows).items()))
    max_amp = max(distribution)
    max_rows = [row for row in rows if row["unoriented_abs_d5_minus_d1_amplitude"] == max_amp]
    branch_block_rows = [row for row in rows if row["contains_branch_supports"]]
    assert len(max_rows) == 1
    assert len(branch_block_rows) == 1
    assert max_rows[0]["full_affine_orbit_index"] == branch_block_rows[0]["full_affine_orbit_index"]
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_FULL_AUT_BLOCK_AMPLITUDE_CHI11_POLARITY_OBSTRUCTION__BLOCK_ONLY_NOT_A_SELECTOR_THEOREM",
        "status": "full-aut-amplitude-uniquely-locates-branch-block-but-not-chi11-polarity",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "dihedral_units": DIHEDRAL_UNITS,
            "automorphism_units": UNITS,
            "d12_orbit_count": len(d_orbits),
            "full_affine_aut_orbit_count": len(full_orbits),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "block_amplitude_summary": {
            "candidate_score": "max_component_abs(h_d5-h_d1) per full-Aut block",
            "full_aut_block_score": True,
            "exports_chi11_polarity": False,
            "max_amplitude": max_amp,
            "maximizer_count": len(max_rows),
            "unique_maximizer_is_branch_full_aut_block": max_rows[0]["contains_branch_supports"],
            "amplitude_distribution": distribution,
            "branch_block_signed_values": branch_block_rows[0]["signed_d5_minus_d1_values_by_component"],
            "branch_block_has_polarity_pair": branch_block_rows[0]["has_polarity_pair"],
        },
        "branch_block_certificate": branch_block_rows[0],
        "branch_rows": branch_rows(d_index, full_index),
        "full_aut_block_rows": rows,
        "exact_proof_certificate": exact_proof_certificate(max_amp, len(max_rows), distribution),
        "interpretation": {
            "honest_positive": "A full-Aut-safe unoriented amplitude uniquely identifies the full-Aut block containing A1/A11/A5/A7.",
            "honest_negative": "The same full-Aut block contains both signs, so the chi_11 polarity remains unavailable without a unit-axis premise.",
            "relation_to_previous_probe": "The previous max-imbalance selector chose a signed D12 generator; this probe isolates the full-Aut block information that survives after forgetting the sign.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used.",
            "No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
            "Result locates a full-Aut branch block but does not select chi_11 polarity or close the strict-core selector obstruction.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    summary = payload["block_amplitude_summary"]
    branch = payload["branch_block_certificate"]
    lines = [
        "# Full-Aut block amplitude and chi_11 polarity obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- D12 orbit count: `{model['d12_orbit_count']}`",
        f"- Full affine Aut orbit count: `{model['full_affine_aut_orbit_count']}`",
        "",
        "## Block amplitude summary",
        "",
        f"- Candidate score: `{summary['candidate_score']}`",
        f"- Full-Aut block score: `{summary['full_aut_block_score']}`",
        f"- Exports chi_11 polarity: `{summary['exports_chi11_polarity']}`",
        f"- Max amplitude: `{summary['max_amplitude']}`",
        f"- Maximizer count: `{summary['maximizer_count']}`",
        f"- Unique maximizer is branch full-Aut block: `{summary['unique_maximizer_is_branch_full_aut_block']}`",
        f"- Amplitude distribution: `{summary['amplitude_distribution']}`",
        f"- Branch block signed values: `{summary['branch_block_signed_values']}`",
        "",
        "## Branch block certificate",
        "",
        f"- Full-Aut orbit index: `{branch['full_affine_orbit_index']}`",
        f"- D12 components: `{branch['dihedral_component_indices']}`",
        f"- Signed d5-d1 values: `{branch['signed_d5_minus_d1_values_by_component']}`",
        f"- Unoriented amplitude: `{branch['unoriented_abs_d5_minus_d1_amplitude']}`",
        f"- Has polarity pair: `{branch['has_polarity_pair']}`",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
