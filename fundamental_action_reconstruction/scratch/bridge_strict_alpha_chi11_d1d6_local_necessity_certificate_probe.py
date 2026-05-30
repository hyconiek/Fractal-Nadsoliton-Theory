#!/usr/bin/env python3
"""Scratch probe: local necessity certificate for the chi_11 d1+d6 selector.

The previous chi_11 shell-exclusion lattice audit showed that, inside the
conditional grammar "cardinality 5 plus forbidden folded shells d1..d6", exactly
one mask selects A5/d5: {d1,d6}.  This probe makes that uniqueness more
proof-like and local by auditing the Hasse-neighborhood of {d1,d6}.

Result:
- the center {d1,d6} has 12 solutions, one dihedral orbit, histogram A5/d5;
- deleting d1 leaves {d6}, which overselects 192 supports over 9 histogram
  classes and includes A1/contiguous;
- deleting d6 leaves {d1}, which overselects 36 supports over 3 gap necklaces;
- adding any one of d2,d3,d4,d5 to {d1,d6} makes the system UNSAT.

Thus, conditional on the chi_11 shell-label premise and the pure shell-exclusion
exact-cover grammar, d1 and d6 are both necessary, and the selected mask is an
isolated singleton: downward neighbors overselect and upward neighbors overkill.

No false pass: this is a local finite certificate inside an explicitly
axiom-augmented grammar.  It does not derive chi_11, shell-exclusion grammar, or
exact-cover clauses from strict geometry; it is not a QW-2191 discharge and not
ToE closure.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_chi11_d1d6_local_necessity_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_chi11_d1d6_local_necessity_certificate_report.md"

N = 12
ACTIVE_COUNT = 5
SHELLS = [1, 2, 3, 4, 5, 6]
CENTER_FORBIDDEN_SHELLS = [1, 6]
CHI11_KERNEL_UNITS = [1, 11]
TARGET_A5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
TARGET_A1_HISTOGRAM = [4, 3, 2, 1, 0, 0]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def cyclic_gap_tuple(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def canonical_necklace(gaps: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def support_satisfies(support: tuple[int, ...], forbidden_shells: set[int]) -> bool:
    return all(folded(right - left) not in forbidden_shells for left, right in combinations(support, 2))


def enumerate_solutions(forbidden_shells: list[int]) -> list[tuple[int, ...]]:
    forbidden = set(forbidden_shells)
    return [support for support in combinations(range(N), ACTIVE_COUNT) if support_satisfies(support, forbidden)]


def translate_support(support: tuple[int, ...], shift: int) -> tuple[int, ...]:
    return tuple(sorted((value + shift) % N for value in support))


def reflect_support(support: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((-value) % N for value in support))


def dihedral_orbits(solutions: list[tuple[int, ...]]) -> list[list[tuple[int, ...]]]:
    remaining = set(solutions)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = set()
        for shift in range(N):
            orbit.add(translate_support(seed, shift))
            orbit.add(translate_support(reflect_support(seed), shift))
        orbit &= set(solutions)
        orbits.append(sorted(orbit))
        remaining -= orbit
    return sorted(orbits, key=lambda orbit: orbit[0])


def classify_histogram(histogram: tuple[int, int, int, int, int, int]) -> str:
    if list(histogram) == TARGET_A5_HISTOGRAM:
        return "A5_d5_target"
    if list(histogram) == TARGET_A1_HISTOGRAM:
        return "A1_contiguous_conjugate"
    return "other"


def system_summary(name: str, forbidden_shells: list[int], relation_to_center: str) -> dict[str, Any]:
    solutions = enumerate_solutions(forbidden_shells)
    orbits = dihedral_orbits(solutions)
    histograms = sorted({distance_histogram(support) for support in solutions})
    gap_necklaces = sorted({canonical_necklace(cyclic_gap_tuple(support)) for support in solutions})
    histogram_rows = [
        {"histogram_d1_to_d6": list(histogram), "classification": classify_histogram(histogram)}
        for histogram in histograms
    ]
    return {
        "name": name,
        "relation_to_center": relation_to_center,
        "forbidden_shells": forbidden_shells,
        "forbidden_shell_labels": [f"d{shell}" for shell in forbidden_shells],
        "solution_count": len(solutions),
        "dihedral_orbit_count": len(orbits),
        "histogram_class_count": len(histograms),
        "histogram_rows": histogram_rows,
        "gap_necklace_count": len(gap_necklaces),
        "gap_necklaces": [list(necklace) for necklace in gap_necklaces],
        "selects_A5_d5_uniquely": len(solutions) == 12 and len(orbits) == 1 and histograms == [tuple(TARGET_A5_HISTOGRAM)],
        "contains_A5_d5_histogram": tuple(TARGET_A5_HISTOGRAM) in histograms,
        "contains_A1_contiguous_histogram": tuple(TARGET_A1_HISTOGRAM) in histograms,
        "is_unsat": len(solutions) == 0,
        "representative_solutions": [list(support) for support in solutions[:6]],
    }


def local_boundary_rows() -> list[dict[str, Any]]:
    center = set(CENTER_FORBIDDEN_SHELLS)
    rows = [system_summary("center_d1_plus_d6", sorted(center), "center")]
    for shell in CENTER_FORBIDDEN_SHELLS:
        neighbor = sorted(center - {shell})
        rows.append(system_summary(f"delete_d{shell}", neighbor, f"delete d{shell}"))
    for shell in SHELLS:
        if shell in center:
            continue
        neighbor = sorted(center | {shell})
        rows.append(system_summary(f"add_d{shell}", neighbor, f"add d{shell}"))
    return rows


def exact_proof_certificate(rows: list[dict[str, Any]]) -> dict[str, str]:
    by_name = {row["name"]: row for row in rows}
    upward = [row for row in rows if row["relation_to_center"].startswith("add")]
    return {
        "finite_domain": "Each local system enumerates all C(12,5)=792 supports exactly under the chi_11-stable shell-exclusion grammar.",
        "center_certificate": "The center {d1,d6} has 12 solutions, one dihedral orbit, and only histogram [0,3,2,1,4,0], so it selects A5/d5.",
        "d1_necessity": f"Deleting d1 leaves only d6 and gives {by_name['delete_d1']['solution_count']} solutions over {by_name['delete_d1']['histogram_class_count']} histogram classes, including A1/contiguous; hence d1 is necessary for specificity.",
        "d6_necessity": f"Deleting d6 leaves only d1 and gives {by_name['delete_d6']['solution_count']} solutions over {by_name['delete_d6']['gap_necklace_count']} gap necklaces; hence d6 is necessary for exact cover uniqueness.",
        "maximality_boundary": f"Adding any one missing shell to {{d1,d6}} gives UNSAT counts {[row['solution_count'] for row in upward]}, so the selector is upward-isolated.",
        "conditional_scope": "Necessity is local and conditional on the chi_11/shell-label premise plus the pure shell-exclusion exact-cover grammar; no strict source theorem is claimed.",
    }


def build_payload() -> dict[str, Any]:
    rows = local_boundary_rows()
    center = next(row for row in rows if row["relation_to_center"] == "center")
    deletion_rows = [row for row in rows if row["relation_to_center"].startswith("delete")]
    addition_rows = [row for row in rows if row["relation_to_center"].startswith("add")]
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_CHI11_D1D6_LOCAL_NECESSITY_CERTIFICATE_PROBE__CONDITIONAL_NOT_A_THEOREM",
        "status": "d1-d6-is-locally-necessary-and-upward-isolated-inside-chi11-shell-exclusion-grammar",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "folded_shells": SHELLS,
            "chi11_kernel_units": CHI11_KERNEL_UNITS,
            "center_forbidden_shells": CENTER_FORBIDDEN_SHELLS,
            "target_A5_histogram": TARGET_A5_HISTOGRAM,
            "target_A1_histogram": TARGET_A1_HISTOGRAM,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "local_boundary_summary": {
            "center_selects_A5_d5_uniquely": center["selects_A5_d5_uniquely"],
            "deletion_neighbor_count": len(deletion_rows),
            "addition_neighbor_count": len(addition_rows),
            "all_deletions_overselect": all((not row["selects_A5_d5_uniquely"]) and row["solution_count"] > center["solution_count"] for row in deletion_rows),
            "all_additions_unsat": all(row["is_unsat"] for row in addition_rows),
            "delete_d1_solution_count": next(row for row in deletion_rows if row["name"] == "delete_d1")["solution_count"],
            "delete_d6_solution_count": next(row for row in deletion_rows if row["name"] == "delete_d6")["solution_count"],
            "addition_solution_counts": {row["name"]: row["solution_count"] for row in addition_rows},
        },
        "local_boundary_rows": rows,
        "exact_proof_certificate": exact_proof_certificate(rows),
        "interpretation": {
            "honest_positive": "Inside the chi_11 shell-exclusion grammar, {d1,d6} is not only globally unique but locally necessary: removing either clause overselects, adding any clause kills satisfiability.",
            "honest_negative": "This does not derive chi_11, cardinality 5, or the shell-exclusion grammar from strict nadsoliton geometry.",
            "relation_to_previous_probe": "The previous 64-mask uniqueness result is sharpened into a local Hasse-neighborhood certificate around the unique A5/d5 mask.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced by this finite local certificate.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.",
            "The result is a conditional finite local necessity certificate inside an axiom-augmented grammar, not a strict selector-origin theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["local_boundary_summary"]
    proof = payload["exact_proof_certificate"]
    lines = [
        "# Strict alpha chi_11 d1+d6 local necessity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{payload['finite_model']['ring']}`",
        f"- Active count: `{payload['finite_model']['active_count']}`",
        f"- chi_11 kernel units: `{payload['finite_model']['chi11_kernel_units']}`",
        f"- Center forbidden shells: `{payload['finite_model']['center_forbidden_shells']}`",
        "",
        "## Local boundary summary",
        "",
        f"- Center selects A5/d5 uniquely: `{summary['center_selects_A5_d5_uniquely']}`",
        f"- Deletion neighbors: `{summary['deletion_neighbor_count']}`",
        f"- Addition neighbors: `{summary['addition_neighbor_count']}`",
        f"- All deletions overselect: `{summary['all_deletions_overselect']}`",
        f"- All additions UNSAT: `{summary['all_additions_unsat']}`",
        f"- delete d1 solution count: `{summary['delete_d1_solution_count']}`",
        f"- delete d6 solution count: `{summary['delete_d6_solution_count']}`",
        f"- addition solution counts: `{summary['addition_solution_counts']}`",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in proof.items():
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
