#!/usr/bin/env python3
"""Scratch probe: chi_11-kernel shell-exclusion lattice uniquely selects A5/d5.

Previous full-Aut audits showed the cost of closing the successful d1+d6
exact-cover clauses under all units: full Aut(Z_12) sends d1 to d5, so the
closure d1+d5+d6 is UNSAT.  This probe asks the complementary finite question:
if the missing chi_11/shell-label premise is supplied, and if the allowed grammar
is exactly "forbid a subset of folded shells d1..d6 plus cardinality 5", is the
old d1+d6 selector arbitrary or unique inside that grammar?

Under the chi_11 kernel {1,11}, every folded shell label d1..d6 is stable, so
there are 2^6=64 possible shell-exclusion masks.  This probe enumerates every
mask, all C(12,5)=792 five-supports for each mask, and the resulting dihedral
orbit/histogram summaries.  Result: exactly one mask selects the A5/d5 orbit
with 12 solutions and histogram [0,3,2,1,4,0], namely {d1,d6}.

No false pass: uniqueness is only inside the explicitly axiom-augmented
chi_11-kernel + shell-exclusion grammar.  The probe does not derive chi_11, the
shell-exclusion grammar, or the exact-cover clauses from strict geometry; it is
not a QW-2191 discharge and not ToE closure.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_chi11_shell_exclusion_lattice_unique_a5_selector_report.json"
OUT_MD = HERE / "bridge_strict_alpha_chi11_shell_exclusion_lattice_unique_a5_selector_report.md"

N = 12
ACTIVE_COUNT = 5
SHELLS = [1, 2, 3, 4, 5, 6]
CHI11_KERNEL_UNITS = [1, 11]
FULL_AUT_UNITS = [1, 5, 7, 11]
TARGET_A5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
TARGET_A1_HISTOGRAM = [4, 3, 2, 1, 0, 0]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def shell_image(unit: int, shell: int) -> int:
    return folded(unit * shell)


def shell_action_rows(units: list[int]) -> list[dict[str, Any]]:
    return [
        {
            "unit": unit,
            "shell_images": {f"d{shell}": f"d{shell_image(unit, shell)}" for shell in SHELLS},
            "fixes_all_folded_shell_labels": all(shell_image(unit, shell) == shell for shell in SHELLS),
        }
        for unit in units
    ]


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


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


def cyclic_gap_tuple(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def canonical_necklace(gaps: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def forbidden_shells_for_mask(mask: int) -> list[int]:
    return [shell for index, shell in enumerate(SHELLS) if mask & (1 << index)]


def mask_row(mask: int) -> dict[str, Any]:
    forbidden_shells = forbidden_shells_for_mask(mask)
    solutions = enumerate_solutions(forbidden_shells)
    orbits = dihedral_orbits(solutions)
    histograms = sorted({distance_histogram(support) for support in solutions})
    gap_necklaces = sorted({canonical_necklace(cyclic_gap_tuple(support)) for support in solutions})
    selects_a5_d5 = len(solutions) == 12 and len(orbits) == 1 and histograms == [tuple(TARGET_A5_HISTOGRAM)]
    selects_a1_contiguous = len(solutions) == 12 and len(orbits) == 1 and histograms == [tuple(TARGET_A1_HISTOGRAM)]
    return {
        "mask": mask,
        "forbidden_shells": forbidden_shells,
        "forbidden_shell_labels": [f"d{shell}" for shell in forbidden_shells],
        "solution_count": len(solutions),
        "dihedral_orbit_count": len(orbits),
        "histograms_d1_to_d6": [list(histogram) for histogram in histograms],
        "gap_necklaces": [list(necklace) for necklace in gap_necklaces],
        "selects_A5_d5": selects_a5_d5,
        "selects_A1_contiguous": selects_a1_contiguous,
        "is_unsat": len(solutions) == 0,
        "representative_solutions": [list(support) for support in solutions[:6]],
    }


def lattice_rows() -> list[dict[str, Any]]:
    return [mask_row(mask) for mask in range(1 << len(SHELLS))]


def minimal_rows(rows: list[dict[str, Any]], predicate_key: str) -> list[dict[str, Any]]:
    selected = [row for row in rows if row[predicate_key]]
    result = []
    for row in selected:
        shell_set = set(row["forbidden_shells"])
        if not any(set(other["forbidden_shells"]) < shell_set for other in selected):
            result.append(row)
    return sorted(result, key=lambda row: (len(row["forbidden_shells"]), row["forbidden_shells"]))


def exact_proof_certificate(rows: list[dict[str, Any]]) -> dict[str, str]:
    a5_rows = [row for row in rows if row["selects_A5_d5"]]
    a1_rows = [row for row in rows if row["selects_A1_contiguous"]]
    unsat_count = sum(row["is_unsat"] for row in rows)
    return {
        "finite_domain": "For every one of the 64 chi_11-stable shell-exclusion masks, all C(12,5)=792 supports are enumerated exactly.",
        "chi11_shell_stability": "Units {1,11} fix every folded shell label d1..d6; unlike full Aut, they do not force d1 and d5 to share a clause orbit.",
        "unique_a5_selector": f"Exactly {len(a5_rows)} mask selects the A5/d5 orbit; its forbidden shells are {a5_rows[0]['forbidden_shells'] if a5_rows else None}.",
        "conjugate_boundary": f"Exactly {len(a1_rows)} mask selects the A1/contiguous orbit; this records the d5+d6 conjugate boundary separately from the A5 selector.",
        "unsat_boundary": f"There are {unsat_count} masks with no 5-support at all, so uniqueness of A5/d5 is not just generic overconstraint.",
        "conditional_scope": "This is uniqueness only after adjoining the chi_11/shell-label premise and restricting the grammar to pure shell-exclusion exact-cover clauses.",
    }


def build_payload() -> dict[str, Any]:
    rows = lattice_rows()
    a5_rows = [row for row in rows if row["selects_A5_d5"]]
    a1_rows = [row for row in rows if row["selects_A1_contiguous"]]
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_CHI11_SHELL_EXCLUSION_LATTICE_UNIQUE_A5_SELECTOR_PROBE__CONDITIONAL_NOT_A_THEOREM",
        "status": "chi11-shell-exclusion-lattice-has-unique-a5-d5-selector-d1-d6-conditional-on-shell-label-premise",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "folded_shells": SHELLS,
            "chi11_kernel_units": CHI11_KERNEL_UNITS,
            "full_Aut_units_for_contrast": FULL_AUT_UNITS,
            "target_A5_histogram": TARGET_A5_HISTOGRAM,
            "target_A1_histogram": TARGET_A1_HISTOGRAM,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "shell_action_audit": {
            "chi11_kernel_shell_action": shell_action_rows(CHI11_KERNEL_UNITS),
            "full_Aut_shell_action_for_contrast": shell_action_rows(FULL_AUT_UNITS),
        },
        "lattice_summary": {
            "total_chi11_shell_masks": len(rows),
            "a5_d5_selecting_mask_count": len(a5_rows),
            "a5_d5_selecting_forbidden_shells": [row["forbidden_shells"] for row in a5_rows],
            "a5_d5_selecting_masks": [row["mask"] for row in a5_rows],
            "a1_contiguous_selecting_mask_count": len(a1_rows),
            "a1_contiguous_selecting_forbidden_shells": [row["forbidden_shells"] for row in a1_rows],
            "unsat_mask_count": sum(row["is_unsat"] for row in rows),
            "satisfiable_mask_count": sum(not row["is_unsat"] for row in rows),
            "minimal_unsat_forbidden_shells": [row["forbidden_shells"] for row in minimal_rows(rows, "is_unsat")],
        },
        "a5_d5_selector_rows": a5_rows,
        "a1_contiguous_selector_rows": a1_rows,
        "all_lattice_rows": rows,
        "exact_proof_certificate": exact_proof_certificate(rows),
        "interpretation": {
            "honest_positive": "Within the chi_11-kernel pure shell-exclusion grammar, d1+d6 is the unique finite selector of the A5/d5 orbit.",
            "honest_negative": "The probe does not derive the chi_11 premise or the shell-exclusion grammar from strict nadsoliton geometry.",
            "relation_to_full_Aut_audits": "Full Aut destroys the selector by adding d5; chi_11 reduction keeps d1 distinct from d5 and makes the unique A5 selector visible.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced by this finite lattice audit.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, or exact-cover clauses from strict nadsoliton geometry.",
            "The result is a conditional finite uniqueness certificate inside an axiom-augmented grammar, not a strict selector-origin theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["lattice_summary"]
    proof = payload["exact_proof_certificate"]
    lines = [
        "# Strict alpha chi_11 shell-exclusion lattice unique A5 selector audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{payload['finite_model']['ring']}`",
        f"- Active count: `{payload['finite_model']['active_count']}`",
        f"- chi_11 kernel units: `{payload['finite_model']['chi11_kernel_units']}`",
        f"- Folded shells: `{payload['finite_model']['folded_shells']}`",
        f"- Target A5/d5 histogram: `{payload['finite_model']['target_A5_histogram']}`",
        "",
        "## Lattice summary",
        "",
        f"- Total chi_11 shell masks: `{summary['total_chi11_shell_masks']}`",
        f"- A5/d5 selecting mask count: `{summary['a5_d5_selecting_mask_count']}`",
        f"- A5/d5 selecting forbidden shells: `{summary['a5_d5_selecting_forbidden_shells']}`",
        f"- A1/contiguous selecting mask count: `{summary['a1_contiguous_selecting_mask_count']}`",
        f"- A1/contiguous selecting forbidden shells: `{summary['a1_contiguous_selecting_forbidden_shells']}`",
        f"- UNSAT mask count: `{summary['unsat_mask_count']}`",
        f"- SAT mask count: `{summary['satisfiable_mask_count']}`",
        f"- Minimal UNSAT forbidden shells: `{summary['minimal_unsat_forbidden_shells']}`",
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
