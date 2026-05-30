#!/usr/bin/env python3
"""Scratch probe: full-Aut closure of the d1+d6 exact-cover clauses is UNSAT.

Earlier finite certificates showed that the cyclic exact-cover system

    x_i in {0,1}, sum_i x_i = 5,
    forbid folded distance d1 pairs,
    forbid folded distance d6 antipodal pairs

has exactly 12 solutions and selects the A5/d5 dihedral orbit.  Later audits
showed the symmetry cost: Aut(Z_12) units 5 and 7 send the folded shell d1 to
d5 while preserving d6.  Therefore the full-Aut orbit-closure of the same
clause family must forbid d1, d5, and d6 pairs.

This probe performs the exact finite SAT-style audit of that closure.  Result:
full-Aut clause closure is unsatisfiable for five active vertices on Z_12.  The
successful d1+d6 exact-cover certificate is therefore a conditional chi_11-kernel
certificate, not a full-Aut-invariant strict-core selector theorem.

No false pass: this is a finite no-go/conditionality audit, not a derivation of
the chi_11 kernel from strict geometry, not a QW-2191 discharge, and not ToE
closure.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_exact_cover_full_aut_clause_closure_unsat_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_exact_cover_full_aut_clause_closure_unsat_audit_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
BASE_FORBIDDEN_SHELLS = [1, 6]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def shell_image(unit: int, shell: int) -> int:
    return folded(unit * shell)


def shell_orbit(shell: int, units: list[int]) -> list[int]:
    return sorted({shell_image(unit, shell) for unit in units})


def closed_forbidden_shells(seed_shells: list[int], units: list[int]) -> list[int]:
    return sorted({image for shell in seed_shells for image in shell_orbit(shell, units)})


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def gap_necklace(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    rotations = [tuple(gaps[index:] + gaps[:index]) for index in range(len(gaps))]
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
            translated = translate_support(seed, shift)
            orbit.add(translated)
            orbit.add(translate_support(reflect_support(seed), shift))
        orbit &= set(solutions)
        orbits.append(sorted(orbit))
        remaining -= orbit
    return sorted(orbits, key=lambda orbit: orbit[0])


def solution_summary(name: str, forbidden_shells: list[int]) -> dict[str, Any]:
    solutions = enumerate_solutions(forbidden_shells)
    orbits = dihedral_orbits(solutions)
    histograms = sorted({distance_histogram(support) for support in solutions})
    necklaces = sorted({gap_necklace(support) for support in solutions})
    return {
        "name": name,
        "forbidden_shells": forbidden_shells,
        "clauses_forbidden_pair_count": sum(N // 2 if shell == 6 else N for shell in forbidden_shells),
        "cardinality_compatible_assignments_checked": 792,
        "solution_count": len(solutions),
        "dihedral_orbit_count": len(orbits),
        "representative_solutions": [list(support) for support in solutions[:6]],
        "distance_histograms_d1_to_d6": [list(histogram) for histogram in histograms],
        "gap_necklaces": [list(necklace) for necklace in necklaces],
        "selects_A5_d5_histogram": histograms == [(0, 3, 2, 1, 4, 0)],
        "is_unsat": len(solutions) == 0,
    }


def unit_shell_action_rows() -> list[dict[str, Any]]:
    return [
        {
            "unit": unit,
            "image_of_d1": shell_image(unit, 1),
            "image_of_d5": shell_image(unit, 5),
            "image_of_d6": shell_image(unit, 6),
            "preserves_base_clause_set_d1_d6": sorted({shell_image(unit, shell) for shell in BASE_FORBIDDEN_SHELLS}) == BASE_FORBIDDEN_SHELLS,
        }
        for unit in UNITS
    ]


def closure_audit() -> dict[str, Any]:
    full_aut_closure = closed_forbidden_shells(BASE_FORBIDDEN_SHELLS, UNITS)
    chi_11_kernel_units = [1, 11]
    chi_11_closure = closed_forbidden_shells(BASE_FORBIDDEN_SHELLS, chi_11_kernel_units)
    conjugate_seed = [5, 6]
    return {
        "base_seed_forbidden_shells": BASE_FORBIDDEN_SHELLS,
        "unit_shell_action": unit_shell_action_rows(),
        "shell_orbit_of_d1_under_full_Aut": shell_orbit(1, UNITS),
        "shell_orbit_of_d6_under_full_Aut": shell_orbit(6, UNITS),
        "full_Aut_closed_forbidden_shells": full_aut_closure,
        "chi_11_kernel_units": chi_11_kernel_units,
        "chi_11_kernel_closed_forbidden_shells": chi_11_closure,
        "base_chi_11_system": solution_summary("base_or_chi_11_closed_d1_plus_d6", chi_11_closure),
        "full_Aut_closed_system": solution_summary("full_Aut_closed_d1_plus_d5_plus_d6", full_aut_closure),
        "aut_conjugate_d5_plus_d6_system": solution_summary("aut_conjugate_d5_plus_d6", conjugate_seed),
    }


def exact_proof_certificate(audit: dict[str, Any]) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 Boolean cardinality-compatible supports are enumerated exactly.",
        "base_certificate": "Forbidding d1 and d6 gives 12 solutions, one dihedral orbit, with histogram [0,3,2,1,4,0] (A5/d5).",
        "unit_closure_step": "Under Aut(Z_12)={1,5,7,11}, the orbit of d1 is {d1,d5}; d6 is fixed. Full-Aut closure of d1+d6 is therefore d1+d5+d6.",
        "unsat_result": "Forbidding d1,d5,d6 has 0 solutions among all 792 five-supports, so full-Aut closure destroys the exact-cover certificate.",
        "conditional_positive": "The d1+d6 certificate survives exactly at the chi_11-kernel level {1,11}; it is not a full-Aut-invariant strict-core selector theorem.",
        "source_obstruction": "A strict derivation would still need an internal source for the chi_11/shell-label reduction before using the successful exact-cover selector.",
    }


def build_payload() -> dict[str, Any]:
    audit = closure_audit()
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_EXACT_COVER_FULL_AUT_CLAUSE_CLOSURE_UNSAT_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "full-aut-closure-of-d1-d6-exact-cover-is-unsat-chi_11-kernel-certificate-only",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": UNITS,
            "base_forbidden_shells": BASE_FORBIDDEN_SHELLS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "closure_audit": audit,
        "exact_proof_certificate": exact_proof_certificate(audit),
        "interpretation": {
            "what_was_proved": "The finite d1+d6 exact-cover selector is computationally correct but only after the d1-vs-d5 shell label has already been chosen.",
            "new_no_go": "If the clause set is closed under full Aut(Z_12), d5 clauses must be added and the five-active exact-cover problem becomes UNSAT.",
            "relation_to_previous_probe": "This strengthens the branch-pair no-go at the level of the original Boolean clauses rather than only at the level of selected branch pairs.",
            "honest_limit": "QW-2191 remains open: no strict-core source for chi_11, unit-axis bit, or shell-label reduction is derived here.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced by this finite combinatorial audit.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, or unit-axis bit from strict nadsoliton geometry.",
            "The result is a finite full-Aut clause-closure UNSAT audit, not a selector-origin theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    audit = payload["closure_audit"]
    base = audit["base_chi_11_system"]
    full = audit["full_Aut_closed_system"]
    conjugate = audit["aut_conjugate_d5_plus_d6_system"]
    lines = [
        "# Strict alpha exact-cover full-Aut clause-closure UNSAT audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{payload['finite_model']['ring']}`",
        f"- Active Boolean variables: `{payload['finite_model']['active_count']}` of 12",
        f"- Cardinality-compatible assignments checked: `{base['cardinality_compatible_assignments_checked']}`",
        f"- Aut units: `{payload['finite_model']['automorphism_units']}`",
        f"- Base forbidden shells: `{payload['finite_model']['base_forbidden_shells']}`",
        "",
        "## Clause closure computation",
        "",
        f"- Orbit of `d1` under full Aut: `{audit['shell_orbit_of_d1_under_full_Aut']}`",
        f"- Orbit of `d6` under full Aut: `{audit['shell_orbit_of_d6_under_full_Aut']}`",
        f"- Full-Aut closed forbidden shells: `{audit['full_Aut_closed_forbidden_shells']}`",
        f"- chi_11-kernel closed forbidden shells: `{audit['chi_11_kernel_closed_forbidden_shells']}`",
        "",
        "## Exact enumeration results",
        "",
        f"- d1+d6 / chi_11-kernel solution count: `{base['solution_count']}`",
        f"- d1+d6 / chi_11-kernel dihedral orbit count: `{base['dihedral_orbit_count']}`",
        f"- d1+d6 selects A5/d5 histogram: `{base['selects_A5_d5_histogram']}`",
        f"- Full-Aut closure d1+d5+d6 solution count: `{full['solution_count']}`",
        f"- Full-Aut closure is UNSAT: `{full['is_unsat']}`",
        f"- Aut-conjugate d5+d6 solution count: `{conjugate['solution_count']}`",
        f"- Aut-conjugate d5+d6 histograms: `{conjugate['distance_histograms_d1_to_d6']}`",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend([
        "",
        "## Hard limits",
        "",
    ])
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
