#!/usr/bin/env python3
"""Scratch probe: full-Aut invariant forbidden-shell lattice audit.

The previous two probes showed that closing the successful d1+d6 exact-cover
selector under full Aut(Z_12) gives forbidden shells {d1,d5,d6}, and that this
closed graph has independence number alpha(G)=3.  This probe checks where that
obstruction sits in the entire lattice of full-Aut-invariant pair-exclusion
systems on folded shells d1..d6.

Full Aut(Z_12) acts on folded shells with orbits

    O15={d1,d5}, O2={d2}, O3={d3}, O4={d4}, O6={d6}.

Hence there are exactly 2^5=32 full-Aut-invariant forbidden-shell masks.  This
probe enumerates all masks, computes independence profiles, identifies all masks
that make five active vertices impossible, and records the inclusion-minimal
five-UNSAT blockers.  The previous closure {d1,d5,d6} is one minimal blocker,
but not the only one; therefore the obstruction is real while the source choice
of that blocker remains conditional on the exact clause origin.

No false pass: this is a finite lattice/classification audit, not a derivation
of the chi_11 kernel or exact-cover clauses from strict geometry, not a QW-2191
discharge, and not ToE closure.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_full_aut_forbidden_shell_lattice_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_full_aut_forbidden_shell_lattice_audit_report.md"

N = 12
UNITS = [1, 5, 7, 11]
SHELLS = [1, 2, 3, 4, 5, 6]
FULL_AUT_SHELL_ORBITS = [[1, 5], [2], [3], [4], [6]]
TARGET_ACTIVE_COUNT = 5
PREVIOUS_CLOSURE = [1, 5, 6]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def shell_image(unit: int, shell: int) -> int:
    return folded(unit * shell)


def compute_shell_orbits() -> list[list[int]]:
    remaining = set(SHELLS)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = sorted({shell_image(unit, seed) for unit in UNITS})
        orbits.append(orbit)
        remaining -= set(orbit)
    return orbits


def support_satisfies(support: tuple[int, ...], forbidden_shells: set[int]) -> bool:
    return all(folded(right - left) not in forbidden_shells for left, right in combinations(support, 2))


def independent_supports(size: int, forbidden_shells: list[int]) -> list[tuple[int, ...]]:
    forbidden = set(forbidden_shells)
    return [support for support in combinations(range(N), size) if support_satisfies(support, forbidden)]


def cyclic_gap_tuple(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def canonical_necklace(gaps: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def mask_forbidden_shells(mask: int) -> list[int]:
    return sorted({shell for index, orbit in enumerate(FULL_AUT_SHELL_ORBITS) if mask & (1 << index) for shell in orbit})


def orbit_labels(mask: int) -> list[str]:
    labels = ["O15_d1_d5", "O2_d2", "O3_d3", "O4_d4", "O6_d6"]
    return [label for index, label in enumerate(labels) if mask & (1 << index)]


def independence_profile(forbidden_shells: list[int]) -> list[int]:
    return [len(independent_supports(size, forbidden_shells)) for size in range(N + 1)]


def maximum_independent_size(profile: list[int]) -> int:
    return max(size for size, count in enumerate(profile) if count > 0)


def five_gap_necklaces(forbidden_shells: list[int]) -> list[list[int]]:
    necklaces = {
        canonical_necklace(cyclic_gap_tuple(support))
        for support in independent_supports(TARGET_ACTIVE_COUNT, forbidden_shells)
    }
    return [list(necklace) for necklace in sorted(necklaces)]


def lattice_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(FULL_AUT_SHELL_ORBITS)):
        forbidden = mask_forbidden_shells(mask)
        profile = independence_profile(forbidden)
        alpha = maximum_independent_size(profile)
        five_count = profile[TARGET_ACTIVE_COUNT]
        rows.append(
            {
                "mask": mask,
                "orbit_labels": orbit_labels(mask),
                "forbidden_shells": forbidden,
                "independence_profile_k0_to_k6": profile[:7],
                "maximum_independent_size": alpha,
                "five_support_count": five_count,
                "five_active_unsat": five_count == 0,
                "five_gap_necklaces": five_gap_necklaces(forbidden),
                "is_previous_full_Aut_closure": forbidden == PREVIOUS_CLOSURE,
            }
        )
    return rows


def minimal_unsat_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    unsat_rows = [row for row in rows if row["five_active_unsat"]]
    result = []
    for row in unsat_rows:
        shell_set = set(row["forbidden_shells"])
        if not any(set(other["forbidden_shells"]) < shell_set for other in unsat_rows):
            result.append(row)
    return sorted(result, key=lambda row: (len(row["forbidden_shells"]), row["forbidden_shells"]))


def sat_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if not row["five_active_unsat"]]


def proof_certificate(rows: list[dict[str, Any]], minimal_rows: list[dict[str, Any]]) -> dict[str, str]:
    previous = next(row for row in rows if row["is_previous_full_Aut_closure"])
    return {
        "full_aut_lattice": "Aut(Z_12) has shell orbits [[1,5],[2],[3],[4],[6]], so exactly 32 full-Aut-invariant forbidden-shell masks are audited.",
        "classification_counts": f"Of 32 masks, {sum(row['five_active_unsat'] for row in rows)} make k=5 UNSAT and {sum(not row['five_active_unsat'] for row in rows)} leave at least one 5-support.",
        "minimal_blockers": f"There are {len(minimal_rows)} inclusion-minimal full-Aut k=5 blockers: {[row['forbidden_shells'] for row in minimal_rows]}.",
        "previous_closure_position": f"The previous closure {PREVIOUS_CLOSURE} is one minimal blocker with alpha={previous['maximum_independent_size']} and k=5 count={previous['five_support_count']}.",
        "non_uniqueness_warning": "Because other minimal full-Aut blockers exist, the full-Aut UNSAT fact alone does not derive the d1+d6 clause origin or chi_11 source.",
        "conditional_selector_reading": "The successful A5/d5 exact-cover selector still requires the non-full-Aut d1-vs-d5 shell-label premise before the full-Aut closure is avoided.",
    }


def build_payload() -> dict[str, Any]:
    rows = lattice_rows()
    minimal_rows = minimal_unsat_rows(rows)
    satisfiable_rows = sat_rows(rows)
    previous = next(row for row in rows if row["is_previous_full_Aut_closure"])
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_FULL_AUT_FORBIDDEN_SHELL_LATTICE_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "full-aut-forbidden-shell-lattice-classified-previous-d1-d5-d6-closure-is-minimal-but-not-unique-unsat-blocker",
        "finite_model": {
            "ring": "Z_12",
            "automorphism_units": UNITS,
            "folded_shells": SHELLS,
            "computed_shell_orbits": compute_shell_orbits(),
            "full_Aut_shell_orbits_used": FULL_AUT_SHELL_ORBITS,
            "target_active_count": TARGET_ACTIVE_COUNT,
            "previous_full_Aut_closure": PREVIOUS_CLOSURE,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "lattice_summary": {
            "total_full_Aut_masks": len(rows),
            "five_active_unsat_mask_count": sum(row["five_active_unsat"] for row in rows),
            "five_active_sat_mask_count": sum(not row["five_active_unsat"] for row in rows),
            "minimal_five_unsat_mask_count": len(minimal_rows),
            "minimal_five_unsat_forbidden_shells": [row["forbidden_shells"] for row in minimal_rows],
            "previous_closure_is_minimal_five_unsat_blocker": any(row["forbidden_shells"] == PREVIOUS_CLOSURE for row in minimal_rows),
            "previous_closure_alpha": previous["maximum_independent_size"],
            "previous_closure_five_support_count": previous["five_support_count"],
        },
        "minimal_five_unsat_rows": minimal_rows,
        "five_active_satisfiable_rows": satisfiable_rows,
        "all_lattice_rows": rows,
        "exact_proof_certificate": proof_certificate(rows, minimal_rows),
        "interpretation": {
            "what_was_added": "The full-Aut invariant exclusion lattice was exhaustively classified, locating d1+d5+d6 among all invariant masks.",
            "honest_positive": "The previous full-Aut closure obstruction is inclusion-minimal and has alpha(G)=3.",
            "honest_negative": "It is not the unique minimal full-Aut k=5 blocker, so this classification still does not derive the missing chi_11/shell-label source.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced by this finite lattice audit.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, or exact-cover clauses from strict nadsoliton geometry.",
            "The result is a finite full-Aut forbidden-shell lattice audit, not a selector-origin theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["lattice_summary"]
    proof = payload["exact_proof_certificate"]
    lines = [
        "# Strict alpha full-Aut forbidden-shell lattice audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{payload['finite_model']['ring']}`",
        f"- Aut units: `{payload['finite_model']['automorphism_units']}`",
        f"- Shell orbits: `{payload['finite_model']['computed_shell_orbits']}`",
        f"- Target active count: `{payload['finite_model']['target_active_count']}`",
        f"- Previous full-Aut closure: `{payload['finite_model']['previous_full_Aut_closure']}`",
        "",
        "## Lattice summary",
        "",
        f"- Total full-Aut masks: `{summary['total_full_Aut_masks']}`",
        f"- k=5 UNSAT masks: `{summary['five_active_unsat_mask_count']}`",
        f"- k=5 SAT masks: `{summary['five_active_sat_mask_count']}`",
        f"- Minimal k=5 UNSAT masks: `{summary['minimal_five_unsat_mask_count']}`",
        f"- Minimal blockers: `{summary['minimal_five_unsat_forbidden_shells']}`",
        f"- Previous closure is minimal blocker: `{summary['previous_closure_is_minimal_five_unsat_blocker']}`",
        f"- Previous closure alpha: `{summary['previous_closure_alpha']}`",
        f"- Previous closure k=5 count: `{summary['previous_closure_five_support_count']}`",
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
