#!/usr/bin/env python3
"""P2739/S1689: sign-torsor quotient no-section certificate.

This strengthens P2738 with a linear-algebra section obstruction.  P2738
exhausted Boolean laws; P2739 proves why the result is forced: the simultaneous
flip quotient has no nonzero section of the anti-equivariant sign line that is
also a well-defined absolute quotient polarity.  The proof is finite and
executable on the 16 sign states.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2739_s1689_sign_torsor_quotient_no_section_certificate.json"
MD = GEN / "p2739_s1689_sign_torsor_quotient_no_section_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P2738_BOOLEAN_EXHAUSTION": GEN / "p2738_s1688_sign_torsor_boolean_source_law_exhaustion.json",
}

CONTENT_PATTERNS = {
    "boolean_exhaustion_boundary": r"Boolean laws|Boolean source-law|2\^16.*Boolean|sign-torsor Boolean",
    "simultaneous_flip_quotient": r"simultaneous-flip quotient|simultaneous flip|global flip quotient|global sign",
    "absolute_polarity_block": r"absolute.*polarity|non-premise.*polarity|lambda/P2721 polarity|P2721 polarity",
    "equivariant_pairing_block": r"equivariant.*paired|flip-paired|sign-paired|lambda -> -lambda|opposite polarity",
}
NEGATIVE_EXPORT_FLAGS = [
    "nonzero_quotient_sign_section_exported",
    "strict_signed_source_law_exported",
    "lambda_fixed",
    "p2721_polarity_selected",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {
        "content_pattern_count": len(CONTENT_PATTERNS),
        "rows": rows,
        "hit_counts": {row["content_lane"]: row["hit_count"] for row in rows},
        "all_patterns_have_hits": all(row["hit_count"] > 0 for row in rows),
    }


def sign_states() -> list[tuple[int, int, int, int]]:
    return list(itertools.product((-1, 1), repeat=4))


def flip(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(-value for value in state)


def rational_rank(matrix: list[list[int]]) -> int:
    rows = [[Fraction(value) for value in row] for row in matrix if any(value != 0 for value in row)]
    rank = 0
    col_count = len(rows[0]) if rows else 0
    for col in range(col_count):
        pivot = None
        for r in range(rank, len(rows)):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][col]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [rows[r][c] - factor * rows[rank][c] for c in range(col_count)]
        rank += 1
        if rank == len(rows):
            break
    return rank


def build_orbits(states: list[tuple[int, int, int, int]]) -> list[tuple[int, int]]:
    index = {state: i for i, state in enumerate(states)}
    seen: set[tuple[int, int]] = set()
    orbits = []
    for state in states:
        pair = tuple(sorted((index[state], index[flip(state)])))
        if pair not in seen:
            seen.add(pair)
            orbits.append(pair)
    return orbits


def row_for_pair(n: int, i: int, j: int, relation: str) -> list[int]:
    row = [0] * n
    row[i] = 1
    if relation == "equal":
        row[j] = -1  # s_i - s_j = 0
    elif relation == "opposite":
        row[j] = 1  # s_i + s_j = 0
    else:
        raise ValueError(relation)
    return row


def section_obstruction() -> dict[str, Any]:
    states = sign_states()
    n = len(states)
    orbits = build_orbits(states)
    invariant_rows = [row_for_pair(n, i, j, "equal") for i, j in orbits]
    anti_rows = [row_for_pair(n, i, j, "opposite") for i, j in orbits]
    combined_rows = invariant_rows + anti_rows
    invariant_rank = rational_rank(invariant_rows)
    anti_rank = rational_rank(anti_rows)
    combined_rank = rational_rank(combined_rows)
    invariant_nullity = n - invariant_rank
    anti_nullity = n - anti_rank
    combined_nullity = n - combined_rank

    # Boolean cross-check: no {±1}-valued section can satisfy both equations.
    boolean_both_count = 0
    boolean_invariant_count = 0
    boolean_anti_count = 0
    for values in itertools.product((-1, 1), repeat=n):
        invariant_ok = all(values[i] == values[j] for i, j in orbits)
        anti_ok = all(values[i] == -values[j] for i, j in orbits)
        if invariant_ok:
            boolean_invariant_count += 1
        if anti_ok:
            boolean_anti_count += 1
        if invariant_ok and anti_ok:
            boolean_both_count += 1

    return {
        "variables": ["lambda_sign", "orientation_sign", "P2721_polarity_sign", "branch_square_flux_sign"],
        "state_count": n,
        "orbit_count": len(orbits),
        "orbits": [{"i": i, "j": j, "state_i": states[i], "state_j": states[j]} for i, j in orbits],
        "linear_system": {
            "invariant_quotient_descent_rank": invariant_rank,
            "invariant_quotient_descent_nullity": invariant_nullity,
            "anti_equivariant_sign_line_rank": anti_rank,
            "anti_equivariant_sign_line_nullity": anti_nullity,
            "combined_rank": combined_rank,
            "combined_nullity": combined_nullity,
        },
        "boolean_cross_check": {
            "total_pm1_sections": 2 ** n,
            "invariant_pm1_sections": boolean_invariant_count,
            "anti_equivariant_pm1_sections": boolean_anti_count,
            "simultaneously_invariant_and_anti_pm1_sections": boolean_both_count,
        },
        "theorem": "A non-premise absolute polarity would have to descend to the simultaneous-flip quotient, so it must satisfy s(x)=s(Fx).  A sourced sign-line response to the missing torsor would have to be anti-equivariant, so it must satisfy s(x)=-s(Fx).  The combined rational linear system has rank 16 and nullity 0, and the {±1}-valued cross-check has zero simultaneous sections.  Thus the quotient has no nonzero sign-line section capable of fixing lambda/P2721.",
    }


def acceptance_matrix(scan: dict[str, Any], obstruction: dict[str, Any]) -> dict[str, Any]:
    linear = obstruction["linear_system"]
    boolean = obstruction["boolean_cross_check"]
    facts = {
        "content_grep_confirms_no_section_lane_is_not_new_source": scan["all_patterns_have_hits"],
        "quotient_orbits_free": obstruction["state_count"] == 2 * obstruction["orbit_count"],
        "invariant_sections_exist": linear["invariant_quotient_descent_nullity"] == obstruction["orbit_count"],
        "anti_equivariant_sections_exist": linear["anti_equivariant_sign_line_nullity"] == obstruction["orbit_count"],
        "combined_nonzero_section_exists": linear["combined_nullity"] > 0,
        "pm1_combined_section_exists": boolean["simultaneously_invariant_and_anti_pm1_sections"] > 0,
        "new_strict_signed_value_supplied": False,
    }
    return {
        "facts": facts,
        "accepted_as_lambda_p2721_source": facts["combined_nonzero_section_exists"] and facts["new_strict_signed_value_supplied"],
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The quotient/no-section calculation proves that quotient descent and anti-equivariant sign response are incompatible on the current sign torsor data.  This supplies a theorem-level obstruction, not a new strict signed value.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    obstruction = payload["section_obstruction"]
    linear = obstruction["linear_system"]
    boolean = obstruction["boolean_cross_check"]
    lines = [
        "# P2739/S1689 sign-torsor quotient no-section certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Linear section obstruction",
        f"- state_count={obstruction['state_count']}",
        f"- orbit_count={obstruction['orbit_count']}",
        f"- invariant_nullity={linear['invariant_quotient_descent_nullity']}",
        f"- anti_equivariant_nullity={linear['anti_equivariant_sign_line_nullity']}",
        f"- combined_rank={linear['combined_rank']}",
        f"- combined_nullity={linear['combined_nullity']}",
        "",
        "## Boolean cross-check",
        f"- total_pm1_sections={boolean['total_pm1_sections']}",
        f"- invariant_pm1_sections={boolean['invariant_pm1_sections']}",
        f"- anti_equivariant_pm1_sections={boolean['anti_equivariant_pm1_sections']}",
        f"- simultaneously_invariant_and_anti_pm1_sections={boolean['simultaneously_invariant_and_anti_pm1_sections']}",
        "",
        "## Theorem statement",
        obstruction["theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    obstruction = section_obstruction()
    acceptance = acceptance_matrix(scan, obstruction)
    payload = {
        "status": "P2739_SIGN_TORSOR_QUOTIENT_NO_SECTION_CERTIFICATE_NO_GO" if not acceptance["accepted_as_lambda_p2721_source"] else "P2739_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "quotient-descending sections of the anti-equivariant sign line over the simultaneous-flip quotient of the four existing sign torsors",
        "content_evidence_scan": scan,
        "section_obstruction": obstruction,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue the existing sign-torsor lane by looking for a quotient section: P2739 proves that quotient descent and anti-equivariant sign response have zero common nonzero section on the current 16-state torsor.  The next proof-grade move must supply an external-to-this-lane but still strict/internal signed value with a source theorem and a P2721 polarity-coupling theorem, or pivot to a genuinely different typed object; otherwise preserve the P2697-P2739 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2739/S1689 sign-torsor quotient no-section certificate", "## P2739/S1689 sign-torsor quotient no-section certificate\n\n`P2739/S1689` strengthens P2738 by replacing Boolean-law enumeration with a finite section theorem.  On the 16 sign states of `lambda`, orientation, `P2721` polarity, and branch-square flux, simultaneous flip gives 8 free orbits.  Quotient descent requires `s(x)=s(Fx)`, while an anti-equivariant sign-line response requires `s(x)=-s(Fx)`.  The combined rational linear system has rank `16` and nullity `0`, and the `{±1}` cross-check finds zero simultaneous sections.  Thus the current sign-torsor quotient has no nonzero section fixing `lambda/P2721`; no strict signed source value, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2739/S1689 quotient no-section Ltotal guard", "## P2739/S1689 quotient no-section Ltotal guard\n\n`P2739/S1689` is a no-section obstruction for the current sign-torsor quotient, not a variational source construction.  It adds no signed source term and does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current sign-torsor quotient no-section guardrail (P2739/S1689, 2026-06-14)", "## Current sign-torsor quotient no-section guardrail (P2739/S1689, 2026-06-14)\n\n- P2739 strengthens P2738 with a finite section theorem on the 16-state sign torsor: simultaneous flip gives 8 free orbits; quotient descent requires `s(x)=s(Fx)`, while anti-equivariant sign-line response requires `s(x)=-s(Fx)`.\n- The combined rational linear system has rank `16` and nullity `0`, and the `{±1}` cross-check finds zero simultaneous sections.  Therefore the current sign-torsor quotient exports no nonzero section fixing `lambda/P2721`, no `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.\n- Do not continue the existing sign-torsor lane by searching for quotient sections or recombinations.  A next admissible move must supply a genuinely new strict/internal signed value with a source theorem and a `P2721` polarity-coupling theorem, or a different typed object outside this closed sign-torsor lane; otherwise preserve the P2697-P2739 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
