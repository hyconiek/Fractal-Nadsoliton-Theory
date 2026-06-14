#!/usr/bin/env python3
"""P2738/S1688: finite sign-torsor Boolean source-law exhaustion.

Next honest step after P2737.  Before proposing another ``strict signed source
law'', this script greps the repository by research content and then exhausts a
finite class of candidate laws built only from already-audited sign-like inputs:
lambda sign, orientation sign, P2721 polarity sign, and branch-square flux sign.

The point is narrow: if all available inputs are themselves unfixed sign torsors,
then a Boolean law can at best be equivariant under the simultaneous sign flip;
it cannot produce a non-premise absolute polarity on the quotient.  This is a
proof-grade no-go for "make the missing sign out of the existing missing signs".
"""
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
OUT = GEN / "p2738_s1688_sign_torsor_boolean_source_law_exhaustion.json"
MD = GEN / "p2738_s1688_sign_torsor_boolean_source_law_exhaustion.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P2737_TOE_READINESS": GEN / "p2737_s1687_lay_toe_potential_readiness_matrix.json",
}

CONTENT_PATTERNS = {
    "unfixed_lambda_sign": r"lambda remains unfixed|lambda.*unfixed|strict law fixing `lambda`|fixing lambda",
    "p2721_polarity_pair": r"P2721 polarity|P2721.*polarity|coupling polarity|polarity pair",
    "branch_square_flux_no_source": r"branch-square.*flux|non-exact flux|frustrated flux|plaquette holonomy",
    "time_or_chiral_sign_pair": r"tau.*sign|time-arrow.*source|chiral.*signed|signed source law",
    "global_flip_or_sign_pairing": r"sign-paired|lambda -> -lambda|\+omega.*-omega|opposite polarity|two polarity",
}

NEGATIVE_EXPORT_FLAGS = [
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


def sign_inputs() -> list[tuple[int, int, int, int]]:
    return list(itertools.product((-1, 1), repeat=4))


def flip(x: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(-v for v in x)  # simultaneous torsor reversal


def law_output(mask: int, index: int) -> int:
    return 1 if (mask >> index) & 1 else -1


def classify_laws() -> dict[str, Any]:
    inputs = sign_inputs()
    index = {x: i for i, x in enumerate(inputs)}
    total_laws = 2 ** len(inputs)
    equivariant = []
    invariant = []
    quotient_constant_positive = []
    quotient_constant_negative = []
    orbit_rows = []

    seen_orbits = set()
    for x in inputs:
        orb = tuple(sorted((x, flip(x))))
        if orb not in seen_orbits:
            seen_orbits.add(orb)
            orbit_rows.append({"representative": x, "flipped": flip(x)})

    for mask in range(total_laws):
        is_equivariant = True
        is_invariant = True
        constant_plus_on_orbits = True
        constant_minus_on_orbits = True
        for x in inputs:
            y = flip(x)
            fx = law_output(mask, index[x])
            fy = law_output(mask, index[y])
            if fy != -fx:
                is_equivariant = False
            if fy != fx:
                is_invariant = False
            if fx != 1 or fy != 1:
                constant_plus_on_orbits = False
            if fx != -1 or fy != -1:
                constant_minus_on_orbits = False
        if is_equivariant:
            equivariant.append(mask)
        if is_invariant:
            invariant.append(mask)
        if constant_plus_on_orbits:
            quotient_constant_positive.append(mask)
        if constant_minus_on_orbits:
            quotient_constant_negative.append(mask)

    accepted_nonpremise = []
    for mask in equivariant:
        # Acceptance would require an odd/equivariant response to sign reversal
        # AND a well-defined absolute + polarity on every quotient orbit.  These
        # requirements are incompatible unless a new fixed sign is supplied.
        if mask in quotient_constant_positive:
            accepted_nonpremise.append(mask)

    return {
        "variables": ["lambda_sign", "orientation_sign", "P2721_polarity_sign", "branch_square_flux_sign"],
        "input_state_count": len(inputs),
        "global_flip_orbit_count": len(orbit_rows),
        "total_boolean_laws": total_laws,
        "equivariant_odd_law_count": len(equivariant),
        "invariant_even_law_count": len(invariant),
        "absolute_plus_law_count": len(quotient_constant_positive),
        "absolute_minus_law_count": len(quotient_constant_negative),
        "accepted_nonpremise_law_count": len(accepted_nonpremise),
        "sample_equivariant_masks": equivariant[:8],
        "orbit_rows": orbit_rows,
        "theorem": "Across all 2^16 Boolean laws on the four existing sign torsors, the odd/equivariant laws come in flip-paired outputs on the 8 quotient orbits, while the only absolute constant laws are not odd/equivariant. Therefore existing unfixed sign torsors cannot manufacture a non-premise absolute lambda/P2721 polarity.",
    }


def acceptance_matrix(exhaustion: dict[str, Any], scan: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_confirms_existing_sign_torsor_boundaries": scan["all_patterns_have_hits"],
        "all_boolean_laws_exhausted": exhaustion["total_boolean_laws"] == 65536,
        "odd_equivariant_laws_exist": exhaustion["equivariant_odd_law_count"] > 0,
        "nonpremise_absolute_polarity_law_found": exhaustion["accepted_nonpremise_law_count"] > 0,
        "new_strict_signed_source_supplied": False,
        "lambda_p2721_fixed": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_signed_source_law": facts["nonpremise_absolute_polarity_law_found"] and facts["new_strict_signed_source_supplied"],
        "blocker": "Finite exhaustion finds equivariant sign-response laws, but every such law remains paired under simultaneous torsor reversal. No new strict signed value is supplied, so no absolute lambda/P2721 polarity is selected.",
        "missing_criteria": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    ex = payload["boolean_source_law_exhaustion"]
    lines = [
        "# P2738/S1688 sign-torsor Boolean source-law exhaustion",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite exhaustion",
        f"- input_state_count={ex['input_state_count']}",
        f"- global_flip_orbit_count={ex['global_flip_orbit_count']}",
        f"- total_boolean_laws={ex['total_boolean_laws']}",
        f"- equivariant_odd_law_count={ex['equivariant_odd_law_count']}",
        f"- accepted_nonpremise_law_count={ex['accepted_nonpremise_law_count']}",
        "",
        "## Theorem statement",
        ex["theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    exhaustion = classify_laws()
    acceptance = acceptance_matrix(exhaustion, scan)
    payload = {
        "status": "P2738_SIGN_TORSOR_BOOLEAN_SOURCE_LAW_EXHAUSTION_NO_GO" if not acceptance["accepted_as_strict_signed_source_law"] else "P2738_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all Boolean source laws from the already-known unfixed sign torsors lambda, orientation, P2721 polarity, and branch-square flux, modulo simultaneous sign reversal",
        "content_evidence_scan": scan,
        "boolean_source_law_exhaustion": exhaustion,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not try to build the missing strict signed source law only by recombining existing unfixed sign torsors.  P2738 exhausts that finite Boolean class and finds no non-premise absolute lambda/P2721 polarity.  The next proof-grade move must introduce one genuinely new strict signed value with an exported source theorem and an explicit coupling to exactly one P2721 polarity, or pivot to a different typed object outside the sign-torsor replay class; otherwise preserve the P2697-P2738 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2738/S1688 sign-torsor Boolean source-law exhaustion", "## P2738/S1688 sign-torsor Boolean source-law exhaustion\n\n`P2738/S1688` exhausts all `2^16 = 65536` Boolean laws on four already-known unfixed sign torsors: `lambda`, orientation, `P2721` polarity, and branch-square flux.  There are odd/equivariant sign-response laws, but on the eight simultaneous-flip quotient orbits every such law remains paired; the absolute constant laws are not odd/equivariant and therefore only impose a premise.  No new strict signed source value, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2738/S1688 sign-torsor source-law Ltotal guard", "## P2738/S1688 sign-torsor source-law Ltotal guard\n\n`P2738/S1688` is a finite no-go for constructing the missing source by recombining existing unfixed sign torsors.  It adds no variational source term and does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current sign-torsor Boolean source-law exhaustion guardrail (P2738/S1688, 2026-06-14)", "## Current sign-torsor Boolean source-law exhaustion guardrail (P2738/S1688, 2026-06-14)\n\n- P2738 exhausts all `2^16 = 65536` Boolean laws built from the already-known unfixed sign torsors `lambda`, orientation, `P2721` polarity, and branch-square flux, after content-first grep confirms these are existing bounded/replay lanes.\n- Odd/equivariant laws exist, but they remain paired on the simultaneous-flip quotient; absolute constant laws merely impose a premise and are not strict sourced sign laws.  No non-premise `lambda/P2721` polarity, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n- Do not replay recombinations of existing unfixed sign torsors as the missing strict signed source.  A next admissible move must introduce one genuinely new strict signed value with an exported source theorem and explicit coupling to exactly one `P2721` polarity, or a different typed object outside this sign-torsor replay class; otherwise preserve the P2697-P2738 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
