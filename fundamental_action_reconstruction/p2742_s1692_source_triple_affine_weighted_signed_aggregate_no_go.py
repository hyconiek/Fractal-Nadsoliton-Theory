#!/usr/bin/env python3
"""P2742/S1692: affine-weighted source-triple signed aggregate no-go.

P2741 blocked a fixed ordered-triple localizer.  This follow-up tests the more
global escape hatch: do not localize a triple, but assign any affine-invariant
orbit weights to the P2740 ordered-triple chirality and ask whether a nonzero
orbit-safe signed aggregate survives.
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
OUT = GEN / "p2742_s1692_source_triple_affine_weighted_signed_aggregate_no_go.json"
MD = GEN / "p2742_s1692_source_triple_affine_weighted_signed_aggregate_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {"P2741_LOCALIZER_NO_GO": GEN / "p2741_s1691_source_triple_localizer_fixed_point_no_go.json"}

CONTENT_PATTERNS = {
    "p2741_localizer_no_go": r"P2741|fixed ordered triple|singleton orbit|source-triple localizer",
    "affine_orbit_weighting_boundary": r"affine ordered orbits|affine unordered source-orbit|Aut\(Z12\).*orbit|orbit-safe signed value|orbit weighting",
    "p2740_chirality_boundary": r"ordered-triple chirality|cyclic-order chirality|source-triple chirality|ordered triples",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "affine_weighted_source_triple_signed_value_exported",
    "orbit_safe_signed_value_nonzero",
    "strict_source_localizer_exported",
    "p2721_polarity_selected",
    "lambda_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]
Z12 = tuple(range(12))
UNITS = (1, 5, 7, 11)
ORIENTATION_REVERSING_UNITS = (5, 7, 11)  # any unit not preserving cyclic order globally; use 11 as canonical involution witness


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


def affine(triple: tuple[int, int, int], unit: int, shift: int) -> tuple[int, int, int]:
    return tuple((unit * x + shift) % 12 for x in triple)  # type: ignore[return-value]


def cyclic_chirality(triple: tuple[int, int, int]) -> int:
    a, b, c = triple
    if len({a, b, c}) != 3:
        raise ValueError(f"triple must be distinct: {triple}")
    b_rel = (b - a) % 12
    c_rel = (c - a) % 12
    return 1 if 0 < b_rel < c_rel < 12 else -1


def affine_orbit(seed: tuple[int, int, int]) -> set[tuple[int, int, int]]:
    return {affine(seed, unit, shift) for unit in UNITS for shift in Z12}


def affine_ordered_orbits(ordered: list[tuple[int, int, int]]) -> list[set[tuple[int, int, int]]]:
    remaining = set(ordered)
    out = []
    while remaining:
        seed = next(iter(remaining))
        orb = affine_orbit(seed)
        out.append(orb)
        remaining -= orb
    return out


def aggregate_audit() -> dict[str, Any]:
    ordered = list(itertools.permutations(Z12, 3))
    orbits = affine_ordered_orbits(ordered)
    rows = []
    for idx, orb in enumerate(orbits):
        plus = sum(1 for triple in orb if cyclic_chirality(triple) == 1)
        minus = sum(1 for triple in orb if cyclic_chirality(triple) == -1)
        signed_sum = plus - minus
        witness_ok = all(cyclic_chirality(affine(triple, 11, 0)) == -cyclic_chirality(triple) for triple in orb)
        rows.append(
            {
                "orbit_index": idx,
                "orbit_size": len(orb),
                "positive_chirality_count": plus,
                "negative_chirality_count": minus,
                "signed_sum_coefficient": signed_sum,
                "unit_11_pairs_chirality_opposites": witness_ok,
                "representative": sorted(orb)[0],
            }
        )
    coefficient_vector = [row["signed_sum_coefficient"] for row in rows]
    nonzero_coefficients = [value for value in coefficient_vector if value != 0]
    return {
        "typed_escape_hatch": "arbitrary affine-invariant orbit weights on ordered source-triple chirality",
        "ordered_distinct_triples": len(ordered),
        "affine_group_size": 48,
        "affine_ordered_orbit_count": len(orbits),
        "affine_orbit_sizes": sorted({len(orb) for orb in orbits}),
        "orbits_with_nonzero_signed_sum_coefficient": len(nonzero_coefficients),
        "signed_sum_linear_map_rank": 0 if not nonzero_coefficients else 1,
        "signed_sum_linear_map_nullity_over_orbit_weights": len(orbits) if not nonzero_coefficients else len(orbits) - 1,
        "one_hot_weight_crosscheck_nonzero_count": sum(1 for value in coefficient_vector if value != 0),
        "all_unit_11_pairing_witnesses_pass": all(row["unit_11_pairs_chirality_opposites"] for row in rows),
        "orbit_rows": rows,
        "finite_theorem": "Every affine ordered orbit of distinct Z12 triples is paired by the inversion unit 11 into opposite cyclic chiralities, so each orbit has signed-sum coefficient zero.  Therefore any affine-invariant orbit weighting, over any coefficient field of characteristic not 2, has total signed aggregate zero.  The P2740 chirality cannot be made into a nonzero orbit-safe strict signed value by affine orbit weights after P2741 blocks fixed localizers.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_weighting_escape_hatch": scan["all_patterns_have_hits"],
        "candidate_domain_is_nonempty": audit["ordered_distinct_triples"] == 1320,
        "affine_orbit_weighted_signed_aggregate_nonzero": audit["orbits_with_nonzero_signed_sum_coefficient"] > 0,
        "strict_source_localizer_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_orbit_safe_signed_source": all(facts.values()),
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "All affine ordered source-triple orbits have zero chirality signed-sum coefficient, so arbitrary affine-invariant weights cannot produce a nonzero orbit-safe signed source or fix lambda/P2721.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["aggregate_audit"]
    lines = [
        "# P2742/S1692 source-triple affine-weighted signed aggregate no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite aggregate audit",
        f"- ordered_distinct_triples={audit['ordered_distinct_triples']}",
        f"- affine_group_size={audit['affine_group_size']}",
        f"- affine_ordered_orbit_count={audit['affine_ordered_orbit_count']}",
        f"- affine_orbit_sizes={audit['affine_orbit_sizes']}",
        f"- orbits_with_nonzero_signed_sum_coefficient={audit['orbits_with_nonzero_signed_sum_coefficient']}",
        f"- signed_sum_linear_map_rank={audit['signed_sum_linear_map_rank']}",
        f"- signed_sum_linear_map_nullity_over_orbit_weights={audit['signed_sum_linear_map_nullity_over_orbit_weights']}",
        f"- one_hot_weight_crosscheck_nonzero_count={audit['one_hot_weight_crosscheck_nonzero_count']}",
        f"- all_unit_11_pairing_witnesses_pass={audit['all_unit_11_pairing_witnesses_pass']}",
        "",
        "## Theorem statement",
        audit["finite_theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    audit = aggregate_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2742_SOURCE_TRIPLE_AFFINE_WEIGHTED_SIGNED_AGGREGATE_NO_GO" if not acceptance["accepted_as_orbit_safe_signed_source"] else "P2742_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "affine-invariant orbit-weighted aggregate of the exact P2740 source-triple cyclic chirality after P2741 localizer no-go",
        "content_evidence_scan": scan,
        "aggregate_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue the P2740/P2741 source-triple lane by affine orbit weighting: P2742 proves every affine ordered orbit has zero chirality signed-sum coefficient, so arbitrary affine-invariant weights still give total signed aggregate zero.  The next proof-grade move must pivot to a genuinely different strict signed observable with a computable nonzero orbit-safe signed value and a P2721 coupling theorem, or else preserve the P2697-P2742 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2742/S1692 source-triple affine-weighted signed aggregate no-go", "## P2742/S1692 source-triple affine-weighted signed aggregate no-go\n\n`P2742/S1692` tests the global escape hatch left after `P2741/S1691`: instead of localizing one ordered triple, assign arbitrary affine-invariant weights to affine ordered source-triple orbits.  The finite computation finds `34` affine ordered orbits with signed-sum coefficient `0` on every orbit, witnessed by the inversion unit `11` pairing opposite cyclic chiralities.  The signed aggregate linear map has rank `0` and nullity `34`, so no affine orbit-weighted nonzero strict signed value, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2742/S1692 affine-weighted source-triple Ltotal guard", "## P2742/S1692 affine-weighted source-triple Ltotal guard\n\n`P2742/S1692` adds no variational source term: arbitrary affine-invariant orbit weights on the `P2740` ordered-triple chirality have zero total signed aggregate because every affine ordered orbit has zero signed-sum coefficient.  Since no nonzero orbit-safe signed value and no `P2721` coupling theorem are exported, this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current source-triple affine-weighted signed aggregate no-go guardrail (P2742/S1692, 2026-06-14)", "## Current source-triple affine-weighted signed aggregate no-go guardrail (P2742/S1692, 2026-06-14)\n\n- P2742 tests the global escape hatch left after P2741: arbitrary affine-invariant weights on affine ordered source-triple orbits for the P2740 cyclic chirality.\n- The finite computation finds `34` affine ordered orbits and zero signed-sum coefficient on every orbit; the inversion unit `11` pairs each orbit into opposite chiralities, so the signed aggregate map from orbit weights has rank `0` and nullity `34`.\n- Do not continue the source-triple chirality lane by affine orbit weighting or representative choice as `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must pivot to a genuinely different strict signed observable with a computable nonzero orbit-safe signed value plus `P2721` coupling theorem, or preserve the P2697-P2742 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
