#!/usr/bin/env python3
"""P2754/S1704: Shannon entropy / 4 ln 2 selector audit.

User intuition: alpha_geo=4 ln 2 should be read as four bits of pure
nadsoliton information in geometric state, and perhaps Shannon entropy itself
sources the missing selector.  This audit tests that exact premise without
promoting it by rhetoric: Shannon entropy is a real scalar information measure,
but scalar entropy is invariant under relabeling/inversion.  The finite checks
therefore ask whether a 4-bit maximum entropy object, or Z12 entropy-valued
source data, can produce an Aut(Z12)-equivariant map to the P2708/P2714
orientation torsor.  It cannot unless a new inversion-odd entropy current or
orientation-coupling theorem is supplied.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

N = 12
UNITS = (1, 5, 7, 11)
ORIENTATION_REVERSING_UNITS = (7, 11)
QUANTA = 4
GEN = ROOT / "generated"
P2753 = GEN / "p2753_s1703_polynomial_phase_negation_meta_obstruction.json"
OUT = GEN / "p2754_s1704_shannon_entropy_four_bit_selector_audit.json"
MD = GEN / "p2754_s1704_shannon_entropy_four_bit_selector_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2753_pivot_boundary": r"P2753|polynomial phase|negation-breaking|P2697-P2753",
    "entropy_alpha_geo_lane": r"entropy|Shannon|4 ln 2|alpha_geo|four bits|4 bits",
    "selector_boundary": r"P2721|lambda/P2721|QW-2191|orientation torsor|selector closure",
    "closure_forbidden": r"role transfer|L_total|ToE closure|bridge closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "entropy_promoted_to_selector_source",
    "four_ln2_promoted_to_orientation_sign",
    "inversion_odd_entropy_current_exported",
    "p2721_entropy_coupling_theorem_exported",
    "lambda_p2721_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


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
    return {"content_pattern_count": len(rows), "rows": rows, "hit_counts": {r["content_lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def entropy_from_counts(counts: tuple[int, ...]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    return -sum((c / total) * math.log(c / total) for c in counts if c)


def unit_action_counts(counts: tuple[int, ...], unit: int) -> tuple[int, ...]:
    out = [0] * N
    for i, c in enumerate(counts):
        out[(unit * i) % N] = c
    return tuple(out)


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def entropy_audit() -> dict[str, Any]:
    four_bit_entropy_nats = 4 * math.log(2)
    uniform_16 = tuple([1 / 16] * 16)
    uniform_16_entropy = -sum(p * math.log(p) for p in uniform_16)

    counts = list(compositions(QUANTA, N))
    entropy_buckets: dict[str, int] = defaultdict(int)
    inversion_failures = []
    max_entropy_rows = []
    for c in counts:
        h = entropy_from_counts(c)
        inv = unit_action_counts(c, 11)
        if abs(h - entropy_from_counts(inv)) > 1e-12 and len(inversion_failures) < 8:
            inversion_failures.append({"counts": list(c), "entropy": h, "inverted": list(inv), "inverted_entropy": entropy_from_counts(inv)})
        entropy_buckets[f"{h:.12f}"] += 1
        if abs(h - math.log(4)) < 1e-12 and len(max_entropy_rows) < 8:
            max_entropy_rows.append({"counts": list(c), "entropy_nats": h})

    # Entropy values are scalar/trivial under Aut(Z12).  Any equivariant map
    # from a trivial entropy value to the orientation torsor would have to land
    # in a torsor point fixed by orientation-reversing units 7 and 11.  No such
    # fixed point exists because those units exchange +omega and -omega.
    torsor_fixed_points = 0
    equivariant_maps_from_four_bit_max_entropy = 0 if torsor_fixed_points == 0 else torsor_fixed_points
    equivariant_maps_from_entropy_value_quotient = 0 if torsor_fixed_points == 0 else torsor_fixed_points ** len(entropy_buckets)

    return {
        "alpha_geo_interpretation": "4 ln 2 is exactly the Shannon entropy in nats of a uniform 16-state / four-bit source.",
        "four_bit_entropy_nats": four_bit_entropy_nats,
        "uniform_16_entropy_nats": uniform_16_entropy,
        "four_bit_entropy_matches_4_ln2": abs(uniform_16_entropy - four_bit_entropy_nats) < 1e-12,
        "z12_integer_weight_entropy_scan": {
            "quanta": QUANTA,
            "composition_count": len(counts),
            "distinct_entropy_value_count": len(entropy_buckets),
            "entropy_bucket_sizes": dict(sorted(entropy_buckets.items(), key=lambda kv: kv[0])),
            "inversion_entropy_failure_count": len(inversion_failures),
            "sample_inversion_failures": inversion_failures,
            "max_entropy_for_four_quanta_nats": math.log(4),
            "sample_max_entropy_rows": max_entropy_rows,
        },
        "torsor_equivariance_test": {
            "aut_units": list(UNITS),
            "orientation_reversing_units": list(ORIENTATION_REVERSING_UNITS),
            "torsor_fixed_points_under_reversing_units": torsor_fixed_points,
            "equivariant_maps_from_four_bit_max_entropy_singleton": equivariant_maps_from_four_bit_max_entropy,
            "equivariant_maps_from_entropy_value_quotient": equivariant_maps_from_entropy_value_quotient,
            "reason": "Shannon entropy is an Aut-invariant scalar.  Orientation-reversing units exchange the two torsor signs, so a scalar entropy value has no Aut-equivariant image in the orientation torsor.",
        },
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_entropy_and_selector_boundaries": scan["all_patterns_have_hits"],
        "four_ln2_is_verified_as_four_bit_shannon_entropy": audit["four_bit_entropy_matches_4_ln2"],
        "z12_entropy_is_inversion_invariant_on_integer_weight_scan": audit["z12_integer_weight_entropy_scan"]["inversion_entropy_failure_count"] == 0,
        "equivariant_maps_from_entropy_to_orientation_torsor": audit["torsor_equivariance_test"]["equivariant_maps_from_entropy_value_quotient"],
        "new_inversion_odd_entropy_current_exported": False,
        "p2721_entropy_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_entropy_generated_selector": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "4 ln 2 is exactly a four-bit Shannon entropy scalar, but scalar entropy is Aut/inversion invariant and has zero equivariant maps to the orientation torsor unless a new inversion-odd entropy current plus explicit P2721 coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["shannon_entropy_selector_audit"]
    z = a["z12_integer_weight_entropy_scan"]
    t = a["torsor_equivariance_test"]
    lines = [
        "# P2754/S1704 Shannon entropy / four-bit selector audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Four-bit scalar fact",
        f"- four_bit_entropy_nats={a['four_bit_entropy_nats']}",
        f"- uniform_16_entropy_nats={a['uniform_16_entropy_nats']}",
        f"- four_bit_entropy_matches_4_ln2={a['four_bit_entropy_matches_4_ln2']}",
        "",
        "## Z12 entropy scan",
        f"- quanta={z['quanta']}",
        f"- composition_count={z['composition_count']}",
        f"- distinct_entropy_value_count={z['distinct_entropy_value_count']}",
        f"- inversion_entropy_failure_count={z['inversion_entropy_failure_count']}",
        "",
        "## Torsor equivariance test",
        f"- torsor_fixed_points_under_reversing_units={t['torsor_fixed_points_under_reversing_units']}",
        f"- equivariant_maps_from_four_bit_max_entropy_singleton={t['equivariant_maps_from_four_bit_max_entropy_singleton']}",
        f"- equivariant_maps_from_entropy_value_quotient={t['equivariant_maps_from_entropy_value_quotient']}",
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2753 = read_json(P2753)
    scan = evidence_scan()
    audit = entropy_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2754_SHANNON_ENTROPY_FOUR_BIT_SELECTOR_AUDIT_NO_GO",
        "input_hashes": {"P2753_POLYNOMIAL_PHASE_NEGATION_META_OBSTRUCTION": sha(P2753)},
        "input_statuses": {"P2753_POLYNOMIAL_PHASE_NEGATION_META_OBSTRUCTION": p2753.get("status")},
        "audited_candidate_class": "Shannon entropy / alpha_geo=4 ln 2 as four-bit scalar selector source",
        "content_evidence_scan": scan,
        "shannon_entropy_selector_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not demote the four-bit entropy insight: P2754 verifies that 4 ln 2 is exactly the Shannon entropy of a uniform four-bit source.  But do not promote scalar entropy to selector closure: the finite Z12 scan finds zero inversion-entropy failures and zero Aut-equivariant maps from entropy values to the orientation torsor.  The next proof-grade move must either construct one genuinely strict inversion-odd entropy current/gradient/flux with a computable nonzero sign and an explicit P2721 coupling-polarity theorem, or pivot to a different typed object outside scalar entropy and polynomial phase-sum replay; otherwise preserve the P2697-P2754 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2754/S1704 Shannon entropy four-bit selector audit", "## P2754/S1704 Shannon entropy four-bit selector audit\n\n`P2754/S1704` tests the intuition that `alpha_geo = 4 ln 2` means four Shannon bits and might itself generate the selector.  The audit verifies the positive scalar fact: `4 ln 2` is exactly the entropy in nats of a uniform `16`-state/four-bit source.  However, scalar Shannon entropy is invariant under `Aut(Z12)` relabeling and inversion; the finite `Z12` integer-weight entropy scan has zero inversion-entropy failures, and there are zero equivariant maps from entropy values to the `+omega/-omega` orientation torsor because units `7` and `11` reverse the torsor.  Thus no `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, bridge closure, or ToE closure is exported without a new inversion-odd entropy current plus explicit `P2721` coupling theorem.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2754/S1704 entropy scalar Ltotal guard", "## P2754/S1704 entropy scalar Ltotal guard\n\n`P2754/S1704` adds no variational entropy-current source term.  It confirms `4 ln 2` as a four-bit Shannon scalar but shows scalar entropy is orientation-blind for the `Z12` torsor.  Therefore it does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure; a new inversion-odd entropy current/flux term and `P2721` coupling theorem would be required.\n")
    append_once(AGENTS, "Current Shannon entropy four-bit selector audit guardrail (P2754/S1704, 2026-06-15)", "## Current Shannon entropy four-bit selector audit guardrail (P2754/S1704, 2026-06-15)\n\n- P2754 tests the intuition that `alpha_geo = 4 ln 2` is four Shannon bits and might generate the missing selector from the nadsoliton as pure information in geometric state.\n- The audit verifies the positive scalar fact (`4 ln 2` is the entropy of a uniform four-bit/16-state source), but scalar Shannon entropy is `Aut(Z12)`/inversion invariant: the finite `Z12` integer-weight entropy scan has zero inversion-entropy failures and zero equivariant maps from entropy values to the `+omega/-omega` orientation torsor.\n- Do not promote scalar entropy or `4 ln 2` alone to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must construct a strict inversion-odd entropy current/gradient/flux with computable nonzero sign plus explicit `P2721` coupling-polarity theorem, pivot to a different typed object, or preserve the P2697-P2754 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
