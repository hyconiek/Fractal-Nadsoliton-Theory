#!/usr/bin/env python3
"""P2740/S1690: Z12 source-triple chirality orbit no-go.

After P2738/P2739 closed recombinations/sections of existing sign torsors, this
tries a different typed object: a signed cyclic-order chirality on ordered
triples of distinct Z12 source labels.  The value is genuinely nonzero on each
ordered triple, but the finite orbit audit checks whether it descends to a
strict, source-localized, non-premise lambda/P2721 polarity source.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any, Iterable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2740_s1690_z12_source_triple_chirality_orbit_no_go.json"
MD = GEN / "p2740_s1690_z12_source_triple_chirality_orbit_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P2739_NO_SECTION": GEN / "p2739_s1689_sign_torsor_quotient_no_section_certificate.json",
}

CONTENT_PATTERNS = {
    "existing_three_cycle_reference_boundary": r"negative 3-cycle|3-cycle witness|3-cycle sign-holonomy|triangle product",
    "strict_source_triple_boundary": r"strict-source triple|source theorem|strict signed value|signed source",
    "source_orbit_replay_boundary": r"source labels/orbits|source orbit|Aut\(Z12\).*orbit|translation-gauge|translation orbit",
    "p2721_polarity_boundary": r"P2721 polarity|lambda/P2721|polarity-coupling|polarity pair",
}
NEGATIVE_EXPORT_FLAGS = [
    "source_triple_chirality_exported_as_strict_source",
    "source_localizer_exported",
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


def cyclic_chirality(triple: tuple[int, int, int]) -> int:
    """Return +1 if b then c appear counterclockwise after a on Z12, else -1."""
    a, b, c = triple
    if len({a, b, c}) != 3:
        raise ValueError(f"triple must be distinct: {triple}")
    b_rel = (b - a) % 12
    c_rel = (c - a) % 12
    return 1 if 0 < b_rel < c_rel < 12 else -1


def translate_ordered(triple: tuple[int, int, int], shift: int) -> tuple[int, int, int]:
    return tuple((x + shift) % 12 for x in triple)  # type: ignore[return-value]


def affine_ordered(triple: tuple[int, int, int], unit: int, shift: int) -> tuple[int, int, int]:
    return tuple((unit * x + shift) % 12 for x in triple)  # type: ignore[return-value]


def canonical_unordered(triple: Iterable[int]) -> tuple[int, int, int]:
    return tuple(sorted(triple))  # type: ignore[return-value]


def translation_unordered_orbits(unordered: list[tuple[int, int, int]]) -> list[set[tuple[int, int, int]]]:
    remaining = set(unordered)
    orbits = []
    while remaining:
        seed = next(iter(remaining))
        orb = {canonical_unordered(translate_ordered(seed, shift)) for shift in Z12}
        orbits.append(orb)
        remaining -= orb
    return orbits


def affine_unordered_orbits(unordered: list[tuple[int, int, int]]) -> list[set[tuple[int, int, int]]]:
    remaining = set(unordered)
    orbits = []
    while remaining:
        seed = next(iter(remaining))
        orb = {canonical_unordered(affine_ordered(seed, unit, shift)) for unit in UNITS for shift in Z12}
        orbits.append(orb)
        remaining -= orb
    return orbits


def signed_sum_for_unordered_family(family: Iterable[tuple[int, int, int]]) -> int:
    total = 0
    for unordered in family:
        for ordered in itertools.permutations(unordered, 3):
            total += cyclic_chirality(ordered)  # three + and three - per unordered triple
    return total


def triple_chirality_audit() -> dict[str, Any]:
    ordered = list(itertools.permutations(Z12, 3))
    unordered = list(itertools.combinations(Z12, 3))
    signs = [cyclic_chirality(t) for t in ordered]
    translation_orbits = translation_unordered_orbits(unordered)
    affine_orbits = affine_unordered_orbits(unordered)
    translation_rows = [
        {"orbit_size": len(orb), "signed_sum": signed_sum_for_unordered_family(orb), "sample": sorted(orb)[:4]}
        for orb in translation_orbits
    ]
    affine_rows = [
        {"orbit_size": len(orb), "signed_sum": signed_sum_for_unordered_family(orb), "sample": sorted(orb)[:4]}
        for orb in affine_orbits
    ]
    return {
        "typed_object": "Z12 ordered source-triple cyclic chirality chi(a,b,c)=sign of cyclic order after translating a to 0",
        "ordered_distinct_triples": len(ordered),
        "unordered_triples": len(unordered),
        "positive_ordered_triples": signs.count(1),
        "negative_ordered_triples": signs.count(-1),
        "translation_unordered_orbit_count": len(translation_orbits),
        "affine_unordered_orbit_count": len(affine_orbits),
        "translation_orbits_with_nonzero_signed_sum": sum(1 for row in translation_rows if row["signed_sum"] != 0),
        "affine_orbits_with_nonzero_signed_sum": sum(1 for row in affine_rows if row["signed_sum"] != 0),
        "translation_orbit_rows": translation_rows,
        "affine_orbit_rows": affine_rows,
        "theorem": "The Z12 ordered-triple chirality is nonzero on every ordered distinct triple and is balanced globally (660 positive, 660 negative).  After forgetting source labels to translation or affine unordered orbits, every orbit has signed sum zero, because each unordered triple contains three positive and three negative orderings.  Thus the object supplies a real chiral sign only after choosing an ordered source triple; current artifacts provide no strict source localizer or P2721 polarity coupling for that choice.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_related_boundaries": scan["all_patterns_have_hits"],
        "candidate_has_nonzero_pointwise_signed_values": audit["positive_ordered_triples"] > 0 and audit["negative_ordered_triples"] > 0,
        "translation_orbit_signed_source_survives": audit["translation_orbits_with_nonzero_signed_sum"] > 0,
        "affine_orbit_signed_source_survives": audit["affine_orbits_with_nonzero_signed_sum"] > 0,
        "strict_source_localizer_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_signed_source": all(facts.values()),
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The ordered-triple chirality has nonzero pointwise signs, but orbit-safe aggregation cancels and no strict source-localizer chooses one ordered triple or couples it to one P2721 polarity.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["triple_chirality_audit"]
    lines = [
        "# P2740/S1690 Z12 source-triple chirality orbit no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite orbit audit",
        f"- ordered_distinct_triples={audit['ordered_distinct_triples']}",
        f"- unordered_triples={audit['unordered_triples']}",
        f"- positive_ordered_triples={audit['positive_ordered_triples']}",
        f"- negative_ordered_triples={audit['negative_ordered_triples']}",
        f"- translation_unordered_orbit_count={audit['translation_unordered_orbit_count']}",
        f"- affine_unordered_orbit_count={audit['affine_unordered_orbit_count']}",
        f"- translation_orbits_with_nonzero_signed_sum={audit['translation_orbits_with_nonzero_signed_sum']}",
        f"- affine_orbits_with_nonzero_signed_sum={audit['affine_orbits_with_nonzero_signed_sum']}",
        "",
        "## Theorem statement",
        audit["theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    audit = triple_chirality_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2740_Z12_SOURCE_TRIPLE_CHIRALITY_ORBIT_NO_GO" if not acceptance["accepted_as_strict_signed_source"] else "P2740_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "Z12 ordered source-triple cyclic chirality as a new signed typed object outside the closed sign-torsor quotient/recombination class",
        "content_evidence_scan": scan,
        "triple_chirality_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote ordered-triple chirality by itself: it is a real nonzero pointwise sign, but without a strict source-localizer for one ordered triple and an exported P2721 polarity-coupling theorem it remains label/order premise data.  The next proof-grade move must either construct that localizer-and-coupling theorem for this exact triple-chirality object, or pivot to a different strict signed observable with nonzero orbit-safe signed value; otherwise preserve the P2697-P2740 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2740/S1690 Z12 source-triple chirality orbit no-go", "## P2740/S1690 Z12 source-triple chirality orbit no-go\n\n`P2740/S1690` tests a genuinely new signed typed object after P2739: cyclic-order chirality on ordered triples of distinct `Z12` source labels.  The sign is nonzero pointwise (`660` positive and `660` negative ordered triples), but every translation or affine unordered source-orbit has signed sum `0`; the sign survives only after choosing an ordered source triple.  Current artifacts export no strict source-localizer for that choice and no `P2721` polarity-coupling theorem.  Therefore no `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2740/S1690 source-triple chirality Ltotal guard", "## P2740/S1690 source-triple chirality Ltotal guard\n\n`P2740/S1690` supplies a bounded no-go for an ordered source-triple chirality candidate.  Because the orbit-safe signed sums cancel and no strict source-localizer/coupling theorem is exported, it adds no variational source term and does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Z12 source-triple chirality orbit no-go guardrail (P2740/S1690, 2026-06-14)", "## Current Z12 source-triple chirality orbit no-go guardrail (P2740/S1690, 2026-06-14)\n\n- P2740 tests a new signed typed object outside the closed sign-torsor recombination/quotient class: cyclic-order chirality on ordered triples of distinct `Z12` source labels.\n- The sign is nonzero pointwise (`660` positive and `660` negative ordered triples), but all translation and affine unordered source-orbit signed sums are zero; the sign requires an ordered source-triple localizer.\n- Do not promote source-triple chirality alone to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must construct a strict source-localizer plus `P2721` polarity-coupling theorem for this exact object, or pivot to a different strict signed observable with nonzero orbit-safe signed value; otherwise preserve the P2697-P2740 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
