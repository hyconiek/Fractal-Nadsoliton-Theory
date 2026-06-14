#!/usr/bin/env python3
"""P2741/S1691: source-triple localizer fixed-point no-go.

P2740 found a pointwise nonzero cyclic chirality on ordered triples, but it
requires a strict localizer for one ordered source triple.  This follow-up tests
that exact missing premise: can translation/affine-safe Z12 data select one
ordered distinct triple without importing a label/order convention?
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
OUT = GEN / "p2741_s1691_source_triple_localizer_fixed_point_no_go.json"
MD = GEN / "p2741_s1691_source_triple_localizer_fixed_point_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {"P2740_TRIPLE_CHIRALITY": GEN / "p2740_s1690_z12_source_triple_chirality_orbit_no_go.json"}

CONTENT_PATTERNS = {
    "p2740_localizer_obligation": r"ordered source-triple localizer|source-localizer|strict source-localizer|ordered-triple chirality",
    "translation_affine_gauge_boundary": r"translation-gauge|translation orbit|affine unordered source-orbit|Aut\(Z12\).*orbit|source labels/orbits",
    "p2721_coupling_boundary": r"P2721 polarity|lambda/P2721|polarity-coupling|coupling theorem",
    "closure_blocks": r"QW-2191|selector closure|role transfer|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "strict_source_triple_localizer_exported",
    "ordered_triple_selected_nonpremise",
    "source_triple_chirality_exported_as_strict_source",
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


def translate(triple: tuple[int, int, int], shift: int) -> tuple[int, int, int]:
    return tuple((x + shift) % 12 for x in triple)  # type: ignore[return-value]


def affine(triple: tuple[int, int, int], unit: int, shift: int) -> tuple[int, int, int]:
    return tuple((unit * x + shift) % 12 for x in triple)  # type: ignore[return-value]


def orbit(seed: tuple[int, int, int], affine_safe: bool) -> set[tuple[int, int, int]]:
    if affine_safe:
        return {affine(seed, u, s) for u in UNITS for s in Z12}
    return {translate(seed, s) for s in Z12}


def orbit_partition(triples: list[tuple[int, int, int]], affine_safe: bool) -> list[set[tuple[int, int, int]]]:
    remaining = set(triples)
    out = []
    while remaining:
        seed = next(iter(remaining))
        orb = orbit(seed, affine_safe)
        out.append(orb)
        remaining -= orb
    return out


def fixed_by_all_translations(triple: tuple[int, int, int]) -> bool:
    return all(translate(triple, s) == triple for s in Z12)


def fixed_by_all_affine(triple: tuple[int, int, int]) -> bool:
    return all(affine(triple, u, s) == triple for u in UNITS for s in Z12)


def localizer_audit() -> dict[str, Any]:
    ordered = list(itertools.permutations(Z12, 3))
    translation_orbits = orbit_partition(ordered, affine_safe=False)
    affine_orbits = orbit_partition(ordered, affine_safe=True)
    translation_fixed = [t for t in ordered if fixed_by_all_translations(t)]
    affine_fixed = [t for t in ordered if fixed_by_all_affine(t)]
    translation_singletons = [orb for orb in translation_orbits if len(orb) == 1]
    affine_singletons = [orb for orb in affine_orbits if len(orb) == 1]
    return {
        "candidate_missing_premise": "strict source-localizer selecting one ordered distinct Z12 source triple for the P2740 chirality",
        "ordered_distinct_triples": len(ordered),
        "translation_group_size": 12,
        "affine_group_size": 48,
        "translation_ordered_orbit_count": len(translation_orbits),
        "affine_ordered_orbit_count": len(affine_orbits),
        "translation_orbit_sizes": sorted({len(orb) for orb in translation_orbits}),
        "affine_orbit_sizes": sorted({len(orb) for orb in affine_orbits}),
        "translation_fixed_ordered_triples": len(translation_fixed),
        "affine_fixed_ordered_triples": len(affine_fixed),
        "translation_singleton_orbits": len(translation_singletons),
        "affine_singleton_orbits": len(affine_singletons),
        "translation_orbit_representatives": [sorted(orb)[0] for orb in translation_orbits[:8]],
        "affine_orbit_representatives": [sorted(orb)[0] for orb in affine_orbits[:8]],
        "finite_theorem": "A source-free localizer compatible with translation or affine source-gauge symmetry would have to select a fixed ordered triple, equivalently a singleton orbit.  No ordered distinct Z12 triple is fixed by all translations, and no affine singleton orbit exists.  Therefore P2740's pointwise chirality cannot be converted into a non-premise strict signed source by an orbit-safe ordered-triple localizer on current data.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_localizer_obligation": scan["all_patterns_have_hits"],
        "candidate_domain_is_nonempty": audit["ordered_distinct_triples"] == 1320,
        "translation_safe_fixed_ordered_triple_exists": audit["translation_fixed_ordered_triples"] > 0,
        "affine_safe_singleton_ordered_orbit_exists": audit["affine_singleton_orbits"] > 0,
        "strict_source_localizer_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_source_localizer": all(facts.values()),
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "No translation/affine-safe fixed ordered source triple exists, so the P2740 chirality still requires an external label/order premise and cannot fix lambda/P2721.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["localizer_audit"]
    lines = [
        "# P2741/S1691 source-triple localizer fixed-point no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite fixed-point audit",
        f"- ordered_distinct_triples={audit['ordered_distinct_triples']}",
        f"- translation_group_size={audit['translation_group_size']}",
        f"- affine_group_size={audit['affine_group_size']}",
        f"- translation_ordered_orbit_count={audit['translation_ordered_orbit_count']}",
        f"- affine_ordered_orbit_count={audit['affine_ordered_orbit_count']}",
        f"- translation_orbit_sizes={audit['translation_orbit_sizes']}",
        f"- affine_orbit_sizes={audit['affine_orbit_sizes']}",
        f"- translation_fixed_ordered_triples={audit['translation_fixed_ordered_triples']}",
        f"- affine_fixed_ordered_triples={audit['affine_fixed_ordered_triples']}",
        f"- translation_singleton_orbits={audit['translation_singleton_orbits']}",
        f"- affine_singleton_orbits={audit['affine_singleton_orbits']}",
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
    audit = localizer_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2741_SOURCE_TRIPLE_LOCALIZER_FIXED_POINT_NO_GO" if not acceptance["accepted_as_strict_source_localizer"] else "P2741_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "translation/affine-safe strict source-localizer for the exact P2740 ordered-triple chirality object",
        "content_evidence_scan": scan,
        "localizer_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue P2740 by merely choosing an ordered source triple: P2741 proves that no translation/affine-safe fixed ordered triple or singleton orbit exists, so such a choice is still a label/order premise.  The next proof-grade move must either supply an independent strict source-localizer theorem that breaks this source-gauge obstruction and a P2721 coupling theorem, or pivot to a different strict signed observable with a nonzero orbit-safe signed value; otherwise preserve the P2697-P2741 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2741/S1691 source-triple localizer fixed-point no-go", "## P2741/S1691 source-triple localizer fixed-point no-go\n\n`P2741/S1691` audits the exact missing premise left by `P2740/S1690`: a strict localizer selecting one ordered distinct `Z12` source triple.  The finite fixed-point/orbit computation finds `1320` ordered triples, `110` translation ordered orbits of size `12`, `34` affine ordered orbits, no translation-fixed ordered triple, and no affine singleton orbit.  Thus an ordered-triple choice remains a source-label/order premise unless a new source-localizer theorem is exported.  No `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2741/S1691 source-triple localizer Ltotal guard", "## P2741/S1691 source-triple localizer Ltotal guard\n\n`P2741/S1691` adds no variational source term: the finite fixed-point audit finds no translation/affine-safe ordered-triple localizer for the `P2740` chirality.  Since the required source-localizer and `P2721` coupling theorem remain unexported, this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current source-triple localizer fixed-point no-go guardrail (P2741/S1691, 2026-06-14)", "## Current source-triple localizer fixed-point no-go guardrail (P2741/S1691, 2026-06-14)\n\n- P2741 audits the exact missing premise left by P2740: a strict source-localizer choosing one ordered distinct `Z12` source triple for cyclic-order chirality.\n- The finite fixed-point/orbit computation finds `1320` ordered triples, `110` translation ordered orbits of size `12`, `34` affine ordered orbits, zero translation-fixed ordered triples, and zero affine singleton orbits; an ordered-triple choice remains a label/order premise on current data.\n- Do not promote source-triple chirality or an externally chosen ordered triple to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must supply an independent strict source-localizer theorem plus `P2721` coupling theorem, or pivot to a different strict signed observable with nonzero orbit-safe signed value; otherwise preserve the P2697-P2741 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
