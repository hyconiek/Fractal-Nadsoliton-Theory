#!/usr/bin/env python3
"""P2748/S1698: absence-of-selector self-synchronization no-go.

This audit addresses the proposed idea that the nadsoliton, as pure
Information, may contain information about the absence of a selector and thereby
synchronize with a selector.  The ontology is respected: there is no lower
informational layer.  The finite question is whether an Aut(Z12)-invariant
"absence/no-selector" datum can itself map to a directed selector torsor.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2747 = GEN / "p2747_s1697_z12_cubic_phase_orbit_signed_observable_audit.json"
OUT = GEN / "p2748_s1698_absence_of_selector_self_synchronization_no_go.json"
MD = GEN / "p2748_s1698_absence_of_selector_self_synchronization_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = (1, 5, 7, 11)
TORSOR = (-1, 1)
ABSENCE_SINGLETON = ("absence_of_selector",)
ABSENCE_BIT = ("selector_absent_false", "selector_absent_true")
CONTENT_PATTERNS = {
    "post_p2747_boundary": r"P2747|cubic phase|P2697-P2747|no-new-live-frontier",
    "absence_selector_language": r"absence of selector|no selector|selector absence|brak selektora|selector",
    "ontology_boundary": r"nadsoliton.*primordial|pure Information|czyst.*Informac|no lower information|informational layer",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "absence_information_promoted_to_selector",
    "self_synchronization_selector_exported",
    "nonpremise_orientation_sign_exported",
    "p2721_polarity_coupling_exported",
    "lambda_fixed",
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
    return {"content_pattern_count": len(rows), "rows": rows, "hit_counts": {r["content_lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def orientation_action(unit: int, value: int) -> int:
    # P2708/P2714 boundary: units 7 and 11 reverse the boundary-cocycle orientation.
    return -value if unit in (7, 11) else value


def trivial_action(_unit: int, state: str) -> str:
    return state


def odd_absence_action(unit: int, value: int) -> int:
    # This is not absence alone: it adds the missing inversion-odd sign representation.
    return orientation_action(unit, value)


def enumerate_maps(domain: tuple[Any, ...], codomain: tuple[int, ...]) -> list[dict[Any, int]]:
    return [dict(zip(domain, values)) for values in product(codomain, repeat=len(domain))]


def equivariant_maps_from_trivial_absence(domain: tuple[str, ...]) -> list[dict[str, int]]:
    good = []
    for mapping in enumerate_maps(domain, TORSOR):
        ok = True
        for unit in UNITS:
            for state in domain:
                if mapping[trivial_action(unit, state)] != orientation_action(unit, mapping[state]):
                    ok = False
        if ok:
            good.append(mapping)
    return good


def equivariant_maps_from_odd_signed_source() -> list[dict[int, int]]:
    good = []
    for mapping in enumerate_maps(TORSOR, TORSOR):
        ok = True
        for unit in UNITS:
            for state in TORSOR:
                if mapping[odd_absence_action(unit, state)] != orientation_action(unit, mapping[state]):
                    ok = False
        if ok:
            good.append(mapping)
    return good


def synchronization_audit() -> dict[str, Any]:
    singleton_maps = enumerate_maps(ABSENCE_SINGLETON, TORSOR)
    bit_maps = enumerate_maps(ABSENCE_BIT, TORSOR)
    singleton_equivariant = equivariant_maps_from_trivial_absence(ABSENCE_SINGLETON)
    bit_equivariant = equivariant_maps_from_trivial_absence(ABSENCE_BIT)
    odd_equivariant = equivariant_maps_from_odd_signed_source()
    relation_count_singleton = 2 ** (len(ABSENCE_SINGLETON) * len(TORSOR))
    relation_count_bit = 2 ** (len(ABSENCE_BIT) * len(TORSOR))
    return {
        "typed_candidate": "absence-of-selector information as a self-synchronizing selector source",
        "aut_units": list(UNITS),
        "orientation_torsor": list(TORSOR),
        "orientation_reversing_units": [7, 11],
        "singleton_absence_map_count": len(singleton_maps),
        "singleton_absence_equivariant_map_count": len(singleton_equivariant),
        "bit_absence_map_count": len(bit_maps),
        "bit_absence_equivariant_map_count": len(bit_equivariant),
        "singleton_relation_count": relation_count_singleton,
        "bit_relation_count": relation_count_bit,
        "odd_signed_source_equivariant_map_count": len(odd_equivariant),
        "odd_signed_source_equivariant_maps": odd_equivariant,
        "finite_theorem": "An Aut(Z12)-invariant absence/no-selector datum has trivial action.  Because units 7 and 11 reverse the orientation torsor, every map from a fixed absence state (or a trivial absence bit) to the selector torsor fails equivariance: the finite counts are 0 equivariant maps from the singleton and 0 from the trivial bit.  Equivariant synchronization becomes possible only after replacing absence by a new inversion-odd signed source; that has 2 equivariant maps, but it is exactly the missing new sign, not information about absence alone.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_current_boundaries": scan["all_patterns_have_hits"],
        "absence_candidate_formalized": True,
        "absence_singleton_exports_equivariant_selector": audit["singleton_absence_equivariant_map_count"] > 0,
        "absence_bit_exports_equivariant_selector": audit["bit_absence_equivariant_map_count"] > 0,
        "new_inversion_odd_signed_source_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_self_synchronizing_selector": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "Absence/no-selector information is Aut-invariant and has no equivariant map to the orientation selector torsor.  A synchronizing map requires adding a new inversion-odd signed source, which current artifacts do not export.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["synchronization_audit"]
    lines = [
        "# P2748/S1698 absence-of-selector self-synchronization no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite synchronization audit",
        f"- orientation_reversing_units={a['orientation_reversing_units']}",
        f"- singleton_absence_map_count={a['singleton_absence_map_count']}",
        f"- singleton_absence_equivariant_map_count={a['singleton_absence_equivariant_map_count']}",
        f"- bit_absence_map_count={a['bit_absence_map_count']}",
        f"- bit_absence_equivariant_map_count={a['bit_absence_equivariant_map_count']}",
        f"- odd_signed_source_equivariant_map_count={a['odd_signed_source_equivariant_map_count']}",
        "",
        "## Theorem statement",
        a["finite_theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2747 = read_json(P2747)
    scan = evidence_scan()
    audit = synchronization_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2748_ABSENCE_OF_SELECTOR_SELF_SYNCHRONIZATION_NO_GO",
        "input_hashes": {"P2747_CUBIC_PHASE_AUDIT": sha(P2747)},
        "input_statuses": {"P2747_CUBIC_PHASE_AUDIT": p2747.get("status")},
        "audited_candidate_class": "absence-of-selector information as self-synchronizing selector after P2747",
        "ontology_note": "The nadsoliton remains primordial pure Information; this audit does not introduce a lower informational layer.",
        "content_evidence_scan": scan,
        "synchronization_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote 'information about absence of selector' to a selector.  P2748 proves that an Aut-invariant absence/no-selector datum has zero equivariant maps to the orientation torsor; synchronization requires an added inversion-odd signed source, which is precisely the missing new object.  The next proof-grade move must either construct that concrete inversion-odd signed source with a P2721 coupling theorem, or pivot to a different typed object; otherwise preserve the P2697-P2748 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2748/S1698 absence-of-selector self-synchronization no-go", "## P2748/S1698 absence-of-selector self-synchronization no-go\n\n`P2748/S1698` audits the idea that the nadsoliton, as primordial pure Information, could contain information about the absence of a selector and thereby synchronize with a selector.  Formalized as an Aut-invariant absence/no-selector datum, the finite map count gives `0` equivariant maps from a singleton absence state and `0` from a trivial absence bit to the orientation torsor, because units `7` and `11` reverse the torsor.  Equivariant maps exist only after adding a new inversion-odd signed source, which is the missing object rather than absence alone.  No `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2748/S1698 absence selector Ltotal guard", "## P2748/S1698 absence selector Ltotal guard\n\n`P2748/S1698` adds no variational source term.  Information about selector absence, treated as an Aut-invariant datum inside the nadsoliton ontology, has no equivariant map to the orientation torsor; a synchronizing map requires an additional inversion-odd signed source and `P2721` coupling theorem.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current absence-of-selector self-synchronization no-go guardrail (P2748/S1698, 2026-06-14)", "## Current absence-of-selector self-synchronization no-go guardrail (P2748/S1698, 2026-06-14)\n\n- P2748 audits the idea that the nadsoliton, as primordial pure Information, could contain information about the absence of a selector and thereby synchronize with a selector, without introducing any lower informational layer.\n- The finite Aut(Z12)-equivariance test finds `0` equivariant maps from a singleton absence state and `0` from a trivial absence bit to the orientation torsor; units `7` and `11` reverse the torsor, so absence/no-selector information alone is orientation-blind.\n- Do not promote absence-of-selector information to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must construct a concrete inversion-odd signed source with explicit `P2721` coupling, pivot to a different typed object, or preserve the P2697-P2748 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
