#!/usr/bin/env python3
"""P2749/S1699: minimal inversion-odd source coupling-polarity audit.

P2748 showed that absence/no-selector information cannot synchronize to the
orientation torsor unless a new inversion-odd signed source is added.  This file
tests that exact next premise: even if such a minimal odd source sign is supplied,
does equivariance alone select the P2721 coupling polarity?
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
P2748 = GEN / "p2748_s1698_absence_of_selector_self_synchronization_no_go.json"
OUT = GEN / "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.json"
MD = GEN / "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = (1, 5, 7, 11)
SIGN = (-1, 1)
CONTENT_PATTERNS = {
    "post_p2748_missing_source": r"P2748|inversion-odd signed source|absence-of-selector|self-synchronization",
    "coupling_polarity_boundary": r"P2721 coupling|coupling polarity|polarity theorem|polarity selection",
    "selector_closure_boundary": r"lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
    "odd_source_prior_boundary": r"inversion-odd|orientation torsor|Aut\(Z12\)|units `?7`? and `?11`?",
}
NEGATIVE_EXPORT_FLAGS = [
    "minimal_odd_source_promoted_to_selector",
    "coupling_polarity_selected_by_equivariance",
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


def odd_action(unit: int, value: int) -> int:
    return -value if unit in (7, 11) else value


def enumerate_maps(domain: tuple[int, ...], codomain: tuple[int, ...]) -> list[dict[int, int]]:
    return [dict(zip(domain, values)) for values in product(codomain, repeat=len(domain))]


def equivariant_maps() -> list[dict[int, int]]:
    good = []
    for mapping in enumerate_maps(SIGN, SIGN):
        ok = True
        for unit in UNITS:
            for source_sign in SIGN:
                if mapping[odd_action(unit, source_sign)] != odd_action(unit, mapping[source_sign]):
                    ok = False
        if ok:
            good.append(mapping)
    return good


def compose_with_p2721_polarity(mapping: dict[int, int], polarity: int) -> dict[int, int]:
    return {key: polarity * value for key, value in mapping.items()}


def coupling_audit() -> dict[str, Any]:
    maps = equivariant_maps()
    normalized = [{str(k): v for k, v in sorted(m.items())} for m in maps]
    polarity_rows = []
    for index, mapping in enumerate(maps):
        for polarity in SIGN:
            coupled = compose_with_p2721_polarity(mapping, polarity)
            polarity_rows.append({"base_map_index": index, "p2721_polarity": polarity, "coupled_map": {str(k): v for k, v in sorted(coupled.items())}})
    unique_coupled = {json.dumps(row["coupled_map"], sort_keys=True) for row in polarity_rows}
    return {
        "typed_candidate": "minimal inversion-odd signed source with equivariant map to the orientation torsor and a P2721 polarity choice",
        "aut_units": list(UNITS),
        "orientation_reversing_units": [7, 11],
        "all_set_maps_count": len(enumerate_maps(SIGN, SIGN)),
        "equivariant_map_count": len(maps),
        "equivariant_maps": normalized,
        "p2721_polarity_rows": polarity_rows,
        "unique_coupled_map_count_after_p2721_polarity": len(unique_coupled),
        "polarity_pairing_witness": "Multiplying the coupling by P2721 polarity -1 exchanges the two equivariant maps; equivariance alone therefore leaves a twofold coupling-polarity ambiguity.",
        "finite_theorem": "A minimal inversion-odd signed source is the correct representation type requested by P2748: there are exactly 2 Aut(Z12)-equivariant maps from the odd source sign to the orientation torsor.  But those two maps are opposite coupling polarities, and composing with P2721 polarity exchanges them.  Therefore supplying an abstract odd sign is not enough: a further strict law must choose the source sign and the coupling polarity.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_p2748_boundary": scan["all_patterns_have_hits"],
        "minimal_odd_source_representation_admitted": audit["equivariant_map_count"] == 2,
        "coupling_polarity_unique": audit["equivariant_map_count"] == 1,
        "strict_source_sign_value_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_selector_source": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The minimal inversion-odd source has the right representation type, but equivariance leaves two opposite coupling polarities and current artifacts export neither a signed source value nor a P2721 polarity-selection theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["coupling_audit"]
    lines = [
        "# P2749/S1699 minimal inversion-odd source coupling-polarity audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite coupling audit",
        f"- orientation_reversing_units={a['orientation_reversing_units']}",
        f"- all_set_maps_count={a['all_set_maps_count']}",
        f"- equivariant_map_count={a['equivariant_map_count']}",
        f"- equivariant_maps={a['equivariant_maps']}",
        f"- unique_coupled_map_count_after_p2721_polarity={a['unique_coupled_map_count_after_p2721_polarity']}",
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
    p2748 = read_json(P2748)
    scan = evidence_scan()
    audit = coupling_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2749_MINIMAL_ODD_SOURCE_COUPLING_POLARITY_GAP",
        "input_hashes": {"P2748_ABSENCE_SELECTOR_NO_GO": sha(P2748)},
        "input_statuses": {"P2748_ABSENCE_SELECTOR_NO_GO": p2748.get("status")},
        "audited_candidate_class": "minimal inversion-odd signed source and P2721 coupling polarity after P2748",
        "content_evidence_scan": scan,
        "coupling_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote an abstract minimal inversion-odd source to selector closure.  P2749 shows the representation type is admissible but leaves exactly two opposite equivariant couplings, exchanged by P2721 polarity.  The next proof-grade move must provide a concrete strict source sign value and a theorem selecting one coupling polarity; otherwise pivot to a different typed object or preserve the P2697-P2749 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2749/S1699 minimal inversion-odd source coupling-polarity audit", "## P2749/S1699 minimal inversion-odd source coupling-polarity audit\n\n`P2749/S1699` tests the exact object requested after `P2748/S1698`: a minimal inversion-odd signed source and its coupling to the orientation torsor.  The finite Aut(Z12)-equivariance count gives exactly `2` maps from the odd source sign to the orientation torsor; they are opposite coupling polarities, and composing with `P2721` polarity exchanges them.  Thus the representation type is admissible, but current artifacts still export no concrete strict source sign value, no theorem selecting one coupling polarity, no `lambda/P2721` fixing, no `QW-2191` discharge, no selector closure, role transfer, `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2749/S1699 minimal odd source Ltotal guard", "## P2749/S1699 minimal odd source Ltotal guard\n\n`P2749/S1699` adds no variational source term.  A minimal inversion-odd source has the right equivariant representation type, but the two possible orientation-torsor couplings are opposite `P2721` polarities and no strict sign value/coupling-polarity theorem is exported.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current minimal inversion-odd source coupling-polarity guardrail (P2749/S1699, 2026-06-14)", "## Current minimal inversion-odd source coupling-polarity guardrail (P2749/S1699, 2026-06-14)\n\n- P2749 tests the exact object requested after P2748: a minimal inversion-odd signed source and its coupling to the orientation torsor.\n- The finite Aut(Z12)-equivariance count gives exactly `2` maps from the odd source sign to the orientation torsor; they are opposite coupling polarities, exchanged by `P2721` polarity.\n- Do not promote an abstract inversion-odd source representation to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a concrete strict source sign value and a theorem selecting one coupling polarity.  Otherwise pivot to a different typed object or preserve the P2697-P2749 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
