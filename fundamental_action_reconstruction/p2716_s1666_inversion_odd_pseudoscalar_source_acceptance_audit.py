#!/usr/bin/env python3
"""P2716/S1666: inversion-odd pseudoscalar source acceptance audit.

P2715 ruled out Aut-trivial scalar sources.  The remaining physically honest
possibility is not another scalar replay, but a genuine inversion-odd
pseudoscalar/chiral source: a signed strict datum whose Aut(Z12) action matches
the orientation torsor and whose nonzero sign is exported by the theory rather
than chosen by convention.  P2716 separates the representation-theoretic
possibility from the missing source-value problem.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.json"
MD = GEN / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = [1, 5, 7, 11]
SIGNS = [-1, 1]

INPUTS = {
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
    "P2715_SCALAR_SOURCE_NO_GO": GEN / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.json",
    "P2711_SIGN_COUPLING": GEN / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.json",
    "P2709_RELEASE_BACKSCAN": GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json",
    "P2706_DAMPING_INTERFACE": GEN / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.json",
}

SOURCE_EXPORT_TOKENS = [
    "strict_pseudoscalar_or_chiral_source_exported\": true",
    "inversion_odd_source_value_exported\": true",
    "nonzero_pseudoscalar_vev_exported\": true",
    "strict chiral source exported",
    "strict pseudoscalar source exported",
]

SOURCE_BLOCK_TOKENS = [
    "no-current-unlock", "no_current_unlock", "no strict source", "no non-premise",
    "premise", "convention", "orientation-blind", "blocked", "no-go", "no_go",
    "does not export", "missing", "insufficient",
]

NEGATIVE_EXPORT_FLAGS = [
    "strict_pseudoscalar_or_chiral_source_exported",
    "inversion_odd_source_value_exported",
    "nonzero_pseudoscalar_vev_exported",
    "pseudoscalar_source_fixes_orientation_torsor",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore") if path.exists() else ""


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sign_action(unit: int, value: int) -> int:
    if unit in (1, 5):
        return value
    if unit in (7, 11):
        return -value
    raise ValueError(f"not a Z12 unit: {unit}")


def equivariant_sign_maps() -> list[dict[str, Any]]:
    """Enumerate maps from a pseudoscalar sign torsor to the orientation torsor."""
    rows: list[dict[str, Any]] = []
    for images in itertools.product(SIGNS, repeat=len(SIGNS)):
        mapping = dict(zip(SIGNS, images))
        failures = []
        for unit in UNITS:
            for source_sign in SIGNS:
                lhs = mapping[sign_action(unit, source_sign)]
                rhs = sign_action(unit, mapping[source_sign])
                if lhs != rhs:
                    failures.append({
                        "unit": unit,
                        "source_sign": "+chi" if source_sign == 1 else "-chi",
                        "lhs": "+omega" if lhs == 1 else "-omega",
                        "rhs": "+omega" if rhs == 1 else "-omega",
                    })
        rows.append({
            "map": {"+chi" if key == 1 else "-chi": "+omega" if value == 1 else "-omega" for key, value in mapping.items()},
            "aut_equivariant": not failures,
            "orientation_fixing_if_source_sign_exported": not failures,
            "still_requires_nonzero_source_sign": not failures,
            "failure_count": len(failures),
            "failure_sample": failures[:4],
        })
    return rows


def artifact_pseudoscalar_source_scan() -> list[dict[str, Any]]:
    rows = []
    for key, path in INPUTS.items():
        text = read_text(path).lower()
        export_hits = sorted(token for token in SOURCE_EXPORT_TOKENS if token in text)
        block_hits = sorted(token for token in SOURCE_BLOCK_TOKENS if token in text)[:10]
        data = read_json(path) if path.suffix == ".json" and path.exists() else {}
        rows.append({
            "artifact": key,
            "path": rel(path),
            "status": data.get("status"),
            "pseudoscalar_export_hits": export_hits,
            "block_hits_sample": block_hits,
            "exports_strict_nonzero_pseudoscalar_source": bool(export_hits) and not block_hits,
        })
    return rows


def acceptance_rows(maps: list[dict[str, Any]], scan: list[dict[str, Any]]) -> list[dict[str, Any]]:
    source_exported = any(row["exports_strict_nonzero_pseudoscalar_source"] for row in scan)
    rows = []
    for row in maps:
        if not row["aut_equivariant"]:
            continue
        rows.append({
            "equivariant_map": row["map"],
            "representation_theoretically_admissible": True,
            "strict_nonzero_source_sign_exported": source_exported,
            "fixes_orientation_torsor_now": source_exported,
            "blocker": "The map would transfer a strict signed pseudoscalar value to the orientation torsor, but current artifacts export no nonzero inversion-odd source value.",
        })
    return rows


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2716/S1666 inversion-odd pseudoscalar source acceptance audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Representation-theoretic maps",
        f"- candidate_maps={payload['finite_summary']['candidate_maps']}",
        f"- equivariant_maps={payload['finite_summary']['equivariant_maps']}",
        "",
        "## Acceptance rows",
    ]
    for row in payload["acceptance_rows"]:
        lines.append(f"- map={row['equivariant_map']}: source_exported={row['strict_nonzero_source_sign_exported']}, fixes_now={row['fixes_orientation_torsor_now']}. {row['blocker']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    maps = equivariant_sign_maps()
    scan = artifact_pseudoscalar_source_scan()
    rows = acceptance_rows(maps, scan)
    equivariant_count = sum(1 for row in maps if row["aut_equivariant"])
    source_exported = any(row["exports_strict_nonzero_pseudoscalar_source"] for row in scan)
    no_unlock = equivariant_count == 2 and not source_exported and all(not row["fixes_orientation_torsor_now"] for row in rows)
    payload = {
        "status": "P2716_PSEUDOSCALAR_REPRESENTATION_ADMISSIBLE_BUT_NO_STRICT_SOURCE_VALUE" if no_unlock else "P2716_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "finite_summary": {
            "source_domain": "inversion-odd pseudoscalar sign torsor {+chi,-chi}",
            "target_torsor": "P2708/P2714 orientation torsor {+omega,-omega}",
            "candidate_maps": len(maps),
            "equivariant_maps": equivariant_count,
            "strict_nonzero_pseudoscalar_source_exported": source_exported,
        },
        "equivariant_map_rows": maps,
        "artifact_pseudoscalar_source_scan": scan,
        "acceptance_rows": rows,
        "decision": {
            "pseudoscalar_representation_can_couple_to_orientation_torsor": equivariant_count == 2,
            "strict_pseudoscalar_or_chiral_source_exported": False,
            "inversion_odd_source_value_exported": False,
            "nonzero_pseudoscalar_vev_exported": False,
            "pseudoscalar_source_fixes_orientation_torsor": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2716 separates two issues.  Representation-theoretically, an inversion-odd pseudoscalar sign torsor can couple equivariantly to the +omega/-omega orientation torsor: two equivariant maps exist.  But current artifacts export no non-premise, nonzero signed pseudoscalar/chiral source value, so the maps have no strict input sign to transport and do not fix lambda or discharge QW-2191.",
            "next_honest_step": "A next admissible move must either construct/export a concrete strict inversion-odd pseudoscalar or chiral source value, then test one of the two equivariant maps as a bounded witness, or pivot to a different new typed object outside the closed lanes.  Without a signed source value, preserve the P2697-P2716 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2716/S1666 inversion-odd pseudoscalar source acceptance audit", "## P2716/S1666 inversion-odd pseudoscalar source acceptance audit\n\n`P2716/S1666` separates representation-theoretic admissibility from source export.  An inversion-odd pseudoscalar sign torsor `{+chi,-chi}` admits exactly two `Aut(Z12)`-equivariant maps to the orientation torsor `{+omega,-omega}`, so this is the correct physical kind of source.  However, current artifacts export no non-premise, nonzero signed pseudoscalar/chiral source value to feed either map.  Therefore no strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2716/S1666 pseudoscalar source Ltotal guard", "## P2716/S1666 pseudoscalar source Ltotal guard\n\n`P2716/S1666` identifies the admissible representation type for a future source: an inversion-odd pseudoscalar/chiral signed datum.  It does not itself export such a nonzero source value, so it is not a variational source construction and does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current inversion-odd pseudoscalar source acceptance guardrail (P2716/S1666, 2026-06-14)", "## Current inversion-odd pseudoscalar source acceptance guardrail (P2716/S1666, 2026-06-14)\n\n- P2716 finds that an inversion-odd pseudoscalar/chiral sign torsor is representation-theoretically the right kind of source: exactly two `Aut(Z12)`-equivariant maps to the `+omega/-omega` orientation torsor exist.\n- Current artifacts do not export a non-premise, nonzero signed pseudoscalar/chiral source value, so the admissible maps have no strict input sign and do not fix `lambda` or discharge `QW-2191`.\n- Do not promote pseudoscalar representation admissibility alone to selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE; a next admissible move must construct/export one concrete strict signed pseudoscalar/chiral source value or pivot to a different genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    main()
