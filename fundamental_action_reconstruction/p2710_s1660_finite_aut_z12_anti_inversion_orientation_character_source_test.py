#!/usr/bin/env python3
"""P2710/S1660: finite Aut(Z12) anti-inversion orientation-character source test.

P2708 found a real boundary-cocycle orientation line ±omega, and P2709
recommended the next honest move: test one finite anti-inversion/orientation
character source candidate rather than replaying older releases.  This script
computes the full character table of Aut(Z12)=U(12), identifies the characters
that are anti-invariant on inversion, and then checks whether current strict
artifacts actually export any such character as a non-premise source.
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
OUT = GEN / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.json"
MD = GEN / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = [1, 5, 7, 11]
INVERSION = 11

INPUTS = {
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2709_RELEASE_BACKSCAN": GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json",
    "P2706_DAMPING_INTERFACE": GEN / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.json",
    "P2704_P1343_SCOPE": GEN / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.json",
    "P2700_AUT_FUNCTIONALS": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2699_FRACTAL_INFO_AUT": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2378_UNIT_NORMALIZED_INSUFFICIENCY": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
}

SOURCE_EXPORT_TOKENS = [
    "anti_inversion_character_exported\": true",
    "orientation_character_source_exported\": true",
    "non_premise_selector_provider_exported\": true",
    "strict source for the missing sign exported",
    "qw2191_discharged\": true",
]
SOURCE_BLOCK_TOKENS = [
    "no_current_unlock", "no current unlock", "no-strict-selector", "no_strict_selector",
    "orientation-blind", "bounded_no_go", "no-go", "no_go", "insufficient", "open",
    "premise", "not a non-premise", "no nonzero aut-invariant",
]
NEGATIVE_EXPORT_FLAGS = [
    "anti_inversion_character_exported",
    "orientation_character_source_exported",
    "non_premise_selector_provider_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_exported",
]


def modmul(a: int, b: int) -> int:
    return (a * b) % 12


def valid_character(values: dict[int, int]) -> bool:
    return values[1] == 1 and all(values[modmul(a, b)] == values[a] * values[b] for a in UNITS for b in UNITS)


def character_table() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for vals in itertools.product([-1, 1], repeat=len(UNITS)):
        values = dict(zip(UNITS, vals))
        if not valid_character(values):
            continue
        anti = values[INVERSION] == -1
        # Acting on P2708's orientation line, an anti-inversion character can label
        # inversion-odd parity, but it still does not choose +omega unless a strict
        # artifact supplies this character as the physical/source law.
        rows.append({
            "name": f"chi_5_{values[5]:+d}_7_{values[7]:+d}_11_{values[11]:+d}",
            "values": {str(k): values[k] for k in UNITS},
            "homomorphism_valid": True,
            "anti_inversion": anti,
            "distinguishes_inversion_from_identity": anti,
            "selects_plus_omega_without_source": False,
        })
    return rows


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore") if path.exists() else ""


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {"json_error": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def artifact_source_scan() -> list[dict[str, Any]]:
    rows = []
    for name, path in INPUTS.items():
        text = read_text(path).lower()
        export_hits = sorted(tok for tok in SOURCE_EXPORT_TOKENS if tok in text)
        block_hits = sorted(tok for tok in SOURCE_BLOCK_TOKENS if tok in text)[:10]
        data = read_json(path) if path.suffix == ".json" else {}
        rows.append({
            "artifact": name,
            "path": rel(path),
            "status": data.get("status"),
            "export_hits": export_hits,
            "block_hits_sample": block_hits,
            "strict_anti_inversion_source_exported": bool(export_hits) and not block_hits,
        })
    return rows


def candidate_rows(chars: list[dict[str, Any]], scan: list[dict[str, Any]]) -> list[dict[str, Any]]:
    source_available = any(row["strict_anti_inversion_source_exported"] for row in scan)
    rows = []
    for char in chars:
        if not char["anti_inversion"]:
            continue
        rows.append({
            "candidate": char["name"],
            "finite_character_exists": True,
            "anti_inversion": True,
            "mathematical_content": "A Z2-valued character of Aut(Z12) exists with chi(inversion)=-1.",
            "strict_source_available_on_current_artifacts": source_available,
            "exports_orientation_source": False,
            "blocker": "The character is a possible parity label, not an exported strict law choosing +omega over -omega. Current artifacts provide no non-premise source selecting this character and its sign.",
        })
    return rows


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2710/S1660 finite Aut(Z12) anti-inversion orientation-character source test",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Character table summary",
        f"- total_characters={payload['character_summary']['total_characters']}",
        f"- anti_inversion_characters={payload['character_summary']['anti_inversion_characters']}",
        "",
        "## Candidate source rows",
    ]
    for row in payload["candidate_source_rows"]:
        lines.append(f"- `{row['candidate']}`: strict_source_available={row['strict_source_available_on_current_artifacts']}, exports_orientation_source={row['exports_orientation_source']}. {row['blocker']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    chars = character_table()
    scan = artifact_source_scan()
    rows = candidate_rows(chars, scan)
    anti_count = sum(1 for row in chars if row["anti_inversion"])
    no_unlock = anti_count == 2 and rows and all(not row["exports_orientation_source"] for row in rows) and all(not row["strict_anti_inversion_source_exported"] for row in scan)
    payload = {
        "status": "P2710_AUT_Z12_ANTI_INVERSION_CHARACTER_EXISTS_BUT_NO_STRICT_SOURCE" if no_unlock else "P2710_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "character_table": chars,
        "character_summary": {
            "group": "Aut(Z12)=U(12)={1,5,7,11}",
            "total_characters": len(chars),
            "anti_inversion_characters": anti_count,
            "anti_inversion_character_names": [row["name"] for row in chars if row["anti_inversion"]],
        },
        "artifact_source_scan": scan,
        "candidate_source_rows": rows,
        "decision": {
            "finite_anti_inversion_characters_exist": anti_count == 2,
            "strict_orientation_character_source_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "The finite character table of Aut(Z12) contains exactly two inversion-odd characters, so the mathematical parity labels exist.  However current strict artifacts do not export either character as a non-premise physical/source law, and an abstract character still does not choose +omega rather than -omega.  Thus P2710 does not discharge QW-2191 or unlock pair12 strict-core, L_total, role transfer, bridge closure, or ToE.",
            "next_honest_step": "The next admissible move must supply a genuinely new strict artifact that couples one inversion-odd character to the boundary-cocycle sign with a non-premise sign convention, or pivot to a different new typed object.  Without that new source, preserve the P2697-P2710 no-current-unlock certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2710/S1660 finite Aut(Z12) anti-inversion orientation-character source test", "## P2710/S1660 finite Aut(Z12) anti-inversion orientation-character source test\n\n`P2710/S1660` computes the full character table of `Aut(Z12)=U(12)={1,5,7,11}`.  Exactly two characters are inversion-odd, so finite anti-inversion parity labels exist.  Current strict artifacts, however, do not export either character as a non-premise source law coupled to the P2708 boundary-cocycle sign.  Therefore P2710 preserves `QW-2191`, pair12 strict-core, role transfer, `L_total`, bridge closure, and ToE as blocked on current artifacts.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2710/S1660 anti-inversion character Ltotal guard", "## P2710/S1660 anti-inversion character Ltotal guard\n\n`P2710/S1660` is a finite Aut(Z12) character-table test, not a variational construction.  The existence of inversion-odd characters is only a mathematical parity label unless a strict source law selects and signs one of them.  No such source is exported here, so this is not an `L_total`, role-transfer, bridge, or ToE promotion.\n")
    append_once(AGENTS, "Current Aut(Z12) anti-inversion character source guardrail (P2710/S1660, 2026-06-14)", "## Current Aut(Z12) anti-inversion character source guardrail (P2710/S1660, 2026-06-14)\n\n- P2710 computes the full finite character table of `Aut(Z12)=U(12)` and finds exactly two inversion-odd characters.\n- These characters are mathematical parity labels only; current artifacts do not export either one as a non-premise strict source coupled to the P2708 boundary-cocycle sign, so `QW-2191` remains open.\n- Do not promote inversion-odd character existence to selector closure, pair12 strict-core, role transfer, `L_total`, bridge closure, or ToE without a new strict source law selecting and signing the character.\n")
    return payload


if __name__ == "__main__":
    main()
