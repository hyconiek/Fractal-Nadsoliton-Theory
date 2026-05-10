#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P789 = GENERATED / "p789_current_strict_alpha_s_normalized_boundary_interface_candidate_probe.json"
IN_F704 = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_F789 = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"

OUT = GENERATED / "p790_current_strict_alpha_s_canonical_anchor_upgrade_probe.json"
OUT_SUMMARY = GENERATED / "p790_current_strict_alpha_s_canonical_anchor_upgrade_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return (
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("-", " ")
        .replace("_", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def find_anchor_artifacts() -> list[str]:
    patterns = [
        "*alpha_s*anchor*",
        "*alpha_s*boundary*anchor*",
        "*canonical*anchor*alpha_s*",
    ]
    hits: set[str] = set()
    for pattern in patterns:
        for path in GENERATED.glob(pattern):
            hits.add(str(path.relative_to(REPO)))
        for path in ROOT.glob(pattern):
            if path.suffix not in {".md", ".json"}:
                continue
            hits.add(str(path.relative_to(REPO)))
    return sorted(
        hit
        for hit in hits
        if "p790_current_strict_alpha_s_canonical_anchor_upgrade_probe" not in hit
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P789, IN_F704, IN_F789, IN_N703]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P790",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p789 = load_json(IN_P789)
    f704 = load_json(IN_F704)
    f789 = load_json(IN_F789)
    n703_text = IN_N703.read_text(encoding="utf-8")

    strongest_family_id = str(p789.get("strongest_family_id"))
    strongest_family = ((p789.get("candidate_families") or {}).get(strongest_family_id)) or {}
    anchor_artifacts = find_anchor_artifacts()

    more_than_one_family_present = len(p789.get("candidate_families") or {}) > 1
    selection_principle_exported = bool(anchor_artifacts)
    semantic_upgrade_exported = contains_all(
        json.dumps(f789, ensure_ascii=True),
        ["canonical alpha_s boundary anchor", "already exported"],
    )
    n703_host_matching_block_present = contains_all(
        n703_text,
        [
            "they are not yet physical masses in gev",
            "they do not constitute a standard model match claim",
        ],
    )
    nf_attachment_exported = bool(strongest_family.get("checks", {}).get("nf_mapping_exported"))

    checks = [
        {
            "id": "strongest_family_present_and_dimensionless",
            "pass": (
                p789.get("candidate_family_constructible") is True
                and strongest_family.get("checks", {}).get("dimensionless_points_present") is True
                and strongest_family.get("checks", {}).get("no_gev_fields_used") is True
            ),
            "details": "The strongest candidate family must already exist and stay fully dimensionless.",
        },
        {
            "id": "multiple_candidate_families_prevent_silent_selection",
            "pass": more_than_one_family_present,
            "details": "There is more than one candidate family, so the repo cannot silently treat one family as canonical without an exported selection principle.",
        },
        {
            "id": "no_exported_anchor_selection_artifact_detected",
            "pass": not selection_principle_exported,
            "details": "No current repo artifact explicitly selects the strongest family as the canonical alpha_s anchor.",
        },
        {
            "id": "n703_meaning_boundary_blocks_semantic_overread",
            "pass": n703_host_matching_block_present,
            "details": "N703 keeps the mass-proxy layer in internal dimensionless meaning scope, blocking silent semantic upgrade into an alpha_s boundary role.",
        },
        {
            "id": "no_nf_attachment_exported_for_strongest_family",
            "pass": not nf_attachment_exported,
            "details": "No exported n_f attachment exists for the strongest candidate family.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    if not failed_checks:
        status = "P790_CURRENT_STRICT_ALPHA_S_CANONICAL_ANCHOR_UPGRADE_BLOCKED_ON_CURRENT_REPO_STATE"
    else:
        status = "P790_REQUIRES_REVIEW"

    narrowed_anchor_blockers = [
        "No exported selection principle chooses f704_max_mode_anchor_family over the other candidate normalized families.",
        "No exported semantic-upgrade rule attaches the chosen F704 anchor candidate to the alpha_s boundary role rather than to generic internal quadratic meaning only.",
        "No exported n_f attachment maps the candidate anchor into the alpha_s boundary-sector interface.",
    ]

    artifact = {
        "stage": "P790",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p789_candidate_probe": rel(IN_P789),
            "f704_mass_observable": rel(IN_F704),
            "f789_interface_target": rel(IN_F789),
            "n703_meaning_theorem": rel(IN_N703),
        },
        "strongest_family_id": strongest_family_id,
        "strongest_family_snapshot": {
            "normalization_rule": strongest_family.get("normalization_rule"),
            "mu0_tilde_candidate": strongest_family.get("mu0_tilde_candidate"),
            "normalized_validation_points_candidate": strongest_family.get("normalized_validation_points_candidate"),
            "strict_input_chain": strongest_family.get("strict_input_chain"),
        },
        "anchor_artifacts_detected": anchor_artifacts,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "P789 already shows that f704_max_mode_anchor_family is the strongest current dimensionless candidate family.",
            "P790 shows that upgrading it into a canonical alpha_s anchor is still blocked, not because the family is numerically absent, but because the repo exports no selection principle, no semantic-upgrade rule, and no n_f attachment for that family.",
            "This narrows the anchor blocker further: the issue is now canonical selection semantics, not construction of dimensionless candidate points.",
        ],
        "narrowed_anchor_blockers": narrowed_anchor_blockers,
        "recommended_next_packet": {
            "id": "F790_CURRENT_STRICT_ALPHA_S_CANONICAL_ANCHOR_SELECTION_TARGET_PACKET",
            "goal": "Freeze the exact missing canonical-anchor selection object needed before any normalized alpha_s boundary family can be promoted.",
            "minimum_fields": [
                "candidate_anchor_family_id",
                "selection_principle_ref",
                "anchor_to_boundary_role_rule_ref",
                "n_f_attachment_rule_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P790",
        "status": status,
        "as_of": AS_OF,
        "strongest_family_id": strongest_family_id,
        "narrowed_anchor_blocker_count": len(narrowed_anchor_blockers),
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
