#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F784 = GENERATED / "f784_current_minimal_strict_sm_gr_bridge_registry_packet.json"
IN_RELEASE7 = REPO / "RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_MODEL_BRIEF.md"
IN_WORKING_NOTE = REPO / "WORKING_NOTE_LEGACY_KEEP_CUT_AND_MINIMAL_STRICT_SM_GR_BRIDGE.md"

OUT = GENERATED / "p785_current_minimal_strict_sm_gr_bridge_registry_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p785_current_minimal_strict_sm_gr_bridge_registry_boundary_audit_probe_summary.json"


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
        .replace("‑", "-")
        .replace("–", "-")
        .replace("—", "-")
        .replace("→", "->")
    )


def contains_any(texts: list[str], needles: list[str]) -> bool:
    hay = normalize("\n".join(texts))
    return any(normalize(needle) in hay for needle in needles)


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F784, IN_RELEASE7, IN_WORKING_NOTE]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P785",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f784 = load_json(IN_F784)
    release7_text = IN_RELEASE7.read_text(encoding="utf-8")
    working_note_text = IN_WORKING_NOTE.read_text(encoding="utf-8")

    targets = f784.get("targets") or {}
    mass = targets.get("mass_ratio_ordering_layer") or {}
    sin2 = targets.get("sin2_theta_w_eff") or {}
    alpha_s = targets.get("alpha_s_boundary_mu0_alpha0") or {}
    ggrav = targets.get("g_dimensionless_mu_ref") or {}

    checks_spec = [
        {
            "id": "packet_hard_limits_no_sm_identification",
            "pass": contains_any(f784.get("hard_limits") or [], ["Standard Model identification"]),
            "details": "F784 hard limits must explicitly deny Standard Model identification.",
        },
        {
            "id": "packet_hard_limits_no_proxy_to_gev_strict_claim",
            "pass": contains_any(f784.get("hard_limits") or [], ["proxy-to-GeV"]),
            "details": "F784 hard limits must explicitly deny strict proxy-to-GeV calibration.",
        },
        {
            "id": "release7_hard_limits_match_boundary",
            "pass": contains_all(
                release7_text,
                [
                    'No Standard Model identification / "match" claim',
                    "No physical-unit proxy->GeV calibration map",
                ],
            ),
            "details": "Release 7 brief must explicitly keep SM identification and proxy->GeV outside strict scope.",
        },
        {
            "id": "mass_layer_avoids_nonstrict_host_matching",
            "pass": (
                mass.get("status_label") == "strict-derived"
                and contains_any(mass.get("hard_limits") or [], ["host matching", "physical mass units"])
                and not contains_any(mass.get("source_artifacts") or [], ["p710_", "sm_mass_targets", "proxy_to_gev"])
            ),
            "details": "Mass layer must remain internal proxy/order only, with no hidden calibration or host-matching artifact.",
        },
        {
            "id": "sin2_entry_keeps_legacy_semantic_blocker",
            "pass": (
                sin2.get("status_label") == "strict-derived"
                and contains_any(sin2.get("promotion_blockers") or [], ["legacy Weinberg", "semantic-transfer"])
                and contains_any(sin2.get("hard_limits") or [], ["legacy Weinberg role"])
            ),
            "details": "sin2 entry must keep the legacy semantic-transfer blocker explicit.",
        },
        {
            "id": "alpha_s_entry_keeps_frozen_ansatz_blocker",
            "pass": (
                alpha_s.get("status_label") == "strict-derived"
                and contains_any(alpha_s.get("promotion_blockers") or [], ["frozen ansatz", "kernel-invariants-only"])
                and contains_any(alpha_s.get("hard_limits") or [], ["full QCD closure"])
            ),
            "details": "alpha_s boundary entry must keep the frozen-ansatz / non-uniqueness blocker explicit.",
        },
        {
            "id": "gravity_entry_marks_external_origin_and_open_internal_origin",
            "pass": (
                ggrav.get("status_label") == "strict-derived-with-external-observable-origin"
                and ggrav.get("external_observable_origin") is True
                and contains_any(ggrav.get("promotion_blockers") or [], ["external_dimensionless_observable", "internal origin"])
            ),
            "details": "Gravity entry must keep external origin explicit and internal origin open.",
        },
        {
            "id": "working_note_cut_rules_respected",
            "pass": contains_all(
                working_note_text,
                [
                    "Cut any silent transfer of legacy physical roles onto `K_strict_gate`.",
                    "proxy -> GeV",
                    "Do not rescan parameters after declaring the bridge target set.",
                ],
            ),
            "details": "Boundary audit expects the working note cut rules to remain present and active.",
        },
        {
            "id": "packet_does_not_claim_toe_closure",
            "pass": contains_any(f784.get("hard_limits") or [], ["ToE closure"]),
            "details": "F784 must explicitly deny ToE closure.",
        },
    ]

    checks = []
    failed = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "pass": bool(item["pass"]),
                "details": item["details"],
            }
        )
        if not item["pass"]:
            failed.append(item["id"])

    if failed:
        status = "P785_CURRENT_MINIMAL_BRIDGE_REGISTRY_BOUNDARY_AUDIT_REQUIRES_REVIEW"
    else:
        status = "P785_CURRENT_MINIMAL_BRIDGE_REGISTRY_BOUNDARY_AUDIT_PASS_WITH_BOUNDARIES_EXPLICIT"

    artifact = {
        "stage": "P785",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f784_packet": rel(IN_F784),
            "release7_brief": rel(IN_RELEASE7),
            "working_note": rel(IN_WORKING_NOTE),
        },
        "checks": checks,
        "failed_checks": failed,
        "current_honest_reading": [
            "F784 does not currently leak proxy-to-GeV or host-matching artifacts into the strict mass bridge layer.",
            "F784 keeps the Weinberg semantic-transfer blocker explicit on the sin2 entry.",
            "F784 keeps the alpha_s frozen-ansatz blocker explicit on the boundary entry.",
            "F784 keeps the gravity external-origin boundary explicit and does not relabel it as internal origin of G.",
        ],
        "recommended_next_move": "Attack the QW-2093 formula layer directly: test whether sin2 and alpha_s can be re-expressed through a kernel-invariants-only map without the current frozen-ansatz / legacy-touchpoint residues.",
        "no_false_pass": True,
    }

    summary = {
        "stage": "P785",
        "status": status,
        "as_of": AS_OF,
        "n_checks": len(checks),
        "n_failed": len(failed),
        "failed_checks": failed,
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
