#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P790 = GENERATED / "p790_current_strict_alpha_s_canonical_anchor_upgrade_probe.json"
IN_F790 = GENERATED / "f790_current_strict_alpha_s_canonical_anchor_selection_target_packet.json"
IN_F447 = ROOT / "F447_CURRENT_STRICT_T169_QW2122_SCALAR_TO_T168_PER_SITE_VALUE_PROVIDER_POWERLAW_ELEMENT_ORDER_REFERENCE_PACKET.md"
IN_N483 = ROOT / "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM.md"
IN_T169 = ROOT / "T169_CURRENT_STRICT_CONSTRAINED_LIFT_FROM_QW2122_SCALAR_TO_T168_PER_SITE_PROVIDER_TARGET_SPEC.md"
IN_P415 = ROOT / "P415_CURRENT_STRICT_SIGMA_INT_TO_THETA_AX20_TYPED_Z12_PHASE_DENSITY_BERRY_SLOTLESS_ADMISSIBILITY_AUDIT_PROBE.md"
IN_P448 = GENERATED / "p448_current_strict_f447_t169_provider_provenance_audit_probe_summary.json"

OUT = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe.json"
OUT_SUMMARY = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe_summary.json"


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
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("-", " ")
        .replace("_", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P790, IN_F790, IN_F447, IN_N483, IN_T169, IN_P415, IN_P448]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P791",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p790 = load_json(IN_P790)
    f790 = load_json(IN_F790)
    p448 = load_json(IN_P448)

    f447_text = IN_F447.read_text(encoding="utf-8")
    n483_text = IN_N483.read_text(encoding="utf-8")
    t169_text = IN_T169.read_text(encoding="utf-8")
    p415_text = IN_P415.read_text(encoding="utf-8")

    target_required_fields = [
        item.get("name")
        for item in (((f790.get("target_object") or {}).get("required_fields")) or [])
        if isinstance(item, dict)
    ]

    template_domain_markers_present = (
        contains_all(f447_text, ["QW-2122", "T168", "vpsi", "g4", "g6"])
        and contains_all(n483_text, ["existence + uniqueness", "sign selection", "K_total"])
        and contains_all(t169_text, ["candidate family", "selection principle", "theorem-level existence + uniqueness"])
    )

    alpha_s_semantics_absent_from_template = not any(
        needle in normalize(f447_text + "\n" + n483_text + "\n" + t169_text)
        for needle in [
            "alpha_s",
            "n_f attachment",
            "n_f_active",
            "qcd",
            "normalized boundary",
            "canonical anchor",
            "boundary role",
        ]
    )

    p790_blocker_mentions_missing_selection = any(
        "selection principle" in normalize(item)
        for item in (p790.get("narrowed_anchor_blockers") or [])
    )

    checks = [
        {
            "id": "repo_exports_real_strict_selection_template_elsewhere",
            "pass": (
                p448.get("provider_classification") == "strict_derived"
                and p448.get("theorem_level_pass") is True
                and template_domain_markers_present
            ),
            "details": "F447/N483/T169 together show that the repo knows how to export a real strict selection/lift object with theorem-level existence and uniqueness.",
        },
        {
            "id": "existing_template_is_domain_specific_not_alpha_s_boundary",
            "pass": template_domain_markers_present and alpha_s_semantics_absent_from_template,
            "details": "The current strict selection template lives on the QW-2122 -> T168 per-site provider domain, not on the alpha_s normalized anchor-family domain.",
        },
        {
            "id": "existing_template_requires_non_alpha_s_input_chain",
            "pass": contains_all(f447_text, ["QW-2122", "R14", "R15", "element-order reference"])
            and contains_all(json.dumps(p790, ensure_ascii=True), ["f704_max_mode_anchor_family", "f789_interface_target"]),
            "details": "The F447/N483 template depends on QW-2122 scalar closure plus R14/R15 host-kernel inputs, whereas the alpha_s anchor lane currently sits on F704/P789/F790 normalized family objects.",
        },
        {
            "id": "p415_blocks_rhetorical_selector_transfer",
            "pass": contains_all(
                p415_text,
                [
                    "selection principles",
                    "not exported as strict typed objectives whose unique optimizers are proved",
                ],
            ),
            "details": "P415 explicitly blocks upgrading verbal selection rhetoric into a strict selector export.",
        },
        {
            "id": "selection_principle_ref_remains_missing_on_alpha_s_lane",
            "pass": p790_blocker_mentions_missing_selection and "selection_principle_ref" in target_required_fields,
            "details": "P790/F790 already classify selection_principle_ref as a still-missing field on the current alpha_s anchor lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    reusable_selection_principle_found = False

    if not failed_checks and not reusable_selection_principle_found:
        status = "P791_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
    else:
        status = "P791_REQUIRES_REVIEW"

    narrowed_selection_blockers = [
        "No current exported strict selection object acts directly on the P789 normalized alpha_s candidate-family domain.",
        "No current exported strict objective/order rule selects f704_max_mode_anchor_family over the other normalized alpha_s candidates.",
        "No current theorem-level uniqueness or finite-residual rule is exported for alpha_s family selection itself.",
    ]

    artifact = {
        "stage": "P791",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p790_canonical_anchor_probe": rel(IN_P790),
            "f790_anchor_selection_target": rel(IN_F790),
            "f447_selection_template_packet": rel(IN_F447),
            "n483_template_theorem": rel(IN_N483),
            "t169_template_target": rel(IN_T169),
            "p415_selector_rhetoric_boundary": rel(IN_P415),
            "p448_template_provenance_summary": rel(IN_P448),
        },
        "template_artifacts_considered": [
            {
                "id": "F447/N483/T169",
                "classification": "real strict selection-template family",
                "why_not_reusable": [
                    "Acts on the QW-2122 scalar-to-per-site provider domain rather than on alpha_s normalized candidate families.",
                    "Uses a different admissible input chain: QW-2122, R14, R15, typed Z12 element-order reference, and kernel-mixing energy.",
                    "Exports no alpha_s boundary-role semantics and no n_f attachment.",
                ],
            }
        ],
        "checks": checks,
        "failed_checks": failed_checks,
        "reusable_selection_principle_found": reusable_selection_principle_found,
        "current_honest_reading": [
            "The repo already knows the formal pattern of a strict selection theorem: declared domain, explicit selection rule, theorem-level existence and uniqueness, and exported value object.",
            "But the current exported template lives on a different semantic and input domain, so it cannot be silently reused as selection_principle_ref for alpha_s canonical anchor selection.",
            "Therefore the next honest move is not reuse by analogy, but freezing an alpha_s-specific selection-principle target object.",
        ],
        "narrowed_selection_blockers": narrowed_selection_blockers,
        "recommended_next_packet": {
            "id": "F791_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_TARGET_PACKET",
            "goal": "Freeze the exact alpha_s-specific selection-principle object still missing before selection_principle_ref can be filled.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "selection_objective_or_order_rule_ref",
                "strict_input_chain_ref",
                "uniqueness_or_finite_residual_rule_ref",
                "selected_family_output_schema",
                "nontransfer_boundary_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P791",
        "status": status,
        "as_of": AS_OF,
        "reusable_selection_principle_found": reusable_selection_principle_found,
        "selection_template_exists_elsewhere": True,
        "narrowed_selection_blocker_count": len(narrowed_selection_blockers),
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
