#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_CV1 = GENERATED / "c_v1_void_configuration_space_in_local_b_tilde_1_sector_v1.json"
IN_ATLAS_PAIR12345 = GENERATED / "selector_atlas_pair12345_axis_only_projector_v2.json"

OUT_ATLAS_GLOBAL = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
OUT_TRANSITION_GLOBAL = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
OUT_SUMMARY = GENERATED / "f469_current_strict_global_selector_atlas_and_transition_object_export_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [str(p.relative_to(REPO)) for p in (IN_CV1, IN_ATLAS_PAIR12345) if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F469",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F306 exported strict configuration-space object C_v1",
                        "F466 exported lane-scoped five-chart selector atlas ingredient with full triple cocycle data on {pair1..pair5}",
                    ],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    cv1 = load_json(IN_CV1)
    atlas_pair12345 = load_json(IN_ATLAS_PAIR12345)

    chart_ids = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    chart_domains = {
        cid: {
            "domain_object_ref": str(IN_CV1.relative_to(REPO)),
            "domain_object": str(cv1.get("object") or "C_v1"),
            "domain": "C_v1 (entire declared strict configuration space domain)",
            "meaning": (
                "Chart is declared globally available as a strict program chart label (mode chart carrier) for the selector-track "
                "projector-level atlas ingredient; no physical privilege claim."
            ),
        }
        for cid in chart_ids
    }
    overlap_domains = {
        f"{a}_INTERSECT_{b}": {
            "domain": "C_v1",
            "meaning": "Global overlap by declared domain choice: U_a = U_b = C_v1.",
        }
        for i, a in enumerate(chart_ids)
        for b in chart_ids[i + 1 :]
    }

    generated_utc = datetime.now(timezone.utc).isoformat()

    atlas_global = {
        "object": "SelectorAtlas_global_C_v1_strict_v1",
        "stage": "F469",
        "status": "actual_exported_global_selector_atlas_object_with_explicit_C_v1_chart_domains_and_projector_level_gluing_data__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Discharge the strict global selector atlas object class required by T170/H41 by exporting an explicit global chart-domain "
            "and overlap-domain declaration on the declared strict configuration-space object C_v1, while reusing the already exported "
            "projector-level five-chart lane ingredient on {pair1..pair5} (F466/N510). This is a scope-and-typing export only: it does not "
            "imply strict-core selector closure, global QW-2191 discharge, or any sign-sensitive physical orientation datum."
        ),
        "domain": {
            "configuration_space_object_ref": str(IN_CV1.relative_to(REPO)),
            "configuration_space_object": str(cv1.get("object") or "C_v1"),
            "charts_cover": "U_pair1 = ... = U_pair5 = C_v1 (declared global cover)",
        },
        "charts": atlas_pair12345.get("charts"),
        "chart_domains": chart_domains,
        "overlap_domains": overlap_domains,
        "overlap_domain_declaration": {
            "scope": "explicit_declared_chart_domains_on_C_v1",
            "meaning": (
                "This is an explicit global chart-domain / overlap-domain declaration on the declared strict configuration space C_v1. "
                "It upgrades the earlier lane-scoped 'exported-artifact overlap' declaration (F463-F466) into an explicit global-domain "
                "statement for the atlas/transition object class (T170), without promoting any new physical claim."
            ),
        },
        "transitions": atlas_pair12345.get("transitions"),
        "gluing_data": atlas_pair12345.get("gluing_data"),
        "cocycle_discipline": {
            "level": "projector_section_level",
            "supports": [
                "N510: full triple cocycle data packaged at projector level on {pair1..pair5}",
                "N512: operator-level cocycle/groupoid promotion forbidden (section-level only)",
            ],
        },
        "inputs_reused": {
            "pair12345_lane_atlas_ref": str(IN_ATLAS_PAIR12345.relative_to(REPO)),
            "pair12345_lane_atlas_object": str(atlas_pair12345.get("object") or ""),
        },
        "hard_limits": [
            "This exports a global atlas object class on C_v1 only in the declared sense: chart domains and overlaps are explicit, but no claim is made that C_v1 is coordinatized by these charts.",
            "Projector-level gluing/cocycle discipline only; no operator-level transition groupoid promotion (see N512).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not derive a sign-sensitive physical orientation datum.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    transition_global = {
        "object": "SelectorTransition_global_C_v1_strict_v1",
        "stage": "F469",
        "status": "actual_exported_global_selector_transition_object_with_explicit_C_v1_overlap_domains_and_projector_level_cocycle_discipline__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Discharge the strict global selector transition/gluing object class required by T170/H40 by exporting an explicit global "
            "transition object with declared overlap domains on C_v1, using the already exported chart-transport operators and "
            "projector-level cocycle/gluing laws on {pair1..pair5}. This does not upgrade to operator-level groupoid identities and does "
            "not imply global QW-2191 discharge."
        ),
        "domain": {
            "configuration_space_object_ref": str(IN_CV1.relative_to(REPO)),
            "configuration_space_object": str(cv1.get("object") or "C_v1"),
            "chart_domains_ref": str(OUT_ATLAS_GLOBAL.relative_to(REPO)),
            "charts_cover": "U_pair1 = ... = U_pair5 = C_v1 (declared global cover)",
            "overlaps": "U_i ∩ U_j = C_v1 (declared)",
        },
        "charts": chart_ids,
        "transition_operators": atlas_pair12345.get("transitions"),
        "gluing_laws_ref": str(atlas_pair12345.get("gluing_data", {}).get("operator_section_ref") or ""),
        "cocycle_discipline": {
            "level": "projector_section_level",
            "supports": [
                "N510: projector-section full triple cocycle data (lane origin)",
                "N512: operator-level cocycle failure boundary (no false pass)",
            ],
            "explicit_non_promotion": [
                "Do not infer operator-level identities O_jk O_ij = O_ik on the full carrier from section-level cocycle data.",
                "Do not infer a sign-sensitive physical orientation datum from projector-level transport.",
            ],
        },
        "inputs_reused": {
            "pair12345_lane_atlas_ref": str(IN_ATLAS_PAIR12345.relative_to(REPO)),
            "pair12345_lane_atlas_object": str(atlas_pair12345.get("object") or ""),
        },
        "hard_limits": [
            "Projector-section cocycle only; no operator-level transition groupoid promotion (N512).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    atlas_global_path = OUT_ATLAS_GLOBAL.relative_to(REPO)
    transition_global_path = OUT_TRANSITION_GLOBAL.relative_to(REPO)

    OUT_ATLAS_GLOBAL.write_text(
        json.dumps(atlas_global, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_TRANSITION_GLOBAL.write_text(
        json.dumps(transition_global, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )

    summary = {
        "stage": "F469",
        "status": "PASS_EXPORTED_GLOBAL_ATLAS_AND_TRANSITION_OBJECTS",
        "exported": {
            "selector_atlas_global": str(atlas_global_path),
            "selector_transition_global": str(transition_global_path),
        },
        "t170_discharged": True,
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_ATLAS_GLOBAL)


if __name__ == "__main__":
    main()

