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
IN_ATLAS_GLOBAL = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
IN_TRANSITION_GLOBAL = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
IN_PROJECTOR_SECTION = GENERATED / "a_12345_pair12345_chart_glued_orientation_projector_operator_section_strict_core_v2.json"

OUT_STATE_GLOBAL = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f470_current_strict_global_projective_selector_state_object_export_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [
        str(p.relative_to(REPO))
        for p in (IN_CV1, IN_ATLAS_GLOBAL, IN_TRANSITION_GLOBAL, IN_PROJECTOR_SECTION)
        if not p.exists()
    ]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F470",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F306 exported strict configuration-space object C_v1",
                        "F469 exported global selector atlas and transition objects on C_v1 (T170 discharge)",
                        "F466 exported five-chart projector section with full cocycle data on {pair1..pair5}",
                    ],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    cv1 = load_json(IN_CV1)
    atlas_global = load_json(IN_ATLAS_GLOBAL)
    transition_global = load_json(IN_TRANSITION_GLOBAL)
    section = load_json(IN_PROJECTOR_SECTION)

    generated_utc = datetime.now(timezone.utc).isoformat()

    chart_ids = ["pair1", "pair2", "pair3", "pair4", "pair5"]
    atlas_charts = atlas_global.get("charts") or {}
    charts: dict[str, Any] = {}
    for cid in chart_ids:
        if isinstance(atlas_charts, dict) and cid in atlas_charts:
            charts[cid] = atlas_charts[cid]

    state = {
        "object": "SelectorState_global_C_v1_projective_strict_v1",
        "stage": "F470",
        "status": "actual_exported_global_projective_selector_state_object_on_C_v1_packaging_projector_section_and_global_atlas_transition__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit global projective/ray-level selector state object on the declared strict configuration space C_v1 "
            "by packaging the exported global selector atlas/transition objects (F469/N515) with the exported five-chart glued "
            "projector operator section on {pair1..pair5} (F466/N510). This is a projector/span (sign-gauge-safe) state object and "
            "does not claim sign-sensitive physical orientation, strict-core selector closure, or global QW-2191 discharge."
        ),
        "domain": {
            "configuration_space_object_ref": str(IN_CV1.relative_to(REPO)),
            "configuration_space_object": str(cv1.get("object") or "C_v1"),
            "charts_cover": str(atlas_global.get("domain", {}).get("charts_cover") or "U_pair1 = ... = U_pair5 = C_v1 (declared)"),
        },
        "state_type": {
            "level": "projective_ray_state",
            "encoding": "rank_one_projector_section_on_pair_charts",
            "sign_gauge": "projector/span semantics (u and -u identified at state level)",
        },
        "global_atlas_ref": str(IN_ATLAS_GLOBAL.relative_to(REPO)),
        "global_atlas_object": str(atlas_global.get("object") or ""),
        "global_transition_ref": str(IN_TRANSITION_GLOBAL.relative_to(REPO)),
        "global_transition_object": str(transition_global.get("object") or ""),
        "projector_section_ref": str(IN_PROJECTOR_SECTION.relative_to(REPO)),
        "projector_section_object": str(section.get("object") or ""),
        "charts": charts,
        "gluing": {
            "projector_section_level": True,
            "gluing_laws_ref": str(transition_global.get("gluing_laws_ref") or ""),
            "cocycle_discipline": transition_global.get("cocycle_discipline"),
        },
        "hard_limits": [
            "Projective/ray-level object only; no sign-sensitive directed orientation datum is derived.",
            "Does not promote section-level cocycle data to operator-level groupoid identities (see N512).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_STATE_GLOBAL.write_text(
        json.dumps(state, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )

    summary = {
        "stage": "F470",
        "status": "PASS_EXPORTED_GLOBAL_PROJECTIVE_SELECTOR_STATE_OBJECT",
        "exported": {"selector_state_global": str(OUT_STATE_GLOBAL.relative_to(REPO))},
        "inputs": {
            "selector_atlas_global": str(IN_ATLAS_GLOBAL.relative_to(REPO)),
            "selector_transition_global": str(IN_TRANSITION_GLOBAL.relative_to(REPO)),
            "projector_section": str(IN_PROJECTOR_SECTION.relative_to(REPO)),
        },
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )

    print(OUT_STATE_GLOBAL)


if __name__ == "__main__":
    main()

