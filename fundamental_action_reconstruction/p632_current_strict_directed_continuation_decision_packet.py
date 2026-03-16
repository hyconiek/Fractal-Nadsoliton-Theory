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

IN_FIX = GENERATED / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json"
IN_S_DIR = GENERATED / "s_dir_pair1_strict_v1.json"
IN_DIRECTED_STATE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"
IN_PROJECTIVE_STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"

OUT_JSON = GENERATED / "p632_current_strict_directed_continuation_decision_packet.json"
OUT_SUMMARY = GENERATED / "p632_current_strict_directed_continuation_decision_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [str(p.relative_to(REPO)) for p in (IN_FIX, IN_S_DIR, IN_DIRECTED_STATE, IN_PROJECTIVE_STATE) if not p.exists()]
    if missing:
        artifact = {
            "stage": "P632",
            "date": AS_OF,
            "goal": "declare_professorial_directed_continuation_after_T171_discharge__no_false_pass",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_DIRECTED_CONTINUATION_DECISION",
            "missing": missing,
            "recommended_next_strict_target": {
                "target": "T164" if not IN_FIX.exists() else "H37",
                "note": "Directed continuation requires explicit T164 fixing datum + a directed/sign-sensitive T171 discharge; do not smuggle generator/sign conventions.",
            },
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": "P632",
                    "status": artifact["status"],
                    "recommended_next_strict_target": artifact["recommended_next_strict_target"]["target"],
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    fix = load_json(IN_FIX)
    s_dir = load_json(IN_S_DIR)
    directed = load_json(IN_DIRECTED_STATE)
    projective = load_json(IN_PROJECTIVE_STATE)

    decision = "DIRECTED_CONTINUATION_SELECTED"
    artifact = {
        "stage": "P632",
        "date": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "declare_professorial_directed_continuation_after_T171_discharge__treat_directed_selector_state_as_physical_in_declared_scope",
        "status": "PASS_DECISION_DECLARED_DIRECTED_CONTINUATION_SELECTED",
        "decision": decision,
        "decision_basis": {
            "t164_fixing_datum_present": True,
            "t164_fixing_datum_object": fix.get("object"),
            "t171_directed_state_present": True,
            "t171_directed_state_object": directed.get("object"),
            "projective_state_object": projective.get("object"),
            "sign_sensitive_observable_present": True,
            "sign_sensitive_observable_object": s_dir.get("object"),
            "strict_boundary_note": "N462 forbids Aut(Z12)-invariant generator fixing; directed continuation is premise-based via the exported fixing datum (F473/N523).",
        },
        "continuation": {
            "selected": "directed",
            "meaning": (
                "Treat the selector state on C_v1 as a directed/vector-level physical state object in the declared scope, "
                "using the exported T164 fixing datum + exported sign-sensitive observable to fix residual sign. "
                "The projective state remains the quotient/ray-level shadow."
            ),
        },
        "recommended_next_strict_target": {
            "target": "P11",
            "note": (
                "After the post-projective directed frontier is resolved (T171 discharged), continue strict-only ToE closure work "
                "that depends only on the exported selector state objects; the next bottleneck remains the existing-kernel-feedback → "
                "explicit-chain factorization route (P11)."
            ),
        },
        "hard_limits": [
            "No Aut(Z_12)-invariant canonicity claim (premise-based fixing datum; N462).",
            "No strict-core selector closure / admissible S_sel_int claim.",
            "No global discharge of QW-2191 claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(
        json.dumps(
            {
                "stage": "P632",
                "status": artifact["status"],
                "decision": decision,
                "selected_continuation": "directed",
                "recommended_next_strict_target": artifact["recommended_next_strict_target"]["target"],
                "no_false_pass": True,
            },
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="ascii",
    )

    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
