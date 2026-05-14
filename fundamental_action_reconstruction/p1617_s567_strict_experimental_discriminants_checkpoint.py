#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1616 = GEN / "p1616_s566_strict_phenomenology_bounded_predictions_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    s16 = _load(IN1616)
    p = s16.get("predictions", {})

    discr = {
        "D1_scalar_response": {
            "observable": "effective scalar-mode response amplitude",
            "linked_prediction": "P1_scalar_mass_proxy",
            "required_relative_precision": 0.05,
            "systematics": ["background subtraction", "calibration drift"],
        },
        "D2_gr_coupling_shift": {
            "observable": "effective gravitational coupling shift",
            "linked_prediction": "P2_gr_coupling_proxy",
            "required_relative_precision": 0.03,
            "systematics": ["curvature modeling", "detector transfer function"],
        },
        "D3_curvature_backreaction": {
            "observable": "curvature-induced scalar backreaction trend",
            "linked_prediction": "P3_mix_backreaction_proxy",
            "required_relative_precision": 0.1,
            "systematics": ["environmental curvature proxy", "model truncation"],
        },
    }

    cl = s16.get("strict_core_closure", {})
    ready = (
        s16.get("status", "").startswith("PASS")
        and not cl.get("missing_exports")
        and not cl.get("missing_witnesses")
        and not cl.get("missing_theorems")
        and all(k in p for k in ["P1_scalar_mass_proxy", "P2_gr_coupling_proxy", "P3_mix_backreaction_proxy"])
    )

    summary = {
        "checkpoint": "P1617_S567",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1617_STRICT_EXPERIMENTAL_DISCRIMINANTS" if ready else "KEEP_OPEN_P1617_DISCRIMINANT_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "discriminants": discr,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": cl.get("missing_exports", []),
            "missing_witnesses": cl.get("missing_witnesses", []),
            "missing_theorems": cl.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Export strict data-assimilation protocol mapping measured discriminants back to kernel parameter posteriors.",
        "lay_summary": "To plan testów: mówimy dokładnie co mierzyć i z jaką precyzją, żeby sprawdzić przewidywania strict modelu.",
    }

    out = GEN / "p1617_s567_strict_experimental_discriminants_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
