#!/usr/bin/env python3
"""P1203: classify no-CAS vs mathematical-negative outcomes for W4 workflow."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1202 = json.loads((GEN / "p1202_w4_cas_expand_bridge_probe_summary.json").read_text(encoding="utf-8"))

    sympy_available = bool(p1202.get("sympy_available", False))
    expand_ok = bool(p1202.get("expand_ok", False))

    if not sympy_available:
        classification = "NO_CAS_INFRA_BLOCK"
    elif sympy_available and not expand_ok:
        classification = "CAS_RUNTIME_FAILURE"
    else:
        classification = "CAS_STEP_OK"

    out = {
        "packet": "P1203",
        "as_of": "2026-05-10",
        "input_sympy_available": sympy_available,
        "input_expand_ok": expand_ok,
        "classification": classification,
        "is_mathematical_negative_result": False,
        "is_infrastructure_block": classification in {"NO_CAS_INFRA_BLOCK", "CAS_RUNTIME_FAILURE"},
        "strict_closure_claim_allowed": False,
        "note": "Separates infrastructure blockers from mathematical witness outcomes.",
    }

    out_path = GEN / "p1203_w4_no_cas_mode_classifier_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1203] classification={classification} wrote {out_path}")


if __name__ == "__main__":
    main()
